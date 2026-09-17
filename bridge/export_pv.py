# -*- coding: utf-8 -*-
"""
export_pv.py —— 数据桥：光伏概率预测 → MATPOWER 微电网算例
================================================================
项目：微电网分布鲁棒优化预研平台（D:\microgrid-dro）

职责
----
1. 重跑 LSTM 分位数推理，拿全 **11 个分位数**
   （`results/predictions_quantile.npz` 只存了 q05/q50/q95 三个，
     11 分位数的张量只存在于模型输出内存，evaluate.py 切片后即丢）
2. 修正分位数交叉 + 截断负值
3. 导出全时段分位数 → data/pv_quantile.csv
4. 聚合成 96 点标幺日内曲线（夜间补零）→ data/pv_daily_profile.csv
   以及 24 点小时版（给 MOST 的 nt=24）→ data/pv_daily_profile_hourly.csv

量级说明
--------
原始数据是 **50MW 电站**，而 IEEE 33 节点总负荷只有 **3.715MW**（差 13 倍）。
所以本脚本**只导出归一化形状**（q50 日内峰值 = 1.0），绝对容量由 MATLAB 侧按
渗透率决定。微电网研究看的是渗透率，不是绝对 MW。

用法（必须用 pv 项目的 Python 环境，因为要 import 已有模型）
--------------------------------------------
    D:\\anaconda\\python.exe D:\\microgrid-dro\\bridge\\export_pv.py

本脚本**不修改** pv 项目里的任何文件，只读取。
"""
import os
import sys
import json

import numpy as np
import pandas as pd

# ---- 控制台按 UTF-8 输出，避免 Windows GBK 把中文打乱 ----
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

# ============================================================
# 路径
# ============================================================
# 上游光伏预测项目（只读）。默认是本机路径，换机器时用环境变量覆盖：
#     set PV_PROJECT_ROOT=D:\somewhere\pv-probabilistic-forecast
PV_ROOT = os.environ.get("PV_PROJECT_ROOT", r"D:\pv-probabilistic-forecast")
OUT_DIR = os.environ.get("MICROGRID_DRO_DATA", r"D:\microgrid-dro\data")

QUANTILES = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
Q_LABELS = ["q{:02.0f}".format(q * 100) for q in QUANTILES]   # ['q05', ..., 'q95']
SLOTS_PER_HOUR = 4          # 15 分钟 → 每小时 4 个时段
SLOTS_PER_DAY = 24 * SLOTS_PER_HOUR


def banner(msg):
    print("=" * 68)
    print(msg)
    print("=" * 68)


# ============================================================
# 1. 推理：拿全 11 个分位数
# ============================================================
def infer_quantiles():
    """重跑推理，返回 (time_str, y_true_kw, pred_kw) —— pred_kw 形状 (n, 11)。

    模板取自 pv 项目 dro_demo.py:33-52（已正确处理逆标准化 *scale + mean）。
    """
    # 让 pv 项目的模块可被 import。config.py 内部用 __file__ 定位，
    # 所以不依赖 cwd，从任何地方跑都行。
    if PV_ROOT not in sys.path:
        sys.path.insert(0, PV_ROOT)

    import torch
    import config
    from data import preprocess
    from models.lstm_quantile import LSTMQuantileForecaster

    data = preprocess.get_dataset()
    X_test = torch.from_numpy(data["X_test"]).float()
    y_test_scaled = data["y_test"]                  # (n, horizon)
    mean, scale = data["target_mean"], data["target_scale"]

    model = LSTMQuantileForecaster(
        X_test.shape[-1], config.HIDDEN_SIZE, config.NUM_LAYERS,
        config.DROPOUT, config.HORIZON, config.NUM_QUANTILES, config.MODEL_TYPE)
    ckpt = os.path.join(config.CHECKPOINT_DIR, "best_quantile.pt")
    model.load_state_dict(torch.load(ckpt, map_location="cpu"))
    model.eval()

    with torch.no_grad():
        pred_scaled = model(X_test).numpy()          # (n, horizon, nq)

    # 反归一化回 kW，取第一个预测步（horizon = 1）
    y_true = y_test_scaled[:, 0] * scale + mean
    pred = pred_scaled[:, 0, :] * scale + mean       # (n, 11)

    assert pred.shape[1] == len(QUANTILES), \
        "模型输出 {} 个分位数，但期望 {}".format(pred.shape[1], len(QUANTILES))

    return data["t_test"], y_true, pred, (mean, scale)


# ============================================================
# 2. 清洗：修分位数交叉 + 截断负值
# ============================================================
def clean_quantiles(pred):
    """返回 (pred_clean, stats)。顺序：先修交叉，再截断负值。"""
    pred = np.asarray(pred, dtype=float)
    n = pred.shape[0]

    # ---- 修分位数交叉：Q(τ) 必须随 τ 单调不减 ----
    # 用 maximum.accumulate（与 calibrate.py:62 同一手法，保持项目一致性）
    before = pred.copy()
    pred = np.maximum.accumulate(pred, axis=1)
    crossed_rows = int((pred != before).any(axis=1).sum())
    max_fix = float(np.abs(pred - before).max()) if crossed_rows else 0.0

    # ---- 截断负值：光伏出力不可能为负 ----
    neg_counts = {Q_LABELS[j]: int((pred[:, j] < 0).sum()) for j in range(pred.shape[1])}
    n_neg_total = int((pred < 0).any(axis=1).sum())
    pred = np.clip(pred, 0.0, None)

    stats = {
        "n_samples": n,
        "quantile_crossing_rows_fixed": crossed_rows,
        "quantile_crossing_max_shift_kw": round(max_fix, 4),
        "rows_with_negative_before_clip": n_neg_total,
        "negative_count_per_quantile": neg_counts,
    }
    return pred, stats


# ============================================================
# 3. 导出全时段分位数
# ============================================================
def export_timeseries(t_str, y_true, pred, out_path):
    df = pd.DataFrame({"time": np.asarray(t_str, dtype=str)})
    for j, lab in enumerate(Q_LABELS):
        df[lab] = pred[:, j]
    df["y_true"] = y_true
    df.to_csv(out_path, index=False, encoding="utf-8")
    return df


# ============================================================
# 4. 聚合成日内标幺曲线
# ============================================================
def build_daily_profile(t_str, pred, peak_kw):
    """按「日内时段」聚合出典型日曲线，夜间补零，再按 peak_kw 归一化。

    peak_kw 必须传**原始逐点序列的全局最大值**（不是日均曲线峰值），理由见下。

    返回 (profile_96, profile_24, slot_counts, meta)
      profile_96 : (96, 11) —— 15 分钟分辨率，标幺值（1.0 = 装机容量）
      profile_24 : (24, 11) —— 小时分辨率
    """
    t = pd.to_datetime(pd.Series(np.asarray(t_str, dtype=str)))
    slot = (t.dt.hour * SLOTS_PER_HOUR + t.dt.minute // 15).to_numpy()   # 0..95

    n_q = pred.shape[1]
    profile_96 = np.zeros((SLOTS_PER_DAY, n_q), dtype=float)
    slot_counts = np.zeros(SLOTS_PER_DAY, dtype=int)

    for s in range(SLOTS_PER_DAY):
        m = (slot == s)
        slot_counts[s] = int(m.sum())
        if m.any():
            profile_96[s, :] = pred[m, :].mean(axis=0)
        # else: 该时段无数据 → 保持 0（夜间，物理正确）

    # ---- 归一化：所有分位数共用同一个分母，否则不确定性带的相对形状会被扭曲 ----
    # 分母取「原始逐点序列的全局最大值」，而不是日均曲线峰值。
    #
    # 为什么：1.0 必须是物理硬上限 —— 1.0 = 装机容量，光伏出力不可能超过装机容量。
    # 若按日均峰值归一，个别晴天（强于平均日）会冲到 1.31 p.u.，且 21% 的逐点值
    # 越界，MATLAB 侧就没法把 1.0 当容量上限用。
    # 代价：日均曲线峰值只有约 0.76 p.u.（典型日出力达装机的 76%）。
    # 这符合实际 —— 光伏电站很少在典型日打到铭牌。
    i_med = QUANTILES.index(0.5)
    peak = float(peak_kw)
    if peak <= 0:
        raise RuntimeError("归一化基准为 0")
    profile_96_n = profile_96 / peak

    # ---- 小时版：对每个小时内的 4 个时段取平均 ----
    profile_24_n = profile_96_n.reshape(24, SLOTS_PER_HOUR, n_q).mean(axis=1)

    meta = {
        "normalization": "原始逐点全局最大值归一到 1.0 → 1.0 = 装机容量(物理硬上限)",
        "peak_reference_kw": round(peak, 2),
        "daily_profile_peak_pu": round(float(profile_96_n.max()), 6),
        "q50_daily_peak_pu": round(float(profile_96_n[:, i_med].max()), 6),
        "note_daily_peak": "日均曲线峰值 < 1.0 属正常：典型日出力低于装机容量",
        "slots_per_day": SLOTS_PER_DAY,
        "hours_covered_by_raw_data": sorted(set((np.asarray(slot) // SLOTS_PER_HOUR).tolist())),
        "hours_zero_filled": sorted(set(range(24)) - set((np.asarray(slot) // SLOTS_PER_HOUR).tolist())),
        "n_days": int(t.dt.normalize().nunique()),
        "date_range": [str(t.min()), str(t.max())],
    }
    return profile_96_n, profile_24_n, slot_counts, meta


def export_profile(profile, counts, out_path, resolution_min):
    """counts 必须已聚合到 profile 自身的分辨率（96 点给 96 个，24 点给 24 个）。"""
    n = profile.shape[0]
    rows = []
    for s in range(n):
        total_min = s * resolution_min
        rows.append([s, total_min // 60, total_min % 60, int(round(counts[s]))]
                    + list(profile[s, :]))
    cols = ["slot", "hour", "minute", "n_samples"] + Q_LABELS
    df = pd.DataFrame(rows, columns=cols)
    # 不用 utf-8-sig：这些 CSV 全 ASCII，BOM 只会让 MATLAB readtable 把首列读成
    # '﻿time'。utf-8 无 BOM 对 MATLAB 和 Excel 都安全。
    df.to_csv(out_path, index=False, encoding="utf-8")
    return df


# ============================================================
# main
# ============================================================
def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    banner("STEP 1/4  重跑推理，拿全 11 个分位数")
    t_str, y_true, pred_raw, (mean, scale) = infer_quantiles()
    print("  样本数        : {}".format(len(t_str)))
    print("  分位数        : {}".format(QUANTILES))
    print("  输出张量      : {}".format(pred_raw.shape))
    print("  归一化参数    : mean={:.2f} kW  scale={:.2f} kW".format(mean, scale))
    print("  时间范围      : {} ~ {}".format(t_str[0], t_str[-1]))

    banner("STEP 2/4  清洗（修分位数交叉 + 截断负值）")
    pred, stats = clean_quantiles(pred_raw)
    print("  分位数交叉修正行数 : {} / {}  (最大调整 {:.4f} kW)".format(
        stats["quantile_crossing_rows_fixed"], stats["n_samples"],
        stats["quantile_crossing_max_shift_kw"]))
    print("  含负值行数(截断前) : {}".format(stats["rows_with_negative_before_clip"]))
    for lab, c in stats["negative_count_per_quantile"].items():
        if c:
            print("    {} : {} 个负值".format(lab, c))
    print("  截断后最小出力     : {:.2f} kW".format(pred.min()))

    banner("STEP 3/4  导出全时段分位数")
    ts_path = os.path.join(OUT_DIR, "pv_quantile.csv")
    df_ts = export_timeseries(t_str, y_true, pred, ts_path)
    print("  已写出: {}  ({} 行 × {} 列)".format(ts_path, len(df_ts), df_ts.shape[1]))

    banner("STEP 4/4  聚合成日内标幺曲线")
    peak_kw = float(pred.max())      # 原始逐点全局最大值 = 归一化基准 = 装机容量
    p96, p24, slot_counts, meta = build_daily_profile(t_str, pred, peak_kw)
    p96_path = os.path.join(OUT_DIR, "pv_daily_profile.csv")
    p24_path = os.path.join(OUT_DIR, "pv_daily_profile_hourly.csv")
    df96 = export_profile(p96, slot_counts.astype(float), p96_path, resolution_min=15)
    counts24 = slot_counts.reshape(24, SLOTS_PER_HOUR).mean(axis=1)
    df24 = export_profile(p24, counts24, p24_path, resolution_min=60)
    print("  已写出: {}  ({} 行, 15 分钟分辨率)".format(p96_path, len(df96)))
    print("  已写出: {}  ({} 行, 小时分辨率)".format(p24_path, len(df24)))
    print("  归一化基准: 原始逐点最大值 = {:.1f} kW → 1.0  (= 装机容量)".format(meta["peak_reference_kw"]))
    print("              日均曲线峰值 = {:.4f} p.u.  (典型日达装机 {:.1f}%)".format(
        meta["daily_profile_peak_pu"], 100 * meta["daily_profile_peak_pu"]))
    print("              其中 q50 日均峰值 = {:.4f} p.u.".format(meta["q50_daily_peak_pu"]))
    print("  原始数据覆盖小时: {}".format(meta["hours_covered_by_raw_data"]))
    print("  补零小时(夜间)  : {}".format(meta["hours_zero_filled"]))
    print("  覆盖天数        : {}".format(meta["n_days"]))

    # ---- 元数据存盘，给 MATLAB 侧读 ----
    meta_out = dict(meta)
    meta_out.update(stats)
    meta_out["quantiles"] = QUANTILES
    meta_out["q_labels"] = Q_LABELS
    meta_out["source_project"] = "companion PV probabilistic-forecasting project (LSTM quantile regression)"
    meta_out["note_scale"] = ("原始为 50MW 电站数据，仅导出归一化形状；"
                              "绝对容量由 MATLAB 侧按渗透率决定")
    meta_path = os.path.join(OUT_DIR, "pv_profile_meta.json")
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(meta_out, f, ensure_ascii=False, indent=2)
    print("  已写出: {}".format(meta_path))

    # ---- 合理性自检 ----
    banner("自检")
    checks = []
    i_lo, i_hi, i_med = 0, len(QUANTILES) - 1, QUANTILES.index(0.5)
    mono_ok = bool((np.diff(pred, axis=1) >= -1e-9).all())
    checks.append(("分位数全时段单调", mono_ok))
    checks.append(("无负值", bool(pred.min() >= 0)))
    raw_n = pred / peak_kw          # 原始逐点序列归一化
    checks.append(("原始逐点峰值 = 1.0 (= 装机容量)", bool(abs(raw_n.max() - 1.0) < 1e-9)))
    checks.append(("原始逐点全 ≤ 1.0 (不超装机)", bool(raw_n.max() <= 1.0 + 1e-9)))
    checks.append(("夜间时段全零", bool(np.allclose(p96[:7 * SLOTS_PER_HOUR, :], 0))))
    checks.append(("上界 ≥ 中位数", bool((pred[:, i_hi] >= pred[:, i_med] - 1e-9).all())))
    checks.append(("中位数 ≥ 下界", bool((pred[:, i_med] >= pred[:, i_lo] - 1e-9).all())))
    checks.append(("标幺曲线 ≤ 1.0 (不超装机)", bool(p96.max() <= 1.0 + 1e-9)))
    checks.append(("日均峰值 < 1.0 (典型日低于装机)", bool(p96.max() < 1.0)))
    for name, ok in checks:
        print("  [{}] {}".format("PASS" if ok else "FAIL", name))
    n_fail = sum(1 for _, ok in checks if not ok)

    print("\n" + "=" * 68)
    print("完成。" + ("全部自检通过。" if n_fail == 0 else "有 {} 项自检未通过！".format(n_fail)))
    print("=" * 68)
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
