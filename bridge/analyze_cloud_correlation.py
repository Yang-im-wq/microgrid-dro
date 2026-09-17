# -*- coding: utf-8 -*-
"""
analyze_cloud_correlation.py —— 从卫星云图序列估计云的时间自相关
================================================================
输入：上游视觉模块输出的云量时间序列
      D:\\pv-probabilistic-forecast\\vision\\output\\cloud_cover_yolo.csv
输出：D:\\microgrid-dro\\data\\cloud_correlation.json

为什么要做这件事
----------------
本项目的 `gen_pv_scenarios.m` 生成光伏场景时是**逐小时独立抽样**的，
等于假设"这一小时的光伏和下一小时完全无关"。

但真实云不是这样的：卫星实测显示相邻小时云量只变化 0.037（中位数 0.011），
滞后 1 小时的自相关高达 0.95。光伏出力直接受云影响，所以也应当是强时间相关的。

把云量序列的实测自相关提取出来喂给场景生成器，就是把"云图"这个数据源
真正接进了优化链路 —— 不是为了好看，而是补上一个会让结论失真的假设。

做法
----
对云量序列拟合一个 AR(1) 过程 c(t) = rho * c(t-1) + eps，
用滞后 1..K 的自相关做最小二乘拟合 rho：
    rho_k(理论) = rho^k      =>   ln(rho_k) = k * ln(rho)
取 K 个滞后做线性回归估计 ln(rho)，再取指数。
这样比只看滞后 1 稳健（滞后 1 的单点估计对噪声敏感）。
"""
import json
import os
import sys

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

CLOUD_CSV = os.environ.get(
    "CLOUD_SERIES",
    r"D:\pv-probabilistic-forecast\vision\output\cloud_cover_yolo.csv")
OUT_PATH = os.path.join(
    os.environ.get("MICROGRID_DRO_DATA", r"D:\microgrid-dro\data"),
    "cloud_correlation.json")

MAX_LAG = 8          # 用滞后 1..8 小时拟合


def spearman_autocorr(x, max_lag=MAX_LAG):
    """秩自相关（Spearman）。用秩相关而不是 Pearson，是因为场景生成器用的是
    高斯 Copula —— Copula 保持的是**秩相关**，不是 Pearson 相关。"""
    s = pd.Series(x).rank().to_numpy()
    s = s[~np.isnan(s)]
    if s.std() < 1e-12:
        return []
    return [float(np.corrcoef(s[:-k], s[k:])[0, 1]) for k in range(1, max_lag + 1)]


def kendall_lag1(x):
    """滞后 1 的 Kendall tau。选它而不是 Pearson/Spearman：
       - Pearson 会被云量的**非线性**扭曲
       - Spearman 在 CLM 这种量化严重（大量并列）的数据上会退化（实测滞后 7-8 直接 NaN）
       - Kendall tau 对并列值稳健，且高斯 Copula 有闭式换算 tau = (2/pi)*arcsin(rho)"""
    try:
        from scipy.stats import kendalltau
    except ImportError:
        return np.nan
    x = np.asarray(x, dtype=float)
    x = x[~np.isnan(x)]
    if len(x) < 20:
        return np.nan
    return float(kendalltau(x[:-1], x[1:]).correlation)


def tau_to_latent(tau):
    """高斯 Copula：tau = (2/pi)*arcsin(rho)  =>  rho = sin(pi*tau/2)"""
    return float(np.sin(np.pi * float(np.clip(tau, -0.999, 0.999)) / 2.0))


def fit_ar1(x, max_lag=MAX_LAG):
    """最小二乘拟合 AR(1) 的 rho；返回 (rho, 各滞后自相关)"""
    x = np.asarray(x, dtype=float)
    x = x[~np.isnan(x)]
    if x.std() < 1e-12:
        return np.nan, []
    ac = [float(np.corrcoef(x[:-k], x[k:])[0, 1]) for k in range(1, max_lag + 1)]
    # 秩相关在并列值很多时可能给出 NaN（云量量化严重），统一兜住
    ac = np.nan_to_num(np.array(ac), nan=1e-6, posinf=1.0, neginf=1e-6)
    ac = np.clip(ac, 1e-6, 1.0)
    k = np.arange(1, max_lag + 1)
    slope = np.polyfit(k, np.log(ac), 1)[0]      # ln(rho_k) = k*ln(rho)
    return float(np.exp(slope)), ac.tolist()


def main():
    if not os.path.exists(CLOUD_CSV):
        print("找不到云量序列: %s" % CLOUD_CSV)
        print("（这是可选输入；缺失时 gen_pv_scenarios 会用内置默认值）")
        return 1

    d = pd.read_csv(CLOUD_CSV)
    d["time"] = pd.to_datetime(d["time"])

    print("=" * 66)
    print("从卫星云图序列估计时间自相关")
    print("=" * 66)
    print("数据: %s" % CLOUD_CSV)
    print("       %d 行, %s ~ %s" % (len(d), d.time.min(), d.time.max()))

    # 载入光伏分位函数（标定潜变量相关要用到真实边缘分布）
    pv_csv = os.path.join(os.environ.get("MICROGRID_DRO_DATA", r"D:\microgrid-dro\data"),
                          "pv_daily_profile_hourly.csv")
    if not os.path.exists(pv_csv):
        print("找不到光伏曲线 %s，无法标定潜变量相关" % pv_csv)
        return 1
    pv = pd.read_csv(pv_csv)
    qcols = ["q%02d" % q for q in [5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95]]
    Q = pv[qcols].to_numpy()                       # nt x 11
    qlev = np.array([0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95])

    res = {}
    for col, label in [("cloud_fraction_true", "CLM 官方云掩膜(权威)"),
                       ("cloud_fraction", "YOLO 检测结果")]:
        if col not in d.columns:
            continue
        rho_p, ac_p = fit_ar1(d[col].values)        # 实测 Pearson 自相关 -> AR(1)
        tau = kendall_lag1(d[col].values)           # 滞后 1 的 Kendall tau
        latent = tau_to_latent(tau)                 # -> 高斯 Copula 潜变量相关
        res[col] = {"ar1_rho_pearson": rho_p, "autocorr": ac_p,
                    "kendall_tau1": tau, "latent_rho": latent, "label": label}
        print("\n%-22s (%s)" % (col, label))
        print("  实测 Pearson 自相关(滞后1..%d): %s" % (MAX_LAG, " ".join("%.3f" % v for v in ac_p)))
        print("  拟合 AR(1) rho(观测)          = %.4f" % rho_p)
        print("  滞后1 Kendall tau             = %.4f" % tau)
        print("  -> 高斯 Copula 潜变量相关     = %.4f" % latent)

    # 用权威的 CLM 结果作为推荐值
    key = "cloud_fraction_true" if "cloud_fraction_true" in res else "cloud_fraction"
    rho_rec = res[key]["latent_rho"]

    # 相邻小时跳变幅度（作为旁证）
    for col in res:
        dd = d[col].diff().abs().dropna()
        res[col]["hourly_jump_mean"] = float(dd.mean())
        res[col]["hourly_jump_median"] = float(dd.median())

    out = {
        "source": CLOUD_CSV,
        "n_samples": int(len(d)),
        "time_range": [str(d.time.min()), str(d.time.max())],
        "sample_interval_hours": 1,
        "fit_method": ("AR(1) 最小二乘拟合 ln(ac_k) = k*ln(rho)，滞后 1..%d（用于观测自相关）；"
                       "潜变量相关由滞后1 Kendall tau 换算 rho = sin(pi*tau/2)" % MAX_LAG),
        "recommended_rho_time": round(rho_rec, 4),
        "recommended_basis": key,
        "parameter_meaning": ("rho_time 是**高斯 Copula 的潜变量相关**，不是直接观测量。"
                              "云量本身（观测）的滞后1 Kendall tau 与它满足 "
                              "tau = (2/pi)*arcsin(rho)。"
                              "选 Kendall tau 而非 Pearson/Spearman 的原因："
                              "Pearson 会被云量的非线性扭曲；Spearman 在 CLM 这种量化严重、"
                              "大量并列的数据上会退化（实测滞后 7-8 直接算不出）。"),
        "per_series": {
            k: {"ar1_rho_observed": round(v["ar1_rho_pearson"], 4),
                "kendall_tau1": round(v["kendall_tau1"], 4),
                "latent_rho": round(v["latent_rho"], 4),
                "autocorr_observed_lag1_8": [round(a, 4) for a in v["autocorr"]],
                "hourly_jump_mean": round(v["hourly_jump_mean"], 4),
                "label": v["label"]}
            for k, v in res.items()},
        "note": ("光伏场景生成器 gen_pv_scenarios.m 默认用 recommended_rho_time "
                 "作为时间相关系数；这是从卫星云图实测得到的，不是拍脑袋的常数。"),
    }

    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    with open(OUT_PATH, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)

    print("\n" + "=" * 66)
    print("实测（观测）时间自相关 rho_obs = %.4f" % res[key]["ar1_rho_pearson"])
    print("标定出的潜变量相关 rho_time  = %.4f  （取自 %s）" % (rho_rec, key))
    print("已写出: %s" % OUT_PATH)
    print("=" * 66)
    return 0


if __name__ == "__main__":
    sys.exit(main())
