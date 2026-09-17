# -*- coding: utf-8 -*-
"""
plot_profile.py —— 数据桥成果可视化（P0 的"看效果"交付物）
================================================================
画两张图：
  左：96 点标幺日内曲线 + 不确定性带（这是喂给 MATPOWER 的输入）
  右：测试集里抽几个真实日，预测分位数带 vs 真实出力（验证形状没跑偏）

用法：D:\\anaconda\\python.exe D:\\microgrid-dro\\bridge\\plot_profile.py
"""
import os
import sys
import json

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

DATA_DIR = r"D:\microgrid-dro\data"
FIG_DIR = r"D:\microgrid-dro\results\figures"
QUANTILES = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
Q_LABELS = ["q{:02.0f}".format(q * 100) for q in QUANTILES]
I05, I10, I25, I50, I75, I90, I95 = 0, 1, 2, 5, 8, 9, 10


def main():
    os.makedirs(FIG_DIR, exist_ok=True)

    T96 = pd.read_csv(os.path.join(DATA_DIR, "pv_daily_profile.csv"))
    Tq = pd.read_csv(os.path.join(DATA_DIR, "pv_quantile.csv"))
    with open(os.path.join(DATA_DIR, "pv_profile_meta.json"), encoding="utf-8") as f:
        meta = json.load(f)
    scale_kw = float(meta["peak_reference_kw"])      # 归一化基准 = 装机容量 (kW)
    M = T96[Q_LABELS].to_numpy()
    hours = T96["hour"].to_numpy() + T96["minute"].to_numpy() / 60.0

    fig, axes = plt.subplots(1, 2, figsize=(15, 5.5))

    # ---------------- 左：标幺日内曲线 + 不确定性带 ----------------
    ax = axes[0]
    ax.fill_between(hours, M[:, I05], M[:, I95], alpha=0.18, color="#1f77b4",
                    label="q05–q95 带")
    ax.fill_between(hours, M[:, I25], M[:, I75], alpha=0.30, color="#1f77b4",
                    label="q20–q80 带")
    ax.plot(hours, M[:, I50], "-", color="#1f77b4", lw=2.5, label="中位数 q50")
    ax.axhline(1.0, color="crimson", ls="--", lw=1.2, label="装机容量 (1.0 p.u.)")
    ax.set_xlim(0, 24)
    ax.set_xticks(range(0, 25, 2))
    ax.set_ylim(0, 1.1)
    ax.set_xlabel("时刻 (h)")
    ax.set_ylabel("归一化出力 (p.u.)")
    ax.set_title("典型日标幺曲线（喂给 MATPOWER 的输入）\n1.0 = 装机容量", fontsize=12)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9, loc="upper left")

    # 标注峰值
    i_pk = int(np.argmax(M[:, I50]))
    ax.annotate("q50 峰值 {:.4f} @ {:.0f}:{:02.0f}".format(M[i_pk, I50], np.floor(hours[i_pk]),
                (hours[i_pk] % 1) * 60),
                xy=(hours[i_pk], M[i_pk, I50]), xytext=(14.5, 0.55),
                arrowprops=dict(arrowstyle="->", color="gray"), fontsize=9, color="gray")

    # ---------------- 右：真实日对照 ----------------
    ax = axes[1]
    t = pd.to_datetime(Tq["time"])
    days = pd.Series(t.dt.normalize()).unique()
    # 挑连续 3 天，且这 3 天点数足够
    pick = [d for d in days[40:60]][:3]
    colors = ["#2ca02c", "#ff7f0e", "#9467bd"]
    for d, c in zip(pick, colors):
        m = (t.dt.normalize() == d).to_numpy()
        if m.sum() < 20:
            continue
        hh = t[m].dt.hour.to_numpy() + t[m].dt.minute.to_numpy() / 60.0
        sub = Tq[m]
        # 真实值换算到标幺（同一分母 = 装机容量）
        ax.fill_between(hh, sub[Q_LABELS[I05]] / scale_kw, sub[Q_LABELS[I95]] / scale_kw,
                        color=c, alpha=0.12)
        ax.plot(hh, sub[Q_LABELS[I50]] / scale_kw, "-", color=c, lw=1.8,
                label="{} 预测q50".format(str(d)[:10]))
        ax.plot(hh, sub["y_true"] / scale_kw, "o", color=c, ms=2.5, alpha=0.6,
                label="{} 真实值".format(str(d)[:10]))
    ax.axhline(1.0, color="crimson", ls="--", lw=1.2)
    ax.set_xlim(6, 20)
    ax.set_ylim(0, 1.1)
    ax.set_xlabel("时刻 (h)")
    ax.set_ylabel("归一化出力 (p.u.)")
    ax.set_title("抽样真实日：预测分位数带 vs 真实出力", fontsize=12)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, ncol=2, loc="upper left")

    fig.suptitle("P0 数据桥输出：光伏概率预测 → 微电网优化输入", fontsize=14)
    fig.tight_layout()
    out = os.path.join(FIG_DIR, "p0_pv_profile.png")
    fig.savefig(out, dpi=140)
    print("已保存: {}".format(out))

    # ---------------- 数值摘要 ----------------
    print("\n=== 标幺日内曲线摘要 ===")
    print("  q50 峰值      = {:.4f} p.u. @ {:02.0f}:{:02.0f}".format(
        M[i_pk, I50], np.floor(hours[i_pk]), (hours[i_pk] % 1) * 60))
    print("  q95 峰值      = {:.4f} p.u.".format(M[:, I95].max()))
    print("  q05 峰值      = {:.4f} p.u.".format(M[:, I05].max()))
    print("  正午不确定性带宽 (q95-q05) @12:00 = {:.4f} p.u.".format(
        M[48, I95] - M[48, I05]))
    h_lo, h_hi = hours[M[:, I50] > 0].min(), hours[M[:, I50] > 0].max()
    print("  日出/日落边界  = {:02.0f}:{:02.0f} / {:02.0f}:{:02.0f}".format(
        np.floor(h_lo), round((h_lo % 1) * 60), np.floor(h_hi), round((h_hi % 1) * 60)))


if __name__ == "__main__":
    main()
