# -*- coding: utf-8 -*-
"""
make_load_profile.py —— 由 IEEE RTS-79 标准数据生成日负荷曲线
================================================================
输出（写到 D:\\microgrid-dro\\data\\）：
  load_profile.csv           当前生效的曲线（默认 = 冬季工作日）
  load_profile.variants.csv  全部 6 种变体并排，供切换

数据来源
--------
IEEE Committee Report, "IEEE Reliability Test System",
IEEE Transactions on Power Apparatus and Systems, Vol. PAS-98, No. 6,
pp. 2047-2054, Nov/Dec 1979.  （即 IEEE RTS-79）

本项目用的是其中 Table 3：**每小时峰值占当日峰值的百分比**。
对应原文的表述是 "Table 3: Hourly peak load in percent of daily peak"。

为什么选它
----------
这是配电网/电力系统研究里最常被引用的日负荷曲线之一，**有明确出处、可复核**。
上游光伏数据集只有光伏与气象列、没有负荷数据，所以负荷形状一直只能靠合成。
换成这条曲线后，至少"形状"这一环有了可引用的依据。

用法
----
    python make_load_profile.py                 # 写出冬季工作日版
    python make_load_profile.py --variant summer_wd
    python make_load_profile.py --list          # 列出所有变体
"""
import argparse
import csv
import os

OUT_DIR = os.environ.get("MICROGRID_DRO_DATA", r"D:\microgrid-dro\data")

# ------------------------------------------------------------------
# IEEE RTS-79 Table 3：每小时峰值占当日峰值的百分比
# 小时顺序为 1..24（对应原文的 12-1am ... 11-12pm）
# ------------------------------------------------------------------
RTS79_TABLE3 = {
    #             h1  h2  h3  h4  h5  h6  h7  h8  h9 h10 h11 h12
    #            h13 h14 h15 h16 h17 h18 h19 h20 h21 h22 h23 h24
    "winter_wd":   [67, 63, 60, 59, 59, 60, 74, 86, 95, 96, 96, 95,
                    95, 95, 93, 94, 99, 100, 100, 96, 91, 83, 73, 63],
    "winter_we":   [78, 72, 68, 66, 64, 65, 66, 70, 80, 88, 90, 91,
                    90, 88, 87, 87, 91, 100, 99, 97, 94, 92, 87, 81],
    "summer_wd":   [64, 60, 58, 56, 56, 58, 64, 76, 87, 95, 99, 100,
                    99, 100, 100, 97, 96, 96, 93, 92, 92, 93, 87, 72],
    "summer_we":   [74, 70, 66, 65, 64, 62, 62, 66, 81, 86, 91, 93,
                    93, 92, 91, 91, 92, 94, 95, 95, 100, 93, 88, 80],
    "springfall_wd": [63, 62, 60, 58, 59, 65, 72, 85, 95, 99, 100, 99,
                      93, 92, 90, 88, 90, 92, 96, 98, 96, 90, 80, 70],
    "springfall_we": [75, 73, 69, 66, 65, 65, 68, 74, 83, 89, 92, 94,
                      91, 90, 90, 86, 85, 88, 92, 100, 97, 95, 90, 85],
}

DESC = {
    "winter_wd": "冬季工作日（晚高峰，光伏已无出力 —— 最难工况）",
    "winter_we": "冬季周末",
    "summer_wd": "夏季工作日（午后高峰，与光伏出力重叠）",
    "summer_we": "夏季周末",
    "springfall_wd": "春秋工作日（早晚双峰）",
    "springfall_we": "春秋周末",
}


def normalize_shape(vals):
    """把 24 个百分数转成日均 = 1.0 的形状系数（模型会再归一化一次）"""
    f = [v / 100.0 for v in vals]          # 先转成"占当日峰值"的比例
    m = sum(f) / len(f)
    return [v / m for v in f]


def write_one(path, shape, variant):
    with open(path, "w", newline="", encoding="utf-8") as fp:
        w = csv.writer(fp)
        w.writerow(["hour", "load_pu"])
        for h, v in enumerate(shape, 1):
            w.writerow([h, "%.6f" % v])
    print("已写出 %s  (%s)" % (path, DESC[variant]))


def write_variants(path, shapes):
    with open(path, "w", newline="", encoding="utf-8") as fp:
        w = csv.writer(fp)
        w.writerow(["hour"] + list(shapes.keys()))
        for h in range(24):
            w.writerow([h + 1] + ["%.6f" % shapes[k][h] for k in shapes])
    print("已写出 %s  （6 种变体并排，供切换）" % path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", default="winter_wd", choices=list(RTS79_TABLE3.keys()))
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()

    if args.list:
        print("可选变体：")
        for k, d in DESC.items():
            print("  %-16s %s" % (k, d))
        return

    shapes = {k: normalize_shape(v) for k, v in RTS79_TABLE3.items()}
    os.makedirs(OUT_DIR, exist_ok=True)

    write_one(os.path.join(OUT_DIR, "load_profile.csv"), shapes[args.variant], args.variant)
    write_variants(os.path.join(OUT_DIR, "load_profile.variants.csv"), shapes)

    s = shapes[args.variant]
    print("\n本次生效曲线的特征：")
    print("  峰值出现在第 %d 小时，峰谷比 %.2f" % (s.index(max(s)) + 1, max(s) / min(s)))
    print("  日均 = %.4f（模型会再归一化，无需手工调整）" % (sum(s) / len(s)))
    print("\n注意：换曲线会改变承载力、网损最低点、鲁棒优化等结论，需重跑 P2~P5 并对比。")


if __name__ == "__main__":
    main()
