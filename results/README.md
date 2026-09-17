# 结果文件索引

所有 CSV 都由 `src/` 下的脚本生成，文件名的**后缀是参数标记** —— 这样同一组扫描
在不同参数下的结果可以并排对比，不会互相覆盖。

---

## 主线（P0–P5）

| 文件 | 生成脚本 | 内容 |
|---|---|---|
| `sweep_pv_vmax105.csv` | `sweep_pv_penetration.m` | 渗透率扫描（Vmax=1.05），51 行 × 16 列 |
| `sweep_pv_vmax110.csv` | 同上 | 同上，Vmax=1.10（用于对比限值敏感性） |
| `sweep_pv_vmax*.meta.json` | 同上 | 扫描参数元数据 |
| `reconfig_radial_pv000_vmax105.csv` | `reconfigure.m` | 径向重构的分支交换收敛轨迹 |
| `reconfig_meshed_pv000_vmax105.csv` | `reconfigure.m` | 32 组联络开关组合的环网运行结果 |
| `ed_sweep_pen200_vmax105.csv` | `run_ed.m` | 单时段 AC OPF 调度（PV 200%），弃光分析 |
| `most_24h_pen150_vmax105.csv` | `run_most_24h.m` | 24h MOST 调度 + 逐时段 AC 校验结果 |
| `most_24h_*.meta.json` | 同上 | 数据来源（光伏=真实 / 负荷=?) |

## DRO 三方对比（P5）

文件名格式：`dro_compare_pen{渗透率}_eps{鲁棒半径×1000}_g{电压惩罚}[_标记].csv`

**基础对比**（同一批 150 个样本外场景，不同设定）：

| 文件 | 设定 | 用途 |
|---|---|---|
| `..._g02000lo.csv` | 电压惩罚 γ=2000 | 低惩罚下的结果（越限严重） |
| `..._g50000hi.csv` | γ=50000 | 高惩罚（越限基本消除）—— **主结果** |
| `..._g50000rts79.csv` | 负荷曲线换成 IEEE RTS-79 后重跑 | 验证结论对负荷曲线的稳健性 |
| `..._g50000cloudtime.csv` | 场景加上云图实测的时间相关性 | 验证结论对时间相关的稳健性 |

**空间相关敏感性**（ρ = 光伏节点间相关系数，其余设定相同）：

| 文件 | ρ | 说明 |
|---|---|---|
| `..._g50000rho100.csv` | 1.0 | 完全相关（一片云盖住全部节点） |
| `..._g50000rho080.csv` | 0.8 | |
| `..._g50000rho050.csv` | 0.5 | |
| `..._g50000rho000.csv` | 0.0 | 完全独立 |

> 这四份一起读，得到"空间解耦降低尾部风险"的结论（CVaR95 降 0.9%）。
> 而横向比较 `..._g50000rho100` 与 `..._g50000cloudtime`，
> 得到"忽略时间相关性会低估尾部风险 4%"的结论。

## 扩展能力

| 文件 | 生成脚本 | 内容 |
|---|---|---|
| `storage_siting_pen150_g042.csv` | `optimize_storage_siting.m` | **选址 + 定容二维扫描**：32 母线 × 3 档容量 = 96 组合 |
| `distflow_model_compare.csv` | `compare_distflow_models.m` | LinDistFlow / SOCP / 全交流 三模型精度对比（同注入） |
| `mp_model_compare.csv` | `compare_mp_models.m` | 多时段模型对比（⚠️ SOCP-MP 未收敛，脚本会中止） |
| `vroot_compare.csv` | `test_vroot.m` | 变电站电压固定 vs 放开 |
| `dro_benefit_scan.csv` | `scan_dro_benefit.m` | **DRO 收益的二维扫描**（4 渗透率 × 4 储能容量 = 16 组，全部无收益） |

---

## 图件

见 `../figures/`。两张总览图：

- **`OVERVIEW.png`** —— 主线 P0→P5 八面板
- **`OVERVIEW_EXT.png`** —— 三项扩展九面板（选址 / SOCP / 云图）
