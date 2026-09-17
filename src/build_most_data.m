function [mpc, xgd, sd, profiles, nt, meta] = build_most_data(varargin)
%BUILD_MOST_DATA  从 case33mg_der 装配 MOST 所需的全部输入
%
%   [mpc, xgd, sd, profiles, nt] = build_most_data()
%   [...] = build_most_data('pv_penetration', 1.0)
%
%   选项：
%     'pv_penetration' 光伏装机渗透率      (默认 1.5)
%     'pv_quantile'    用哪个分位数做日前预测 (默认 0.5，即中位数)
%     'hourly'         是否用小时profile    (默认 true，nt=24)
%     'vmax'/'vmin'                         (默认 1.05 / 0.90)
%     'load_shape'     负荷日曲线（24 点）  (默认内置典型日，见下)
%
%   输出：
%     mpc       MATPOWER 算例（含 PV + 储能）
%     xgd       xGenData（此处返回 []，由 loadmd 建默认值即可 —— 我们不设
%               备用/爬坡报价，纯经济调度）
%     sd        StorageData（**必须提供**：mpc 含 iess 时 loadmd 会检查它非空，
%               且 sd.UnitIdx 必须与 mpc.iess 对应）
%     profiles  MOST profile 数组（光伏可用率 + 负荷曲线）
%     nt        时段数
%     meta      元数据（记录数据来源，便于溯源）
%
%   ★ 数据来源说明（务必注意）
%     光伏曲线：**真实数据**，来自配套的光伏概率预测项目（LSTM 分位数回归）的
%               输出，经 bridge/export_pv.py 聚合成 24 点标幺日内曲线。
%     负荷曲线：**合成数据**。上游数据集只有光伏/气象列，没有负荷。
%               本函数用一条内置的典型日双峰曲线（早高峰 + 晚高峰），
%               归一化到日均 1.0。`meta.load_source` 会标明这一点。
%               拿它做定量结论前必须先替换成真实负荷数据。
%
%   See also run_most_24h, case33mg_der, loadmd, loadstoragedata.

%% ---------- 选项 ----------
opt = struct('pv_penetration', 1.5, 'pv_quantile', 0.5, 'hourly', true, ...
             'vmax', 1.05, 'vmin', 0.90, 'load_shape', [], ...
             'tou', false, 'tou_shape', [], ...
             'quadratic_cost', true, 'c2', 0.5, 'c1', 20, 'c0', 0);
for k = 1:2:numel(varargin)
    name = varargin{k};
    if ~ischar(name) || ~isfield(opt, name)
        error('build_most_data: 未知选项 ''%s''', num2str(name));
    end
    opt.(name) = varargin{k + 1};
end

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'cases'));
define_constants;

%% ---------- 1. 算例 ----------
mpc = case33mg_der('pv_penetration', opt.pv_penetration, ...
                   'vmax', opt.vmax, 'vmin', opt.vmin, ...
                   'ramp', 'pmax');   % case33mg 爬坡率为 0，多时段会不可行
nt = 24;

%% ---------- 2. 光伏可用率曲线（真实数据）----------
% pv_daily_profile_hourly.csv：24 行，q05..q95 列为标幺值（1.0 = 装机容量）
csvp = fullfile(root, 'data', 'pv_daily_profile_hourly.csv');
T = readtable(csvp, 'Delimiter', ',');
qlab = sprintf('q%02.0f', opt.pv_quantile * 100);
if ~ismember(qlab, T.Properties.VariableNames)
    error('build_most_data: 曲线文件里没有分位数列 %s', qlab);
end
pv_pu = T.(qlab);                       % 24 x 1，标幺值
if numel(pv_pu) ~= nt
    error('build_most_data: 期望 %d 个小时点，实际 %d', nt, numel(pv_pu));
end

%% ---------- 3. 负荷曲线（合成，见函数头说明）----------
if isempty(opt.load_shape)
    load_pu = default_load_shape();
    load_src = 'synthetic typical-day two-peak shape (normalized to mean 1.0)';
else
    load_pu = opt.load_shape(:);
    load_src = 'user-provided';
end
load_pu = load_pu / mean(load_pu);      % 归一化到日均 1.0

%% ---------- 4. profiles ----------
% idx_ct 的输出位置必须对齐：CT_TGEN 在第 5 位、CT_ROW 第 10、CT_COL 第 11、
% CT_CHGTYPE 第 12、CT_REP 第 13、CT_REL 第 14、CT_TLOAD 第 17、
% CT_LOAD_ALL_PQ 第 19、CT_TGENCOST 第 25、CT_MODCOST_F 第 27
% （见 most/examples/ex_load_profile.m 的完整列表）。
% 用 ~ 跳过的位置也要占位，否则后面的变量会错位。
[~, ~, ~, ~, CT_TGEN, ~, ~, ~, ~, CT_ROW, CT_COL, CT_CHGTYPE, CT_REP, ...
 CT_REL, ~, ~, CT_TLOAD, ~, CT_LOAD_ALL_PQ, ~, ~, ~, ~, ~, CT_TGENCOST, ...
 ~, CT_MODCOST_F] = idx_ct;

% ★★ 关键写法（踩过坑，务必注意拼接方向）：
%   1) 每台光伏**单独一个 profile**，rows 用标量、values 用二维
%      （官方 most_ex5_mpopf / ex_wind_profile_d 就是这么写的）
%   2) **必须拼成列向量**（用 [a; b]，不是 [a, b]）
%      因为 loadmd.m:274 是用 `nprof = size(profiles,1)` 来数 profile 个数的
%      —— 看的是**行数**，不是 numel。拼成 1xN 行向量的话 nprof 恒等于 1，
%      于是**只有第一个 profile 会被处理，其余全部静默丢弃**（没有任何报错）。
%      实测表现：contab 里只剩第一个 profile 的行，于是"负荷曲线不生效"或
%      "光伏曲线不生效"，取决于哪个写在前面。这个坑极难查。
profiles = [];
for k = 1:numel(mpc.isolar)
    pk = struct('type', 'mpcData', 'table', CT_TGEN, ...
        'rows', mpc.isolar(k), 'col', PMAX, 'chgtype', CT_REL, 'values', pv_pu);
    profiles = [profiles; pk];                                          %#ok<AGROW>
end

% 负荷：整条馈线按日曲线缩放（rows = 0 表示全部母线）
ldprof = struct('type', 'mpcData', 'table', CT_TLOAD, ...
    'rows', 0, 'col', CT_LOAD_ALL_PQ, 'chgtype', CT_REL, 'values', load_pu);
profiles = [profiles; ldprof];

% 分时电价：对平衡机（gen 第 1 行）的发电成本做纵向缩放
% ⚠️ 实测**不生效**（2026-09-17）：即使把乘数全设成 3.0，MOST 的目标值也纹丝不动
%    （1254.14 -> 1254.14）。`apply_changes.m:182-202` 确实实现了 CT_TGENCOST +
%    CT_MODCOST_F 的处理，但这条链路在 MOST 里走不通。**所以默认关闭。**
%    让储能动起来的有效手段是上面的**二次成本**（边际成本递增 -> 削峰有价值）。
%    保留这段代码是为了记录这个坑，不是为了使用。
if opt.tou
    if isempty(opt.tou_shape)
        tou = default_tou_shape();
        tou_src = 'synthetic TOU: off-peak 0.6 / flat 1.0 / peak 1.8';
    else
        tou = opt.tou_shape(:);
        tou_src = 'user-provided';
    end
    tou = tou(:) / mean(tou);        % 归一化到日均 1.0，不改变日均电价水平
    touprof = struct('type', 'mpcData', 'table', CT_TGENCOST, ...
        'rows', 1, 'col', CT_MODCOST_F, 'chgtype', CT_REL, 'values', tou);
    profiles = [profiles; touprof];
else
    tou_src = 'none (flat $20/MWh)';
end

%% ---------- 5. 平衡机成本（决定储能有没有用武之地）----------
% case33mg 的 gencost 是 **纯线性** [2 0 0 3 0 20 0]（20 $/MW，边际成本恒定）。
% 线性成本下**削峰没有任何价值** —— 少发 1 MWh 就省 20 元，多发 1 MWh 就多花 20 元，
% 储能充放电的套利空间为 0，于是它完全不动（实测 throughput 仅 0.358 MWh）。
%
% 真实机组的边际成本随出力递增，用二次成本更贴切，也让储能真正有用武之地：
% 低谷时多充、高峰时放电替代高边际成本的平衡机出力。
% 默认 c2 = 0.5 $/MW^2·h（在 3.7 MW 出力下边际成本 ≈ 20 + 0.5*2*3.7 ≈ 23.7 $/MWh）。
if opt.quadratic_cost
    mpc.gencost(1, :) = [POLYNOMIAL 0 0 3, opt.c2, opt.c1, opt.c0];
end

%% ---------- 6. StorageData（必须提供）----------
pcap = mpc.gen(mpc.iess, PMAX);         % 放电上限 MW
ecap = 3.0;                             % 与 case33mg_der 默认一致：3 MWh
eff  = 0.95;                            % 单程效率

sdt = struct();
sdt.colnames = {'InitialStorage', 'InitialStorageLowerBound', ...
    'InitialStorageUpperBound', 'InitialStorageCost', 'TerminalStoragePrice', ...
    'MinStorageLevel', 'MaxStorageLevel', 'OutEff', 'InEff', 'LossFactor', 'rho'};
sdt.data = [ ...
    0.5*ecap, 0, ecap, 20, 20, ...      % 初始 SOC 50%，初/终值储能作价 20 $/MWh
    0.10*ecap, ecap, ...                % SOC 允许区间 [10%, 100%]
    eff, eff, 0.001, 0];                % 效率、自损耗、rho=0(按期望)
sd = loadstoragedata(sdt, mpc.gen(mpc.iess, :));
% 注意：第二个参数必须只传**储能那几行**（不是整个 mpc.gen）。
% loadstoragedata 要求 data 表行数 == 传入 GEN 的行数；传整个 mpc 会报
% "# of rows in 'data' table (1) do not equal rows in GEN (6)"。
% addstorage 内部也是这么传的：loadstoragedata(sd_table, storage.gen)。
sd.UnitIdx = mpc.iess;                  % loadmd 会校验它和 mpc.iess 长度一致

xgd = [];                               % 交给 loadmd 建默认值（纯 ED，无备用报价）

%% ---------- 6. 元数据 ----------
meta = struct();
meta.pv_source = 'REAL - LSTM quantile forecast via bridge/export_pv.py';
meta.pv_quantile = opt.pv_quantile;
meta.pv_csv = csvp;
meta.load_source = load_src;
meta.tou_source = tou_src;
meta.pv_penetration = opt.pv_penetration;
meta.pv_cap_total_MW = opt.pv_penetration * sum(case33mg().bus(:, PD));
meta.load_mean_MW = sum(case33mg().bus(:, PD));
meta.nt = nt;
% profiles 是列数组：[4 个光伏; 1 个负荷(: 可能还有 TOU)]。
% 用索引而不是硬编码 1/2 —— 改成列数组后 profiles(2) 是第二个光伏，不是负荷（踩过）。
meta.i_pv_profile = 1;
meta.i_load_profile = numel(mpc.isolar) + 1;
meta.storage_MWh = ecap;
meta.storage_MW = pcap;
meta.storage_eff = eff;

end


function s = default_tou_shape()
%DEFAULT_TOU_SHAPE  分时电价乘数曲线（合成）
%   谷段 0.6（0-5 时、23 时）、平段 1.0（7-16 时）、峰段 1.8（18-21 时）。
%   合成数据 —— 上游数据集没有电价列。形状参考常见工商业峰谷分时电价。
s = [0.6 0.6 0.6 0.6 0.6 0.7 0.9 1.0 1.0 1.0 1.0 1.0 ...
     1.0 1.0 1.0 1.0 1.2 1.8 1.8 1.8 1.8 1.2 0.8 0.6]';
end


function s = default_load_shape()
%DEFAULT_LOAD_SHAPE  典型日双峰负荷曲线（早高峰 + 晚高峰）
%   合成数据 —— 上游数据集没有负荷列。形状参考常见的工商业混合日负荷：
%   凌晨低谷、8-11 点早高峰、18-20 点晚高峰（最高）。
s = [0.62 0.58 0.55 0.54 0.55 0.62 0.72 0.82 0.89 0.92 0.92 0.90 ...
     0.87 0.86 0.88 0.92 1.00 1.12 1.20 1.22 1.13 1.00 0.84 0.71]';
end
