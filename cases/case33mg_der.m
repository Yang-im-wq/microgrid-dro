function [mpc, der] = case33mg_der(varargin)
%CASE33MG_DER  IEEE 33 节点微电网算例：光伏 + 储能 + 支路热限值
%
%   在 MATPOWER 自带的 case33mg（Kashem et al., 33 节点配电系统）基础上叠加：
%     1. 光伏机组  —— 挂在馈线末端节点，单位功率因数，可弃光（Pmin=0）
%     2. 电池储能  —— 以"负 PMIN 机组"形式接入（见下方说明），可充可放
%     3. 支路热限值 —— 原始 case33mg 的 RATE_A 全为 0（无容量约束）
%
%   mpc               = case33mg_der()
%   mpc               = case33mg_der('name', value, ...)
%   [mpc, der]        = case33mg_der(...)
%
%   选项（name-value）：
%     'pv_penetration'  光伏装机 / 总负荷            (默认 0.50)
%     'pv_buses'        光伏接入母线                 (默认 [18 22 25 33])
%     'pv_available'    光伏可用出力系数 0~1         (默认 1.0)
%                       —— 潮流用：Pg = 装机 × 系数
%                       —— OPF 用：Pmax = 装机（Pmin=0，由成本/约束决定弃光）
%     'storage'         是否加储能                   (默认 true)
%     'storage_bus'     储能接入母线                 (默认 18)
%     'storage_power'   储能功率容量 MW              (默认 1.0)
%     'storage_energy'  储能能量容量 MWh             (默认 3.0)
%     'thermal'         是否加支路热限值             (默认 true)
%     'i_rated'         导体载流量 A（用于算 RATE_A）(默认 300)
%     'vmax'            电压上限 p.u.                 (默认 [] = 保持 1.1)
%     'vmin'            电压下限 p.u.                 (默认 [] = 保持 0.9)
%     'ramp'            机组爬坡率                   (默认 [] = 保持原值)
%                       —— []     保持原值（case33mg 全为 0，见下方说明）
%                       —— 标量   所有机组 RAMP_AGC/10/30 设为该值
%                       —— 'pmax' 每台机组 RAMP_30 = 其 PMAX
%     'load_scale'      负荷整体缩放系数             (默认 1.0)
%
%   输出：
%     mpc  MATPOWER 算例结构体（可直接喂 runpf / runopf）
%     der  派生信息结构体（供 MOST / DRO 层使用）：
%            .pv_buses .pv_cap_each .pv_cap_total .pv_penetration
%            .storage (struct 或 []) .i_rated .rate_a_mva
%
%   ★ 回归闸门用法（用于验证没改坏标准算例）：
%       mpc = case33mg_der('pv_penetration', 0, 'storage', false, 'thermal', false);
%       此时 bus / branch / gen / gencost 应与 case33mg 逐字段一致。
%
%   ★ 关于热限值（重要说明）
%       原始 case33mg 的 RATE_A 全为 0，即**没有任何容量约束**。后果是 DC 潮流下
%       零成本光伏必然顶到 Pmax，永远看不到阻塞与弃光。为了让"阻塞/弃光"这个现象
%       能够出现，本算例按导体载流量折算出一个**工程化的热限值**：
%           RATE_A = sqrt(3) × 12.66 kV × I_rated        (baseMVA=1，故 MVA 即 p.u.)
%       默认 I_rated = 300 A → 6.578 MVA，此时基态最大支路（branch 1）负载率 84%。
%       **这不是 IEEE 33 节点的文献标准值**，是对标准算例的刻意扩展，引用时须注明。
%
%   ★ 关于电压限值（重要说明）
%       case33mg 的 Vmax = 1.1 p.u.，这是**输电系统的默认值**，对 12.66 kV 配电网
%       过于宽松 —— ANSI C84.1 Range A 对配电网只允许到 1.05。用 1.1 会让光伏反送
%       导致的过电压出现得明显偏晚（实测本馈线约 240% 渗透率才越 1.1，而越 1.05
%       只需约 170%）。**配电网研究建议显式传 'vmax', 1.05**。
%       默认仍保持 1.1 以便与标准算例逐字段对标（回归闸门依赖这一点）。
%
%   ★ 关于储能参数的选择
%       储能容量必须相对负荷重新设定：MOST 自带示例 ex_storage.m 用 200 MWh，
%       相对本算例 3.715 MW 的总负荷是 54 倍，照抄会严重失真。
%       本算例默认 3 MWh / 1 MW（3 小时时长），约为峰荷的 27%。
%       储能在 MATPOWER 里就是一行普通机组，只是 PMIN 为负：
%           PMAX = +Pcap  → 放电上限
%           PMIN = -Pcap  → 充电上限（负出力 = 吸收功率）
%       MOST 另需 xGenData / StorageData 两张表来描述容量、效率、SOC，
%       由 src/build_most_data.m 从本函数返回的 der 结构体生成（见 P4）。
%
%   ★ 关于爬坡率（用 MOST 做时序调度时必读）
%       case33mg 的 RAMP_AGC / RAMP_10 / RAMP_30 **全为 0**。单时段潮流/OPF 无所谓，
%      但多时段调度（MOST）里爬坡率 0 意味着机组在一个时段内**完全无法调整出力**，
%      而负荷在变 —— LP 直接不可行，MOST 只报一句 "Numerically Failed"，
%      不告诉你原因。这是用 case33mg 做时序调度的隐藏陷阱（实测踩过）：
%          RAMP_30 = 0  -> 失败
%          RAMP_30 = 1  -> 失败（爬不过来）
%          RAMP_30 = 3  -> 成功
%      配电网研究里变电站侧通常不受爬坡限制，所以 build_most_data 用 'pmax'。
%
%   See also case33mg, addgen2mpc, runpf, runopf.

%% ---------- 选项解析 ----------
opt = struct( ...
    'pv_penetration', 0.50, ...
    'pv_buses',       [18 22 25 33], ...
    'pv_available',   1.00, ...
    'storage',        true, ...
    'storage_bus',    18, ...
    'storage_power',  1.00, ...
    'storage_energy', 3.00, ...
    'thermal',        true, ...
    'i_rated',        300, ...
    'vmax',           [], ...
    'vmin',           [], ...
    'ramp',           [], ...
    'load_scale',     1.00);

if mod(numel(varargin), 2) ~= 0
    error('case33mg_der: 选项必须成对给出 (name, value)');
end
for k = 1:2:numel(varargin)
    name = varargin{k};
    if ~ischar(name) || ~isfield(opt, name)
        error('case33mg_der: 未知选项 ''%s''', num2str(name));
    end
    opt.(name) = varargin{k + 1};
end

define_constants;

%% ---------- 从标准算例出发 ----------
mpc = case33mg;

if opt.load_scale ~= 1
    mpc.bus(:, [PD QD]) = mpc.bus(:, [PD QD]) * opt.load_scale;
end
total_load_MW = sum(mpc.bus(:, PD));

%% ---------- 0. 电压限值（可选；默认保持 case33mg 原值 1.1 / 0.9） ----------
if ~isempty(opt.vmax)
    mpc.bus(:, VMAX) = opt.vmax;
end
if ~isempty(opt.vmin)
    mpc.bus(:, VMIN) = opt.vmin;
end

%% ---------- 1. 光伏 ----------
n_pv       = numel(opt.pv_buses);
pv_total   = opt.pv_penetration * total_load_MW;   % MW，装机总量
pv_each    = pv_total / n_pv;                      % MW，每节点装机

if pv_each > 0
    if exist('addgen2mpc', 'file') ~= 2
        error(['case33mg_der: 找不到 addgen2mpc（属于 MOST）。' ...
               '请先运行 install_matpower 把 most/lib 加入路径。']);
    end
    gen = zeros(n_pv, 21);
    gen(:, GEN_BUS)     = opt.pv_buses(:);
    gen(:, PG)          = pv_each * opt.pv_available;  % 潮流用的实际注入
    gen(:, QG)          = 0;
    gen(:, QMAX)        = 0;      % 单位功率因数：不参与调压
    gen(:, QMIN)        = 0;
    gen(:, VG)          = 1.0;
    gen(:, MBASE)       = 100;
    gen(:, GEN_STATUS)  = 1;
    gen(:, PMAX)        = pv_each;   % 装机容量
    gen(:, PMIN)        = 0;         % 可弃光
    gencost = repmat([POLYNOMIAL 0 0 2 0 0], n_pv, 1);   % 零边际成本

    mpc = addgen2mpc(mpc, gen, gencost, 'solar');   % 自动建 mpc.isolar
end

%% ---------- 2. 储能 ----------
stor = [];
if opt.storage
    if exist('addgen2mpc', 'file') ~= 2
        error('case33mg_der: 找不到 addgen2mpc（属于 MOST）。');
    end
    pcap = opt.storage_power;
    ecap = opt.storage_energy;

    gen = zeros(1, 21);
    gen(1, [GEN_BUS GEN_STATUS VG MBASE]) = [opt.storage_bus 1 1 100];
    gen(1, PG)              = 0;
    gen(1, [QMAX QMIN])     = [0 0];
    gen(1, [PMAX PMIN])     = [pcap -pcap];   % 放电正 / 充电负
    gencost = [POLYNOMIAL 0 0 2 0 0];         % 零成本（储能的套利行为由 SOC 约束产生）

    mpc = addgen2mpc(mpc, gen, gencost, 'ess');   % 自动建 mpc.iess

    stor = struct('bus', opt.storage_bus, 'pcap', pcap, 'ecap', ecap, ...
                  'out_eff', 1.0, 'in_eff', 1.0, 'loss_factor', 0.0);
end

%% ---------- 3. 支路热限值 ----------
rate_a = 0;
if opt.thermal
    Vbase_kV = mpc.bus(1, BASE_KV);                     % 12.66 kV
    rate_a   = sqrt(3) * Vbase_kV * opt.i_rated / 1000; % MVA（baseMVA=1 → 即 p.u.）
    mpc.branch(:, RATE_A) = rate_a;                     % 联络开关同样赋值，无副作用
end

%% ---------- 4. 爬坡率 ----------
% case33mg 的 RAMP_* 全为 0，多时段调度会不可行（见文件头说明）
ramp_setting = 'unchanged';
if ~isempty(opt.ramp)
    if ischar(opt.ramp) && strcmpi(opt.ramp, 'pmax')
        mpc.gen(:, [RAMP_AGC RAMP_10 RAMP_30]) = repmat(mpc.gen(:, PMAX), 1, 3);
        ramp_setting = 'pmax';
    elseif isnumeric(opt.ramp) && isscalar(opt.ramp)
        mpc.gen(:, [RAMP_AGC RAMP_10 RAMP_30]) = opt.ramp;
        ramp_setting = sprintf('scalar %g', opt.ramp);
    else
        error('case33mg_der: ''ramp'' 只能是 []、标量，或字符串 ''pmax''');
    end
end

%% ---------- 派生信息 ----------
der = struct();
der.pv_buses       = opt.pv_buses(:)';
der.pv_cap_each    = pv_each;
der.pv_cap_total   = pv_total;
der.pv_penetration = opt.pv_penetration;
der.pv_available   = opt.pv_available;
der.storage        = stor;
der.i_rated        = opt.i_rated;
der.rate_a_mva     = rate_a;
der.vmax           = mpc.bus(1, VMAX);
der.vmin           = mpc.bus(1, VMIN);
der.ramp_setting   = ramp_setting;
der.total_load_MW  = total_load_MW;
der.load_scale     = opt.load_scale;

end
