function dat = make_day_data(varargin)
%MAKE_DAY_DATA  构造 24h 调度所需的数据结构（供 LinDistFlow 优化层使用）
%
%   dat = make_day_data()
%   dat = make_day_data('pv_penetration', 1.5, 'gmax', 4.2)
%
%   ★ 数据一致性
%     负荷曲线与光伏曲线**直接复用 build_most_data 生成的 profile**，
%     保证 P4（MOST）与 P5（DRO）用的是同一套输入，不会各写一套然后悄悄分叉。
%
%   字段：
%     nt                  时段数
%     Pd  (nb x nt)       有功负荷 MW（随时间变）
%     Qd  (nb x 1)        无功负荷 MVAr（不随时间变）
%     price (nt x 1)      购电电价 $/MWh（分时电价）
%     pv_bus / pv_cap     光伏接入点 / 各点装机 MW
%     ess_*               储能参数（接入点/功率/能量/初值/效率）
%     gmax                变电站进口上限 MW
%     pv_profile (nt x 1) 光伏可用率标幺曲线（=1.0 时满发）
%
%   ★ 关于电价
%     P4 里想用 MOST 的 CT_TGENCOST profile 做分时电价，实测**不生效**。
%     这里在 LP 目标函数里**直接**写分时电价 —— 自己的模型自己建，没有任何
%     中间层会把它吃掉。所以 P5 里分时电价是真的在工作。
%
%   See also opt_dispatch_lindistflow, dro_24h, build_most_data.

opt = struct('pv_penetration', 1.5, 'gmax', 4.2, 'vmax', 1.05, 'vmin', 0.90, ...
             'ess_emax', 3.0, 'ess_pmax', 1.0, 'use_tou', true);
for k = 1:2:numel(varargin)
    name = varargin{k};
    if ~ischar(name) || ~isfield(opt, name)
        error('make_day_data: 未知选项 ''%s''', num2str(name));
    end
    opt.(name) = varargin{k + 1};
end

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'src'), fullfile(root, 'cases'));
define_constants;

[mpc, ~, ~, profiles, nt, meta] = build_most_data( ...
    'pv_penetration', opt.pv_penetration, 'vmax', opt.vmax, 'vmin', opt.vmin);

load_pu  = profiles(meta.i_load_profile).values(:,1,1);    % nt x 1
pv_pu    = profiles(meta.i_pv_profile).values(:,1,1);      % nt x 1

dat = struct();
dat.nt      = nt;
dat.nb      = size(mpc.bus, 1);
dat.Pd      = mpc.bus(:, PD) * load_pu';        % nb x nt
dat.Qd      = mpc.bus(:, QD);                   % nb x 1
dat.pv_bus  = mpc.gen(mpc.isolar, GEN_BUS)';    % 1 x npv
dat.pv_cap  = mpc.gen(mpc.isolar, PMAX)';       % 1 x npv（各点装机 MW）
dat.pv_profile = pv_pu;
dat.ess_bus = mpc.gen(mpc.iess, GEN_BUS);
dat.ess_pmax = min(opt.ess_pmax, mpc.gen(mpc.iess, PMAX));
dat.ess_emax = opt.ess_emax;
dat.ess_emin = 0.10 * opt.ess_emax;
dat.ess_e0   = 0.50 * opt.ess_emax;             % 初值 50% SOC
dat.ess_eta  = 0.95;
dat.gmax     = opt.gmax;

if opt.use_tou
    tou = tou_price_shape();
    dat.price = 20 * tou(:) / mean(tou);        % 均值仍为 20 $/MWh
    dat.price_source = 'TOU (12/20/36 $/MWh), implemented directly in the LP';
else
    dat.price = 20 * ones(nt, 1);
    dat.price_source = 'flat 20 $/MWh';
end

% 光伏可用出力上限（MW），(nt x npv)
dat.pv_avail_max = pv_pu * dat.pv_cap;
end


function s = tou_price_shape()
%TOU_PRICE_SHAPE  分时电价乘数（合成）：谷 0.6 / 平 1.0 / 峰 1.8
%   形状与 build_most_data.default_tou_shape 一致（那里因为 MOST 吃不下而被弃用）。
s = [0.6 0.6 0.6 0.6 0.6 0.7 0.9 1.0 1.0 1.0 1.0 1.0 ...
     1.0 1.0 1.0 1.0 1.2 1.8 1.8 1.8 1.8 1.2 0.8 0.6]';
end
