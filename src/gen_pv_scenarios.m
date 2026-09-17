function [Xi, qlev, Q, meta] = gen_pv_scenarios(root, nt, npv, N, seed, rho_space, rho_time)
%GEN_PV_SCENARIOS  从 11 分位数的概率预测生成光伏场景（逆变换采样 + 时空相关）
%
%   Xi = gen_pv_scenarios(root, nt, npv, N, seed)
%   Xi = gen_pv_scenarios(root, nt, npv, N, seed, rho_space)
%   Xi = gen_pv_scenarios(root, nt, npv, N, seed, rho_space, rho_time)
%
%   输出：
%     Xi    (nt x npv x N) 各场景各节点各时段的光伏可用出力系数 0~1
%     qlev  11 个分位数水平
%     Q     (nt x 11) 各时段的预测分位函数
%     meta  相关性设置与来源说明
%
%   ★ 空间相关 rho_space
%     决定**总出力**的波动有多大：
%       完全相关 (1.0)：总波动 = 单节点波动的 npv 倍（一片云盖住全部）
%       完全独立 (0.0)：总波动 = 单节点波动的 sqrt(npv) 倍（分散化效应）
%     对 4 个节点，独立时总波动只有完全相关时的一半。
%
%   ★ 时间相关 rho_time（本版新增）
%     默认值从**卫星云图实测**得到（bridge/analyze_cloud_correlation.py）：
%     上游视觉模块输出的云量序列显示相邻小时云量只变化 0.037（中位数 0.011），
%     拟合 AR(1) 得 rho = 0.956。光伏出力直接受云影响，所以也应当强时间相关。
%
%     之前这里是**逐小时独立抽样**，等于假设"这一小时和下一小时的光伏无关"，
%     与实测严重不符 —— 会让场景在小时之间乱跳，高估储能的调节压力。
%     现在默认用实测值；若 data/cloud_correlation.json 不存在则退回独立（rho_time=0）。
%
%   ★ 构造方法：高斯 Copula + 时空可分离相关
%     潜变量 Z 的相关阵取可分离形式
%         R = R_time (x) R_space        （Kronecker 积）
%     其中 R_time(t,t') = rho_time^|t-t'|，R_space(k,k') = rho_space (k≠k')。
%     Z = L * eps（L 为 R 的 Cholesky 因子），U = Phi(Z)，再在分位函数上逆变换。
%     **边缘分布严格等于预测的 11 档分位数**，相关性只改联合结构。
%
%   See also dro_24h, build_most_data, analyze_cloud_correlation.

root_data = fullfile(root, 'data');
if nargin < 6 || isempty(rho_space), rho_space = 1.0; end

% ---- 时间相关：优先读云图分析结果 ----
rho_src = 'independent (rho_time = 0, 逐小时独立)';
if nargin < 7 || isempty(rho_time)
    jf = fullfile(root_data, 'cloud_correlation.json');
    if exist(jf, 'file')
        j = jsondecode(fileread(jf));
        rho_time = j.recommended_rho_time;
        rho_src = sprintf('从卫星云图实测 (AR(1) 拟合, %s)', j.recommended_basis);
    else
        rho_time = 0;
    end
else
    rho_src = 'user-provided';
end

qlev = [0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95]';
qlab = arrayfun(@(q) sprintf('q%02.0f', q*100), qlev, 'UniformOutput', false);

csvp = fullfile(root_data, 'pv_daily_profile_hourly.csv');
T = readtable(csvp, 'Delimiter', ',');
Q = T{:, qlab};                      % nt x 11
assert(size(Q,1) == nt, 'gen_pv_scenarios: 曲线文件时段数与 nt 不一致');

meta = struct('rho_space', rho_space, 'rho_time', rho_time, ...
              'rho_time_source', rho_src, 'nt', nt, 'npv', npv, 'N', N);

% ---- 时空相关阵 R = R_time (x) R_space ----
if rho_time == 0 && rho_space == 1
    % 退化情形：空间完全相关、时间独立 —— 与旧版行为一致，省掉大矩阵
    rng(seed);
    U1 = rand(nt, N);
    Xi1 = zeros(nt, npv, N);
    for t = 1:nt
        Xi1(t,:,:) = repmat(interp1(qlev, Q(t,:), U1(t,:), 'linear', 'extrap'), [npv 1]);
    end
else
    Rsp = (1-rho_space) * eye(npv) + rho_space * ones(npv);
    Rtm = rho_time .^ abs((1:nt)' - (1:nt));         % Toeplitz
    R   = kron(Rtm, Rsp);                            % (nt*npv) x (nt*npv)
    R   = (R + R') / 2;                              % 数值对称化
    [L, p] = chol(R, 'lower');
    if p ~= 0                                        % 半正定保护
        R = R + 1e-10 * eye(size(R));
        L = chol(R, 'lower');
    end
    rng(seed);
    Zv = L * randn(nt*npv, N);                       % 时空相关的高斯潜变量
    Uv = normcdf(Zv);
    % ★ 复原顺序必须与 kron 的排列一致。
    %   MATLAB 的 kron(Rtm, Rsp) 里**第一个因子慢变、第二个快变**，
    %   所以线性索引是 (t-1)*npv + k —— 即"空间最快、时间最慢"。
    %   因此要先 reshape 成 (npv, nt, N) 再 permute 成 (nt, npv, N)。
    %   最初写成 reshape(Uv, nt, npv, N)（时间最快），顺序错位，
    %   表现为"设了空间相关 1.0 却只测出 0.63"。
    U  = permute(reshape(Uv, npv, nt, N), [2 1 3]);

    Xi1 = zeros(nt, npv, N);
    for t = 1:nt
        for j = 1:npv
            Xi1(t,j,:) = interp1(qlev, Q(t,:), squeeze(U(t,j,:)), 'linear', 'extrap');
        end
    end
end

Xi = min(max(Xi1, 0), 1);        % 可用率物理上落在 [0,1]
end
