function [Xi, qlev, Q] = gen_pv_scenarios(root, nt, npv, N, seed, rho)
%GEN_PV_SCENARIOS  从 11 分位数的概率预测生成光伏场景（逆变换采样 + 空间相关）
%
%   Xi = gen_pv_scenarios(root, nt, npv, N, seed)
%   Xi = gen_pv_scenarios(root, nt, npv, N, seed, rho)
%
%   输入：
%     root  项目根目录
%     nt    时段数（24）
%     npv   光伏节点数
%     N     场景数
%     seed  随机种子（可复现）
%     rho   节点间相关系数，标量 0~1 或 npv×npv 相关矩阵（默认 1，即完全相关）
%
%   输出：
%     Xi    (nt x npv x N) 各场景各节点各时段的光伏可用出力系数 0~1
%     qlev  11 个分位数水平
%     Q     (nt x 11) 各时段的预测分位函数
%
%   ★ 空间相关性：为什么重要
%     光伏节点之间的相关性直接决定**总出力**的波动有多大：
%       完全相关 (rho = 1)：总波动 = 单个节点波动的 npv 倍（一条云盖住全部）
%       完全独立 (rho = 0)：总波动 = 单个节点波动的 sqrt(npv) 倍（分散化效应）
%     对 4 个节点，独立时总波动只有完全相关时的一半 ——
%     **波动减半意味着储能没那么必要、鲁棒优化的价值更小。**
%     所以 rho 不是个技术细节，它会直接改变"要不要做鲁棒优化"这个结论。
%
%   ★ 构造方法：高斯 Copula（保证边缘分布仍是那 11 档分位数）
%       1. 抽公共因子 Zc ~ N(0,1) 与独立因子 Zk ~ N(0,1)
%       2. 混合： Zk' = sqrt(rho)*Zc + sqrt(1-rho)*Zk      -> 节点间相关系数 = rho
%       3. 转均匀： Uk = Phi(Zk')   （Phi 是标准正态 CDF）
%       4. 逆变换： xi_k = Q_t(Uk)  （在预测的分位函数上插值）
%     这样各节点的**边缘分布**与预测完全一致，只是**联合分布**受 rho 控制。
%
%   See also dro_24h, build_most_data.

if nargin < 6 || isempty(rho), rho = 1.0; end

qlev = [0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95]';
qlab = arrayfun(@(q) sprintf('q%02.0f', q*100), qlev, 'UniformOutput', false);

csvp = fullfile(root, 'data', 'pv_daily_profile_hourly.csv');
T = readtable(csvp, 'Delimiter', ',');
Q = T{:, qlab};                      % nt x 11
assert(size(Q,1) == nt, 'gen_pv_scenarios: 曲线文件时段数与 nt 不一致');

% ---- 相关矩阵 ----
if isscalar(rho)
    if rho < -1e-9 || rho > 1+1e-9
        error('gen_pv_scenarios: 标量 rho 必须落在 [0,1]');
    end
    R = (1-rho) * eye(npv) + rho * ones(npv);
else
    R = rho;
    if size(R,1) ~= npv || size(R,2) ~= npv
        error('gen_pv_scenarios: 相关矩阵维度应为 %d x %d', npv, npv);
    end
end
% Cholesky 分解（相关矩阵需半正定；数值上给对角线一点保护）
[L, p] = chol(R, 'lower');
if p ~= 0
    R = R + 1e-10 * eye(npv);
    L = chol(R, 'lower');
end

rng(seed);                           % 固定种子，保证结果可复现

if abs(rho - 1) < 1e-12
    % ---- 完全相关：与旧行为一致（单条公共因子驱动全部节点）----
    U = rand(nt, N);
    Xi1 = zeros(nt, npv, N);
    for t = 1:nt
        Xi1(t, :, :) = repmat(interp1(qlev, Q(t,:), U(t,:), 'linear', 'extrap'), [npv 1]);
    end
else
    % ---- 一般情形：高斯 Copula ----
    % Z = Zind * L'  使得 Cov(Z) = L*L' = R（L 是相关矩阵的 Cholesky 下三角）
    Zi = randn(nt * N, npv);                     % 独立标准正态
    Z2 = Zi * L';                                % (nt*N) x npv，列间相关系数 = R
    Z  = permute(reshape(Z2, nt, N, npv), [1 3 2]);
    U  = normcdf(Z);                             % (nt x npv x N) 均匀边缘

    Xi1 = zeros(nt, npv, N);
    for t = 1:nt
        for j = 1:npv
            Xi1(t, j, :) = interp1(qlev, Q(t,:), squeeze(U(t,j,:)), 'linear', 'extrap');
        end
    end
end

% 越界截断：可用率物理上不可能为负，也不可能超过装机
Xi = min(max(Xi1, 0), 1);
end
