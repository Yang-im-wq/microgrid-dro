function [Xi, qlev, Q] = gen_pv_scenarios(root, nt, npv, N, seed)
%GEN_PV_SCENARIOS  从 11 分位数的概率预测生成光伏场景（逆变换采样）
%
%   Xi = gen_pv_scenarios(root, nt, npv, N, seed)
%
%   输入：
%     root  项目根目录
%     nt    时段数（24）
%     npv   光伏节点数
%     N     场景数
%     seed  随机种子（可复现）
%
%   输出：
%     Xi    (nt x npv x N) 各场景各节点各时段的光伏可用出力系数 0~1
%     qlev  11 个分位数水平
%     Q     (nt x 11) 各时段的预测分位函数
%
%   方法：逆变换采样。对每个时段 t 抽 U~Uniform(0,1)，再在分位函数 Q_t(·) 上
%   线性插值得到该场景该时段的可用率。
%
%   ★ 空间相关性说明
%     本函数让**同一场景内所有光伏节点用同一个 U**，即假设空间完全相关
%     （一片云同时盖住四个节点）。这是保守且常用的简化：系统级功率平衡
%     面对的是总出力的波动，而空间完全相关会让总波动最大。
%     若要做空间解耦，给每个节点独立抽 U 即可。
%
%   See also dro_24h, build_most_data.

qlev = [0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95]';
qlab = arrayfun(@(q) sprintf('q%02.0f', q*100), qlev, 'UniformOutput', false);

csvp = fullfile(root, 'data', 'pv_daily_profile_hourly.csv');
T = readtable(csvp, 'Delimiter', ',');
Q = T{:, qlab};                      % nt x 11
assert(size(Q,1) == nt, 'gen_pv_scenarios: 曲线文件时段数与 nt 不一致');

rng(seed);                           % 固定种子，保证结果可复现
U = rand(nt, N);

Xi1 = zeros(nt, N);
for t = 1:nt
    Xi1(t,:) = interp1(qlev, Q(t,:), U(t,:), 'linear', 'extrap');
end
% 越界截断：可用率物理上不可能为负，也不可能超过装机
Xi1 = min(max(Xi1, 0), 1);

% 同一场景内所有 PV 节点共用同一条曲线（空间完全相关）
Xi = repmat(Xi1, [1 1 npv]);
Xi = permute(Xi, [1 3 2]);           % -> (nt x npv x N)
end
