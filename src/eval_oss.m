function R = eval_oss(L, dat, x_first, Xi, varargin)
%EVAL_OSS  样本外评估：给定日前调度，在一批**没参与优化**的场景上评估表现
%
%   R = eval_oss(L, dat, x_first, Xi)
%
%   输入：
%     x_first  第一阶段决策（储能的 d / c / e，来自 solve_multi_scenario）
%     Xi       (nt x npv x N) 评估用场景（应与优化用场景不同）
%
%   输出 R 结构体：
%     .cost     (N x 1) 每场景总成本
%     .shed     (N x 1) 每场景切负荷量 MWh
%     .vmin     (N x 1) 每场景最低电压（LinDistFlow 口径，软约束）
%     .shed_exp 期望切负荷 MWh
%     .cost_exp 期望成本 $/day
%     .cost_cvar95  成本 95% CVaR（尾部风险的常用度量）
%     .viol_prob  电压越限场景占比
%     .worst_cost   最坏场景成本
%
%   ★ 为什么用"物理响应"而不是再解一次优化
%     日前调度一旦定下（储能怎么充放），实时就只能按实际光伏出力被动平衡：
%       需要进口 G_t = 负荷 - 光伏实际 + 充电 - 放电
%       进口不足以满足时 -> 切负荷
%     所以这里按**确定的物理规则**算，不需要再解 LP。这样评估是透明的，
%     而且三种策略（DET/SAA/DRO）用的是同一把尺子。
%
%   ★ 切负荷的分配规则
%     按各节点负荷比例分摊。这是一条规则，不是最优分配 —— 对三方对比是公平的，
%     但它不代表运行者一定会这么做。
%
%   See also dro_24h, solve_multi_scenario.

opt = struct('gamma', 2000, 'voll', 500, 'vmax', 1.05, 'vmin', 0.90);
for k = 1:2:numel(varargin)
    name = varargin{k};
    if ~isfield(opt, name), error('eval_oss: 未知选项 %s', num2str(name)); end
    opt.(name) = varargin{k + 1};
end

nt = dat.nt; nb = dat.nb; npv = numel(dat.pv_bus);
N  = size(Xi, 3);
d  = x_first(1:nt);
c  = x_first(nt + (1:nt));

R.cost = nan(N,1);  R.shed = nan(N,1);  R.vmin = nan(N,1);
R.vmax = nan(N,1);  R.curt = nan(N,1);

for s = 1:N
    tot = 0; shd = 0; curt = 0;
    vmin_s = inf; vmax_s = -inf;
    for t = 1:nt
        % 光伏可用出力（本场景）
        pv_avail_k = Xi(t,:,s)' .* dat.pv_cap(:);        % npv x 1
        pv_avail   = sum(pv_avail_k);

        % 需要的进口量（不含弃光；弃光在本模型里没有收益，故不主动弃）
        G = sum(dat.Pd(:,t)) - pv_avail + c(t) - d(t);

        g_t   = min(G, dat.gmax);                        % 变电站最多送这么多
        shed_t = max(G - g_t, 0);                        % 不够就切负荷
        if shed_t > sum(dat.Pd(:,t))
            shed_t = sum(dat.Pd(:,t));                   % 最多切光
        end

        % 按负荷比例分摊切负荷
        if shed_t > 0 && sum(dat.Pd(:,t)) > 0
            shed_vec = dat.Pd(:,t) / sum(dat.Pd(:,t)) * shed_t;
        else
            shed_vec = zeros(nb,1);
        end

        % 组装注入，跑 LinDistFlow 拿电压
        p_inj = zeros(nb,1);
        p_inj(1) = g_t;
        p_inj(dat.ess_bus) = p_inj(dat.ess_bus) + d(t) - c(t);
        for k = 1:npv
            p_inj(dat.pv_bus(k)) = p_inj(dat.pv_bus(k)) + pv_avail_k(k);
        end
        p_inj = p_inj - (dat.Pd(:,t) - shed_vec);
        [~, V] = distflow_eval(L, p_inj, -dat.Qd, 1.0);

        vmin_s = min(vmin_s, min(V));
        vmax_s = max(vmax_s, max(V));

        % 成本：购电 + 切负荷（VOLL） + 电压越限惩罚
        v_over  = max(0, max(V) - opt.vmax);
        v_under = max(0, opt.vmin - min(V));
        tot = tot + dat.price(t)*g_t + opt.voll*shed_t ...
                  + opt.gamma*(v_over^2 + v_under^2);
        shd = shd + shed_t;
    end
    R.cost(s) = tot;  R.shed(s) = shd;  R.vmin(s) = vmin_s;
    R.vmax(s) = vmax_s;  R.curt(s) = curt;
end

% ---- 汇总指标 ----
R.cost_exp   = mean(R.cost);
R.shed_exp   = mean(R.shed);
R.worst_cost = max(R.cost);
sc = sort(R.cost, 'descend');
k95 = max(1, ceil(0.05 * N));
R.cost_cvar95 = mean(sc(1:k95));            % 最差 5% 的平均
R.viol_prob = mean(R.vmin < opt.vmin - 1e-6 | R.vmax > opt.vmax + 1e-6);
R.vmin_exp  = mean(R.vmin);
end
