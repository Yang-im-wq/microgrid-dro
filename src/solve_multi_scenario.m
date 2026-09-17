function [x, info] = solve_multi_scenario(L, dat, Xi, mode, varargin)
%SOLVE_MULTI_SCENARIO  多场景两阶段 LP（确定性 / 场景平均 / 分布鲁棒 DRO）
%
%   [x, info] = solve_multi_scenario(L, dat, Xi, mode)
%   [x, info] = solve_multi_scenario(L, dat, Xi, 'dro', 'eps', 0.02)
%
%   输入：
%     L     distflow_lin 返回的模型
%     dat   make_day_data 生成
%     Xi    (nt x npv x N) 各场景的光伏可用出力系数
%     mode  'det' 只用第 1 个场景 | 'saa' 场景平均 | 'dro' 平均 + ε·(最坏−最好)
%
%   ★ DRO 的模糊集与闭式化简
%     在概率权重上取 L1 球（Wasserstein 风格的常用形式）：
%         Q_ε = { q >= 0, sum(q)=1, || q - 1/N ||_1 <= 2ε }
%     内层最坏期望有闭式解：
%         max_{q in Q_ε} Σ q_i Q_i = mean_i Q_i + ε·( max_i Q_i − min_i Q_i )
%     推导：令 q = 1/N + Δ，则 ΣΔ = 0、‖Δ‖₁ ≤ 2ε，
%           max_Δ Σ Δ_i Q_i = ε(max Q − min Q)（预算全压到最高/最低两档之间转移）。
%     于是 min-max 退化成**标准 LP**：两个 epigraph 变量
%         T ≥ Q_i ∀i （逼 max）、Tm ≤ Q_i ∀i （逼 min），
%     目标加 ε(T − Tm) 即可。
%
%   ★ 性能：约束矩阵用**三元组索引 + sparse()** 组装
%     最初写成 `Aeq = [Aeq; row]` 逐行拼接，那是 O(n²) —— 12 场景下矩阵约
%     2 万行 × 1.16 万列，直接卡死（实测跑 10 分钟没出来）。
%     三元组写法毫秒级完成。
%
%   See also dro_24h, opt_dispatch_lindistflow, make_day_data.

opt = struct('eps', 0.02, 'gamma', 2000, 'voll', 500, 'curt_pen', 0.01, ...
             'vmax', 1.05, 'vmin', 0.90);
for k = 1:2:numel(varargin)
    name = varargin{k};
    if ~ischar(name) || ~isfield(opt, name)
        error('solve_multi_scenario: 未知选项 %s', num2str(name));
    end
    opt.(name) = varargin{k + 1};
end

nt = dat.nt;  nb = dat.nb;  npv = numel(dat.pv_bus);
N  = size(Xi, 3);
if strcmp(mode, 'det'), N = 1; end
v2min = opt.vmin^2;  v2max = opt.vmax^2;

%% ---------- 线性映射 ----------
Mv    = 2 * L.A_p * diag(L.r) * L.A_p';
dvQ   = 2 * L.A_p * (L.x .* (-L.A_p' * (-dat.Qd)));

B_slack = zeros(nb,1); B_slack(1) = 1;
B_ess   = zeros(nb,1); B_ess(dat.ess_bus) = 1;
B_pv    = zeros(nb, npv);
for k = 1:npv, B_pv(dat.pv_bus(k), k) = 1; end
cg = Mv*B_slack;  ce = Mv*B_ess;  cpv = Mv*B_pv;

%% ---------- 变量布局 ----------
iD = 1:nt;  iC = nt+(1:nt);  iE = 2*nt+(1:nt);
nS = nt + npv*nt + nb*nt + nt + nt;         % g, u, ls, sm, sx
base = 3*nt;
offG = 0;  offU = nt;  offLS = nt + npv*nt;  offSM = offLS + nb*nt;  offSX = offSM + nt;
iSbase = @(i) base + (i-1)*nS;
idxG  = @(i,t)   iSbase(i) + offG  + t;
idxU  = @(i,k,t) iSbase(i) + offU  + (t-1)*npv + k;
idxLS = @(i,j,t) iSbase(i) + offLS + (t-1)*nb  + j;
idxSM = @(i,t)   iSbase(i) + offSM + t;
idxSX = @(i,t)   iSbase(i) + offSX + t;

iT  = base + N*nS + 1;
iTm = iT + 1;
nvar = iTm;

%% ---------- 目标函数 ----------
f = zeros(nvar, 1);
cs = 1/N;
for i = 1:N
    b = iSbase(i);
    f(b + offG  + (1:nt))            = cs * dat.price(:);
    f(b + offLS + (1:(nb*nt)))       = cs * opt.voll;
    f(b + offSM + (1:nt))            = cs * opt.gamma;
    f(b + offSX + (1:nt))            = cs * opt.gamma;
    f(b + offU  + (1:(npv*nt)))      = cs * opt.curt_pen;
end
if strcmp(mode, 'dro')
    f(iT) = opt.eps;  f(iTm) = -opt.eps;
end

%% ---------- 约束：三元组累积 ----------
nrow_eq  = (nt + 1) + N*nt;
nrow_in  = N*nt*nb*2 + 2*(strcmp(mode,'dro'))*N;
nz_eq    = (nt+1)*3 + N*nt*(3 + npv + nb);
nz_in    = N*nt*nb*2*(3 + npv + nb) + 2*N*(nt + npv*nt + nb*nt + 2*nt);
% 留 20% 余量
eqI = zeros(round(nz_eq*1.2),1);  eqJ = eqI;  eqV = eqI;
inI = zeros(round(nz_in*1.2),1);  inJ = inI;  inV = inI;
eqN = 0; eqZ = 0;  inN = 0; inZ = 0;
beq = zeros(nrow_eq,1);  bineq = zeros(nrow_in,1);

eta = dat.ess_eta;

% ---- (1) SOC 动态 + 日循环 ----
for t = 1:nt
    eqN = eqN + 1;
    [eqZ, eqI, eqJ, eqV] = push(eqZ, eqI, eqJ, eqV, eqN, ...
        [iE(t), iC(t), iD(t)], [1, -eta, 1/eta]);
    if t > 1
        [eqZ, eqI, eqJ, eqV] = push(eqZ, eqI, eqJ, eqV, eqN, iE(t-1), -1);
    end
    beq(eqN) = (t == 1) * dat.ess_e0;
end
eqN = eqN + 1;
[eqZ, eqI, eqJ, eqV] = push(eqZ, eqI, eqJ, eqV, eqN, iE(nt), 1);
beq(eqN) = dat.ess_e0;

% ---- (2) 各场景功率平衡： g - u + d - c + shed = load - pmax ----
for i = 1:N
    for t = 1:nt
        eqN = eqN + 1;
        cols = [idxG(i,t), iD(t), iC(t), ...
                idxU(i,1:npv,t), idxLS(i,1:nb,t)];
        vals = [1, 1, -1, -ones(1,npv), ones(1,nb)];
        [eqZ, eqI, eqJ, eqV] = push(eqZ, eqI, eqJ, eqV, eqN, cols, vals);
        beq(eqN) = sum(dat.Pd(:,t)) - sum(squeeze(Xi(t,:,i))' .* dat.pv_cap(:));
    end
end

% ---- (3) 电压软约束 ----
for i = 1:N
    for t = 1:nt
        p_fix = -dat.Pd(:,t);
        for k = 1:npv
            p_fix = p_fix + B_pv(:,k) * Xi(t,k,i) * dat.pv_cap(k);
        end
        vfix = 1.0 - dvQ + Mv * p_fix;
        for j = 1:nb
            % 下限: -(Mv p_dec)_j - sm_t <= vfix_j - v2min
            inN = inN + 1;
            cols = [idxG(i,t), iD(t), iC(t), idxU(i,1:npv,t), ...
                    idxLS(i,1:nb,t), idxSM(i,t)];
            vals = [-cg(j), -ce(j), ce(j), cpv(j,:), -Mv(j,:), -1];
            [inZ, inI, inJ, inV] = push(inZ, inI, inJ, inV, inN, cols, vals);
            bineq(inN) = vfix(j) - v2min;

            % 上限: (Mv p_dec)_j - sx_t <= v2max - vfix_j
            inN = inN + 1;
            vals = [ cg(j),  ce(j), -ce(j), -cpv(j,:),  Mv(j,:), -1];
            [inZ, inI, inJ, inV] = push(inZ, inI, inJ, inV, inN, cols, vals);
            bineq(inN) = v2max - vfix(j);
        end
    end
end

% ---- (4) DRO epigraph： T ≥ cost_i, Tm ≤ cost_i ----
% 写成 "<= 0" 的标准形式：
%   T  ≥ cost_i  <=>  cost_i - T ≤ 0    -> 系数取 +cpos，T 上取 -1
%   Tm ≤ cost_i  <=>  Tm - cost_i ≤ 0   -> 系数取 -cpos，Tm 上取 +1
% ★ 这两个符号极易搞反。搞反后 Tm 没有上界 -> 目标 eps*(T-Tm) 趋于 -inf
%   -> LP 无界 -> linprog 返回空解（表现为"解不存在"但其实是建模错）。实测踩过。
if strcmp(mode, 'dro')
    for i = 1:N
        b = iSbase(i);
        ccols = [b + offG + (1:nt), b + offLS + (1:(nb*nt)), ...
                 b + offSM + (1:nt), b + offSX + (1:nt), b + offU + (1:(npv*nt))];
        cpos = [dat.price(:)', opt.voll*ones(1,nb*nt), ...
                opt.gamma*ones(1,nt), opt.gamma*ones(1,nt), ...
                opt.curt_pen*ones(1,npv*nt)];
        inN = inN + 1;
        [inZ, inI, inJ, inV] = push(inZ, inI, inJ, inV, inN, [ccols, iT], [cpos, -1]);
        bineq(inN) = 0;
        inN = inN + 1;
        [inZ, inI, inJ, inV] = push(inZ, inI, inJ, inV, inN, [ccols, iTm], [-cpos, 1]);
        bineq(inN) = 0;
    end
end

Aeq   = sparse(eqI(1:eqZ), eqJ(1:eqZ), eqV(1:eqZ), nrow_eq, nvar);
Aineq = sparse(inI(1:inZ), inJ(1:inZ), inV(1:inZ), nrow_in, nvar);

%% ---------- 变量上下界 ----------
lb = zeros(nvar, 1);
ub =  inf(nvar, 1);
ub(iD) = dat.ess_pmax;   ub(iC) = dat.ess_pmax;
ub(iE) = dat.ess_emax;   lb(iE) = dat.ess_emin;
for i = 1:N
    for t = 1:nt
        ub(idxG(i,t)) = dat.gmax;
        for k = 1:npv
            ub(idxU(i,k,t)) = Xi(t,k,i) * dat.pv_cap(k);
        end
        for j = 1:nb
            ub(idxLS(i,j,t)) = dat.Pd(j,t);
        end
    end
end
lb(iTm) = -inf;

%% ---------- 求解 ----------
lpopts = optimoptions('linprog', 'Display', 'off', 'OptimalityTolerance', 1e-9, ...
                      'ConstraintTolerance', 1e-8);
[x, fval, exitflag, output] = linprog(f, Aineq, bineq, Aeq, beq, lb, ub, lpopts);

cost_scen = nan(N,1);
if exitflag == 1
    for i = 1:N
        b = iSbase(i);
        cost_scen(i) = dat.price(:)' * x(b+offG+(1:nt)) ...
            + opt.voll * sum(x(b+offLS+(1:(nb*nt)))) ...
            + opt.gamma * (sum(x(b+offSM+(1:nt))) + sum(x(b+offSX+(1:nt)))) ...
            + opt.curt_pen * sum(x(b+offU+(1:(npv*nt))));
    end
end

info = struct('status', exitflag, 'obj', fval, 'cost_scen', cost_scen, ...
              'mode', mode, 'nvar', nvar, 'N', N, ...
              'x_first', x(iD(1):iE(end)), ...
              'idx', struct('D',iD,'C',iC,'E',iE), ...
              'iSbase', iSbase, 'off', struct('G',offG,'U',offU,'LS',offLS,'SM',offSM,'SX',offSX), ...
              'nt', nt, 'nb', nb, 'npv', npv);
if exitflag ~= 1
    warning('solve_multi_scenario(%s): exitflag = %d (%s)', mode, exitflag, output.message);
else
    fprintf('    [%s] nvar=%d, eq=%d, ineq=%d, exitflag=1, obj=%.2f\n', ...
        mode, nvar, nrow_eq, nrow_in, fval);
end
end


function [nz, I, J, V] = push(nz, I, J, V, row, cols, vals)
% 往三元组里追加一行
k = numel(cols);
if nz + k > numel(I)
    grow = max(numel(I), 1e5);
    I = [I; zeros(grow,1)];  J = [J; zeros(grow,1)];  V = [V; zeros(grow,1)];
end
I(nz+1:nz+k) = row;
J(nz+1:nz+k) = cols(:);
V(nz+1:nz+k) = vals(:);
nz = nz + k;
end
