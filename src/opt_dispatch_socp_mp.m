function [x, o] = opt_dispatch_socp_mp(L, dat, Xi, varargin)
%OPT_DISPATCH_SOCP_MP  多时段二阶锥 (SOCP) 调度（含储能）
%
%   [x, o] = opt_dispatch_socp_mp(L, dat, Xi)
%
%   输入：
%     L     distflow_lin 返回的模型
%     dat   make_day_data 生成
%     Xi    (nt x npv) 各时段各光伏节点的可用出力系数
%
%   选项：'gamma' 电压越限惩罚 $/p.u.² (默认 2000)
%         'voll'  切负荷惩罚 $/MWh        (默认 500)
%         'vmax'/'vmin'                   (默认 1.05 / 0.90)
%         'vroot' 'fixed'(默认) | 'free'  —— 变电站电压是否作为决策变量
%         'curt_pen' 弃光惩罚 $/MWh       (默认 0.01)
%
%   ★ 与单时段版 (opt_dispatch_socp) 的关系
%     网络方程完全一致（含线损项 + 二阶锥），区别是：
%       1. 每个时段复制一套网络变量
%       2. 增加跨时段耦合的储能变量 d/c/e（单时段版没有储能，因为退化）
%     这样优化就能**看见线损**。LinDistFlow 看不见，它的调度会系统性偏乐观，
%     尤其在高渗透率反送、线损占比大的时候。
%
%   ★ 'vroot' = 'free' 的意义
%     单时段版把变电站电压固定为 1.0。但真实变电站（有载调压变压器）可以调压，
%     交流最优潮流会把它顶到上限去降损。放开它是一个自由度，代价是：
%     优化器可能用"抬高全网电压"来减少电流、降低线损 —— 这是真实存在的手段，
%     但也会压缩电压上限的裕度。
%
%   See also opt_dispatch_socp, optimize_storage_siting, make_day_data.

opt = struct('gamma', 2000, 'voll', 500, 'curt_pen', 0.01, ...
             'vmax', 1.05, 'vmin', 0.90, 'vroot', 'fixed', 'sm_max', 0.2);
for k = 1:2:numel(varargin)
    name = varargin{k};
    if ~ischar(name) || ~isfield(opt, name)
        error('opt_dispatch_socp_mp: 未知选项 %s', num2str(name));
    end
    opt.(name) = varargin{k + 1};
end

nt = dat.nt;  nb = L.nb;  nbr = L.nbr;  root = L.root;
npv = numel(dat.pv_bus);
v2min = opt.vmin^2;  v2max = opt.vmax^2;
free_vr = strcmpi(opt.vroot, 'free');

%% ---------- 拓扑定向（与单时段版同一套逻辑）----------
fdir = L.fb;  tdir = L.tb;
seen = false(nb,1);  seen(root) = true;
parent_of = zeros(nb,1);
queue = root;
while ~isempty(queue)
    v = queue(1); queue(1) = [];
    for k = 1:nbr
        if fdir(k) == v
            w = tdir(k);  rev = false;
        elseif tdir(k) == v
            w = fdir(k);  rev = true;
        else
            continue
        end
        if ~seen(w)
            if rev, fdir(k) = v; tdir(k) = w; end
            seen(w) = true;  parent_of(w) = k;
            queue(end+1) = w;                  %#ok<AGROW>
        end
    end
end
assert(all(seen), 'opt_dispatch_socp_mp: 网络不连通');
children_of = cell(nb,1);
for k = 1:nbr, children_of{fdir(k)}(end+1) = k; end
r = L.r;  x = L.x;

%% ---------- 变量布局 ----------
iD = 1:nt;  iC = nt+(1:nt);  iE = 2*nt+(1:nt);
nS = 3*nbr + nb + 1 + npv + nb + 2;        % 每时段的网络块
base = 3*nt;
offP = 0; offQ = nbr; offL = 2*nbr; offV = 3*nbr;
offG = 3*nbr + nb;  offU = offG + 1;  offLS = offU + npv;  offSM = offLS + nb;
offSX = offSM + 1;
iVR = base + nt*nS + 1;                     % 变电站电压平方（可选）
nvar = iVR;

blk  = @(t) base + (t-1)*nS + (1:nS);
iPf  = @(t,k) base + (t-1)*nS + offP + k;
iQf  = @(t,k) base + (t-1)*nS + offQ + k;
iLf  = @(t,k) base + (t-1)*nS + offL + k;
iVf  = @(t,j) base + (t-1)*nS + offV + j;
iGf  = @(t)   base + (t-1)*nS + offG + 1;
iUf  = @(t,k) base + (t-1)*nS + offU + k;
iLSf = @(t,j) base + (t-1)*nS + offLS + j;
iSMf = @(t)   base + (t-1)*nS + offSM + 1;
iSXf = @(t)   base + (t-1)*nS + offSX + 1;

%% ---------- 目标 ----------
% 注意：多时段版的变量都按时段索引（iGf(t) 而不是 iG），不能照抄单时段版
f = zeros(nvar,1);
for t = 1:nt
    f(iGf(t)) = dat.price(t);
    f(iLSf(t,1:nb)) = opt.voll;
    f(iUf(t,1:npv)) = opt.curt_pen;
    f([iSMf(t), iSXf(t)]) = opt.gamma;
end

%% ---------- 约束：三元组 ----------
% 等式：储能动态 nt 条 + 日循环 1 条 + 每时段 [P平衡(nb-1) + 根(1) + Q平衡(nb-1)
%       + 电压(nb-1) + 根电压(1)] = 3(nb-1)+2 条
% 预分配用估计值，循环结束后按**实际条数**截断，避免计数写错导致维度不匹配
%（最初写成 5(nb-1)+2，结果 Aeq 与 beq 行数对不上，coneprog 直接报错）。
n_eq = (nt + 1) + nt*(3*(nb-1) + 2);
nz_eq = nt*(nb-1)*2*(4+nb+npv) + nt*40 + 1;
eqI = zeros(nz_eq,1); eqJ = eqI; eqV = eqI;  neqN = 0; nzE = 0;
beq = zeros(n_eq,1);
eqrow = 0;

% 储能动态
for t = 1:nt
    eqrow = eqrow + 1;
    cols = [iE(t), iC(t), iD(t)];
    vals = [1, -dat.ess_eta, 1/dat.ess_eta];
    if t > 1, cols(end+1) = iE(t-1); vals(end+1) = -1; end
    [nzE, eqI, eqJ, eqV] = push(nzE, eqI, eqJ, eqV, eqrow, cols, vals);
    beq(eqrow) = (t == 1) * dat.ess_e0;
end
% 日循环
eqrow = eqrow + 1;
[nzE, eqI, eqJ, eqV] = push(nzE, eqI, eqJ, eqV, eqrow, iE(nt), 1);
beq(eqrow) = dat.ess_e0;

% 逐时段的网络方程
for t = 1:nt
    % 有功平衡
    for j = 1:nb
        if j == root, continue; end
        kp = parent_of(j);
        eqrow = eqrow + 1;
        cols = [iPf(t,kp), iLf(t,kp), iLSf(t,j), iD(t), iC(t)];
        vals = [1, -r(kp), 1, 1, -1];
        for k = children_of{j}, cols(end+1) = iPf(t,k); vals(end+1) = -1; end
        rhs = dat.Pd(j,t);
        for k = 1:npv
            if dat.pv_bus(k) == j
                cols(end+1) = iUf(t,k);  vals(end+1) = -1;
                rhs = rhs - Xi(t,k) * dat.pv_cap(k);
            end
        end
        [nzE, eqI, eqJ, eqV] = push(nzE, eqI, eqJ, eqV, eqrow, cols, vals);
        beq(eqrow) = rhs;
    end
    % 根母线
    eqrow = eqrow + 1;
    cols = []; vals = [];
    for k = children_of{root}, cols(end+1) = iPf(t,k); vals(end+1) = 1; end
    cols(end+1) = iGf(t);  vals(end+1) = -1;
    [nzE, eqI, eqJ, eqV] = push(nzE, eqI, eqJ, eqV, eqrow, cols, vals);
    beq(eqrow) = 0;
    % 无功平衡
    for j = 1:nb
        if j == root, continue; end
        kp = parent_of(j);
        eqrow = eqrow + 1;
        cols = [iQf(t,kp), iLf(t,kp)];
        vals = [1, -x(kp)];
        for k = children_of{j}, cols(end+1) = iQf(t,k); vals(end+1) = -1; end
        [nzE, eqI, eqJ, eqV] = push(nzE, eqI, eqJ, eqV, eqrow, cols, vals);
        beq(eqrow) = dat.Qd(j);
    end
    % 电压降
    for j = 1:nb
        if j == root, continue; end
        kp = parent_of(j);
        eqrow = eqrow + 1;
        cols = [iVf(t,j), iVf(t,fdir(kp)), iPf(t,kp), iQf(t,kp), iLf(t,kp)];
        vals = [1, -1, 2*r(kp), 2*x(kp), -(r(kp)^2 + x(kp)^2)];
        [nzE, eqI, eqJ, eqV] = push(nzE, eqI, eqJ, eqV, eqrow, cols, vals);
        beq(eqrow) = 0;
    end
    % 根电压：固定 1.0 或作为变量
    eqrow = eqrow + 1;
    if free_vr
        cols = [iVf(t,root), iVR];  vals = [1, -1];
        [nzE, eqI, eqJ, eqV] = push(nzE, eqI, eqJ, eqV, eqrow, cols, vals);
        beq(eqrow) = 0;
    else
        [nzE, eqI, eqJ, eqV] = push(nzE, eqI, eqJ, eqV, eqrow, iVf(t,root), 1);
        beq(eqrow) = 1.0;
    end
end

% 电压软约束
n_in = nt*nb*2;
inI = zeros(round(n_in*12*1.3),1); inJ = inI; inV = inI; ninN = 0; nzI = 0;
bin = zeros(n_in,1);
for t = 1:nt
    for j = 1:nb
        ninN = ninN + 1;
        [nzI, inI, inJ, inV] = push(nzI, inI, inJ, inV, ninN, ...
            [iVf(t,j), iSMf(t)], [-1, -1]);
        bin(ninN) = -v2min;
        ninN = ninN + 1;
        [nzI, inI, inJ, inV] = push(nzI, inI, inJ, inV, ninN, ...
            [iVf(t,j), iSXf(t)], [1, -1]);
        bin(ninN) = v2max;
    end
end

Aeq   = sparse(eqI(1:nzE), eqJ(1:nzE), eqV(1:nzE), eqrow, nvar);
beq   = beq(1:eqrow);                    % 以实际条数为准，防计数写错
Aineq = sparse(inI(1:nzI), inJ(1:nzI), inV(1:nzI), n_in, nvar);
assert(size(Aeq,1) == numel(beq) && size(Aineq,1) == numel(bin), ...
    'opt_dispatch_socp_mp: 约束矩阵与右端项行数不一致');

%% ---------- 二阶锥 ----------
soc = [];
for t = 1:nt
    for k = 1:nbr
        A = zeros(3, nvar);
        A(1, iPf(t,k)) = 2;
        A(2, iQf(t,k)) = 2;
        A(3, iLf(t,k)) = 1;  A(3, iVf(t,fdir(k))) = -1;
        d = zeros(1, nvar);  d(iLf(t,k)) = 1;  d(iVf(t,fdir(k))) = 1;
        sk = secondordercone(A, zeros(3,1), d, 0);
        if isempty(soc), soc = sk; else, soc(end+1,1) = sk; end   %#ok<AGROW>
    end
end

%% ---------- 上下界 ----------
lb = zeros(nvar,1);
ub = inf(nvar,1);
for t = 1:nt
    lb(iPf(t,1:nbr)) = -inf;      % 允许反向潮流（光伏反送时 P<0）
    lb(iQf(t,1:nbr)) = -inf;
    % ★ 电压越限量加物理上界。原本留成无穷大，导致可行域极大、目标极浅，
    %   内点法容易在里面游荡不收敛。0.2 p.u.² 已经对应约 0.1 p.u. 的电压越限，
    %   再大就完全不是"可接受的运行点"了，限制它没有副作用。
    ub(iSMf(t)) = opt.sm_max;
    ub(iSXf(t)) = opt.sm_max;
end
ub(iD) = dat.ess_pmax;   ub(iC) = dat.ess_pmax;
ub(iE) = dat.ess_emax;   lb(iE) = dat.ess_emin;
for t = 1:nt
    ub(iGf(t)) = dat.gmax;
    for k = 1:npv, ub(iUf(t,k)) = Xi(t,k) * dat.pv_cap(k); end
    for j = 1:nb,  ub(iLSf(t,j)) = dat.Pd(j,t); end
end
if free_vr, lb(iVR) = v2min;  ub(iVR) = v2max; end

%% ---------- 求解 ----------
opts = optimoptions('coneprog', 'Display', 'off', ...
    'OptimalityTolerance', 1e-8, 'ConstraintTolerance', 1e-7, ...
    'MaxIterations', 400);
t0 = tic;
[x, fval, exitflag, output] = coneprog(f, soc, Aineq, bin, Aeq, beq, lb, ub, opts);
tw = toc(t0);

o = struct('exitflag', exitflag, 'fval', fval, 'output', output, ...
           'nvar', nvar, 'ncone', numel(soc), 'solve_time', tw, ...
           'xraw', x, ...                       % 原始解向量（供外部做交流校验用）
           'model', 'SOCP-MP', 'vroot_free', free_vr, ...
           'idx', struct('D',iD,'C',iC,'E',iE,'base',base,'nS',nS, ...
                         'offP',offP,'offQ',offQ,'offL',offL,'offV',offV, ...
                         'offG',offG,'offU',offU,'offLS',offLS, ...
                         'offSM',offSM,'offSX',offSX,'iVR',iVR), ...
           'nb', nb, 'nbr', nbr, 'nt', nt, 'npv', npv);
%% ---------- 解的内置校验（重要）----------
% coneprog 在本题上会报 exitflag = -7（未达指定容差）。**这个标志不能忽略**：
% 实测它返回的解可能严重违反功率平衡（最大残差 4.16 MW），表现为
% "高峰时段自报进口量为 0，而负荷明明有 4.5 MW"。
% 所以这里直接算一遍全网平衡残差，把结果放进 o.balance_res_*，并在超标时大声警告 ——
% 免得下游脚本把一个物理上不成立的解当结果用。
if any(~isnan(x))
    res = zeros(nt,1);
    for t = 1:nt
        b = base + (t-1)*nS;
        g   = x(b + offG + 1);
        u   = sum(x(b + offU + (1:npv)));
        ls  = sum(x(b + offLS + (1:nb)));
        dpv = sum(Xi(t,:) .* dat.pv_cap) - u;
        gen = g + dpv + x(iD(t)) - x(iC(t)) + ls;
        loss = sum(r .* x(b + offL + (1:nbr)));      % Σ r·l = 总有功线损
        res(t) = gen - (sum(dat.Pd(:,t)) + loss);
    end
    o.balance_res_max  = max(abs(res));
    o.balance_res_mean = mean(abs(res));
else
    o.balance_res_max = NaN;  o.balance_res_mean = NaN;
end
o.converged = (exitflag == 1) && (o.balance_res_max < 1e-4);

fprintf('  [SOCP-MP] nvar=%d, cones=%d, exitflag=%d, %.1fs, obj=%.2f\n', ...
    nvar, numel(soc), exitflag, tw, fval);
fprintf('            平衡残差 max=%.4f MW mean=%.4f MW  -> %s\n', ...
    o.balance_res_max, o.balance_res_mean, ...
    tern(o.converged, 'CONVERGED (可用)', 'NOT CONVERGED (解不可用)'));
if ~o.converged
    fprintf(['    ⚠ coneprog 未收敛，返回的解违反功率平衡，**不能当作结果使用**。\n' ...
             '      这是求解器层面的问题（缩放/收敛判据），不是模型本身的问题 ——\n' ...
             '      同一套方程在单时段"潮流模式"下已被交流潮流验证为电压精确\n' ...
             '      （见 compare_distflow_models）。多时段优化需调缩放或换求解器。\n' ...
             '      msg: %s\n'], output.message);
end
end


function s = tern(c, a, b)
if c, s = a; else, s = b; end
end


function [nz, I, J, V] = push(nz, I, J, V, row, cols, vals)
k = numel(cols);
if nz + k > numel(I)
    g = max(numel(I), 2e5);
    I = [I; zeros(g,1)];  J = [J; zeros(g,1)];  V = [V; zeros(g,1)];
end
I(nz+1:nz+k) = row;  J(nz+1:nz+k) = cols(:);  V(nz+1:nz+k) = vals(:);
nz = nz + k;
end
