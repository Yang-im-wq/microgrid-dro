function [x, o, DBG] = opt_dispatch_socp(L, dat, pv_avail, varargin)
%OPT_DISPATCH_SOCP  单时段配电网调度：二阶锥 (SOCP) DistFlow 模型
%
%   [x, o] = opt_dispatch_socp(L, dat, pv_avail)
%   [x, o] = opt_dispatch_socp(L, dat, pv_avail, 'hour', 19)
%
%   输入：
%     L         distflow_lin 返回的模型（借用其支路数据与拓扑）
%     dat       make_day_data 生成
%     pv_avail  (1 x npv) 各光伏节点的可用出力系数 0~1（单时段）
%
%   选项：'hour' 取哪一时段的负荷与电价（默认 19，即晚高峰）
%         'gamma' 电压越限惩罚 $/p.u.² (默认 2000)
%         'voll'  切负荷惩罚 $/MWh        (默认 500)
%         'vmax'/'vmin'                   (默认 1.05 / 0.90)
%         'curt_pen' 弃光惩罚 $/MWh       (默认 0.01)
%
%   ★ 与 LinDistFlow 的关键差别
%     LinDistFlow 把电压对注入线性化，代价是**丢掉线损项**：
%         v_j = v_i - 2(rP + xQ)                      <- LinDistFlow
%     完整的 DistFlow 多一项：
%         v_j = v_i - 2(rP + xQ) + (r²+x²)·l          <- l = |I|² 是支路电流平方
%     把 l 也当成变量、并用二阶锥 l ≥ (P²+Q²)/v 把它与 (P,Q) 绑起来，
%     就得到**凸**问题，且 Farivar & Low (2013) 证明对辐射网该松弛是**紧的**。
%
%     所以：LinDistFlow = 忽略线损的线性近似；SOCP = 精确的凸模型。
%     本项目用它量出 LinDistFlow 那个系统性偏差到底有多大。
%
%   变量：P,Q,l (nbr)；v (nb)；g；u (npv)；ls (nb)；sm, sx
%
%   See also distflow_lin, opt_dispatch_lindistflow, compare_distflow_models.

opt = struct('hour', 19, 'gamma', 2000, 'voll', 500, 'curt_pen', 0.01, ...
             'vmax', 1.05, 'vmin', 0.90, 'fix', []);
for k = 1:2:numel(varargin)
    name = varargin{k};
    if ~ischar(name) || ~isfield(opt, name)
        error('opt_dispatch_socp: 未知选项 %s', num2str(name));
    end
    opt.(name) = varargin{k + 1};
end
% fix 非空时进入"潮流模式"：把 g / u / ls 固定成给定值，目标改为极小化 Σl。
% 极小化 Σl 是为了让二阶锥**绷紧**（l 取到它的下界 (P²+Q²)/v），
% 这样得到的电压才等价于一次真实潮流计算。若目标为 0，l 可以取任意大的值，
% 电压就会算错。
powerflow = ~isempty(opt.fix);

nb = L.nb;  nbr = L.nbr;  root = L.root;
npv = numel(dat.pv_bus);
h = opt.hour;
Pd = dat.Pd(:, h);  Qd = dat.Qd;
price = dat.price(h);
v2min = opt.vmin^2;  v2max = opt.vmax^2;

%% ---------- 拓扑定向：让 fdir 指向上游、tdir 指向下游 ----------
% ★ 翻转必须放在 `~seen(w)` 判断之内。若放在外面，一条支路会从两端各被访问一次
%   而被**翻转两次**、方向又反回去 —— 结果是根母线的 children_of 变空、
%   根平衡方程退化成 "-g=0"、电压方程整体错位，最终 coneprog 报不可行。
%   实测踩过：当时 max 残差 -4.7577 恰好等于 -g，就是这条线索定位的。
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
        if ~seen(w)                       % 只在真正"向下走"时才翻转
            if rev, fdir(k) = v;  tdir(k) = w; end
            seen(w) = true;
            parent_of(w) = k;
            queue(end+1) = w;             %#ok<AGROW>
        end
    end
end
assert(all(seen), 'opt_dispatch_socp: 网络不连通');
children_of = cell(nb,1);
for k = 1:nbr, children_of{fdir(k)}(end+1) = k; end
assert(~isempty(children_of{root}), 'opt_dispatch_socp: 根母线没有子支路（拓扑定向出错）');

r = L.r;  x = L.x;

%% ---------- 变量布局 ----------
iP  = 1:nbr;
iQ  = nbr + (1:nbr);
iL  = 2*nbr + (1:nbr);
iV  = 3*nbr + (1:nb);
iG  = 3*nbr + nb + 1;
iU  = iG + (1:npv);
iLS = iG + npv + (1:nb);
iSM = iG + npv + nb + 1;
iSX = iSM + 1;
nvar = iSX;

%% ---------- 目标 ----------
f = zeros(nvar,1);
if powerflow
    f(iL) = 1e-6;                     % 潮流模式：极小化 Σl 让锥绷紧
else
    f(iG)  = price;
    f(iLS) = opt.voll;
    f(iSM) = opt.gamma;
    f(iSX) = opt.gamma;
    f(iU)  = opt.curt_pen;
end

%% ---------- 等式约束 ----------
Aeq = zeros(0, nvar);  beq = zeros(0,1);

% (1) 有功平衡（每个非根母线）
%     P_par - r_par*l_par - Σ_children P_k + (PV出力 - 负荷 + 切负荷) = 0
for j = 1:nb
    if j == root, continue; end
    kp = parent_of(j);
    row = zeros(1,nvar);
    row(iP(kp)) = 1;
    row(iL(kp)) = -r(kp);
    for k = children_of{j}, row(iP(k)) = row(iP(k)) - 1; end
    row(iLS(j)) = 1;                       % 切负荷 = 减少需求
    rhs = Pd(j);
    for k = 1:npv
        if dat.pv_bus(k) == j
            row(iU(k)) = -1;               % 弃光 = 减少出力
            rhs = rhs - pv_avail(k) * dat.pv_cap(k);
        end
    end
    Aeq(end+1,:) = row;  beq(end+1,1) = rhs;                     %#ok<AGROW>
end

% (2) 根母线： Σ_children P_k - g = 0
row = zeros(1,nvar);
for k = children_of{root}, row(iP(k)) = 1; end
row(iG) = -1;
Aeq(end+1,:) = row;  beq(end+1,1) = 0;

% (3) 无功平衡（同结构，只有负荷）
for j = 1:nb
    if j == root, continue; end
    kp = parent_of(j);
    row = zeros(1,nvar);
    row(iQ(kp)) = 1;
    row(iL(kp)) = -x(kp);
    for k = children_of{j}, row(iQ(k)) = row(iQ(k)) - 1; end
    Aeq(end+1,:) = row;  beq(end+1,1) = Qd(j);                   %#ok<AGROW>
end

% (4) 电压降： v_j - v_par + 2(r P_par + x Q_par) - (r²+x²) l_par = 0
for j = 1:nb
    if j == root, continue; end
    kp = parent_of(j);
    row = zeros(1,nvar);
    row(iV(j)) = 1;  row(iV(fdir(kp))) = -1;
    row(iP(kp)) = row(iP(kp)) + 2*r(kp);
    row(iQ(kp)) = row(iQ(kp)) + 2*x(kp);
    row(iL(kp)) = row(iL(kp)) - (r(kp)^2 + x(kp)^2);
    Aeq(end+1,:) = row;  beq(end+1,1) = 0;                       %#ok<AGROW>
end

% (5) 根电压固定 1.0（与 LinDistFlow 模型保持一致，便于公平对比）
row = zeros(1,nvar); row(iV(root)) = 1;
Aeq(end+1,:) = row;  beq(end+1,1) = 1.0;

%% ---------- 电压软约束 ----------
Ain = zeros(0, nvar);  bin = zeros(0,1);
for j = 1:nb
    row = zeros(1,nvar); row(iV(j)) = -1; row(iSM) = -1;   % -v - sm <= -v2min
    Ain(end+1,:) = row;  bin(end+1,1) = -v2min;                  %#ok<AGROW>
    row = zeros(1,nvar); row(iV(j)) = 1;  row(iSX) = -1;   %  v - sx <= v2max
    Ain(end+1,:) = row;  bin(end+1,1) = v2max;                   %#ok<AGROW>
end

%% ---------- 二阶锥： l_k >= (P²+Q²)/v_f ----------
% 等价于 || [2P; 2Q; l - v_f] ||_2 <= l + v_f
% 注意 secondordercone 是函数不是可 .empty 的类，先建第一个再逐个增长。
soc = [];
for k = 1:nbr
    A = zeros(3, nvar);
    A(1, iP(k)) = 2;
    A(2, iQ(k)) = 2;
    A(3, iL(k)) = 1;  A(3, iV(fdir(k))) = -1;
    d = zeros(1, nvar);  d(iL(k)) = 1;  d(iV(fdir(k))) = 1;
    sk = secondordercone(A, zeros(3,1), d, 0);
    if isempty(soc), soc = sk; else, soc(k,1) = sk; end          %#ok<AGROW>
end

%% ---------- 上下界 ----------
% ★ P、Q 必须允许为负：光伏反送时支路潮流方向翻转，P < 0。
%   最初写成 lb = zeros（全部非负），结果渗透率 0 能解、一加光伏就报不可行 ——
%   因为反送时 P 必须是负的。这个 bug 只在"有光伏"时暴露，很容易漏掉。
lb = zeros(nvar,1);
lb(iP) = -inf;
lb(iQ) = -inf;
ub = inf(nvar,1);
ub(iG)  = dat.gmax;
for k = 1:npv, ub(iU(k)) = pv_avail(k) * dat.pv_cap(k); end
for j = 1:nb,  ub(iLS(j)) = Pd(j); end

% 潮流模式：把注入固定成给定值（放在最后，覆盖上面的界）
if powerflow
    lb(iG) = opt.fix.g;   ub(iG) = opt.fix.g;
    for k = 1:npv
        lb(iU(k)) = opt.fix.u(k);   ub(iU(k)) = opt.fix.u(k);
    end
    for j = 1:nb
        lb(iLS(j)) = opt.fix.ls(j);  ub(iLS(j)) = opt.fix.ls(j);
    end
end

%% ---------- 求解 ----------
opts = optimoptions('coneprog', 'Display', 'off', ...
    'OptimalityTolerance', 1e-9, 'ConstraintTolerance', 1e-8);
[x, fval, exitflag, output] = coneprog(f, soc, Ain, bin, Aeq, beq, lb, ub, opts);

% coneprog 失败时可能返回空解，先保护再组装输出
if isempty(x)
    x = nan(nvar, 1);
end
% 调试用：把约束矩阵与拓扑一并返回，便于用"已知可行点"回代检查
DBG = struct('Aeq', Aeq, 'beq', beq, 'Ain', Ain, 'bin', bin, ...
             'idx', struct('P',iP,'Q',iQ,'L',iL,'V',iV,'G',iG,'U',iU,'LS',iLS, ...
                           'SM',iSM,'SX',iSX), ...
             'nvar', nvar, 'fdir', fdir, 'tdir', tdir, 'parent_of', parent_of, ...
             'children_of', {children_of}, 'nb', nb, 'nbr', nbr, 'root', root, ...
             'r', r, 'x', x, 'Pd', Pd, 'Qd', Qd, 'npv', npv, 'pv_bus', dat.pv_bus, ...
             'pv_cap', dat.pv_cap, 'pv_avail', pv_avail);
o = struct('exitflag', exitflag, 'fval', fval, 'output', output, ...
           'idx', struct('P',iP,'Q',iQ,'L',iL,'V',iV,'G',iG,'U',iU,'LS',iLS, ...
                         'SM',iSM,'SX',iSX), ...
           'nvar', nvar, 'model', 'SOCP', 'hour', h, ...
           'V', sqrt(max(x(iV), 0)));
if exitflag ~= 1
    fprintf('  [SOCP] exitflag=%d  msg=%s\n', exitflag, output.message);
end
end
