function [x, out] = opt_dispatch_lindistflow(L, dat, pv_avail, varargin)
%OPT_DISPATCH_LINDISTFLOW  基于 LinDistFlow 的日前调度 LP（单场景）
%
%   [x, out] = opt_dispatch_lindistflow(L, dat, pv_avail)
%   [x, out] = opt_dispatch_lindistflow(L, dat, pv_avail, 'gamma', 2000, 'voll', 500)
%
%   这是在**简化模型**（LinDistFlow）上做优化的核心。之后 DRO / 随机 / 确定性
%   都调用它，最后统一用 MATPOWER 全交流潮流校验。
%
%   输入：
%     L         distflow_lin 返回的模型
%     dat       make_day_data 生成的数据结构
%     pv_avail  (nt x npv) 各时段各光伏节点的**可用出力系数** 0~1（不确定量 ξ）
%
%   选项：
%     'gamma'    电压越限惩罚 $/p.u.^2·h (默认 2000)
%     'voll'     切负荷惩罚 $/MWh        (默认 500)
%     'curt_pen' 弃光惩罚 $/MWh          (默认 0.01，仅为消除退化)
%     'vmax'/'vmin'                      (默认 1.05 / 0.90)
%
%   决策变量：
%     d(nt)      储能放电 MW            c(nt)      储能充电 MW
%     e(nt)      SOC MWh                g(nt)      变电站进口 MW
%     u(npv,nt)  各光伏节点弃光 MW
%     ls(nb,nt)  各节点切负荷 MW（<= 该节点负荷，故平衡母线自动为 0）
%     sm(nb,nt)  电压下限越限量 p.u.^2  sx(nb,nt)  电压上限越限量 p.u.^2
%
%   ★ 为什么必须有切负荷
%     本算例晚高峰负荷 5.34 MW，而变电站容量上限 gmax + 储能功率 = 5.2 MW，
%     **顶不上**。没有切负荷的话 LP 直接不可行（实测 exitflag = -2）。
%     切负荷以 VOLL 计价，正好成为"不确定性代价"的载体 ——
%     最坏场景下要切多少负荷，取决于日前把储能调度得保不保守。
%     这就是 DRO 的意义所在。
%
%   ★ 关于电压约束的建模选择
%     v_root 固定为 1.0（与 case33mg 平衡母线一致），即**不把平衡母线电压
%     当作调节手段**。交流 OPF 里它会顶到 1.05 去降损，这里没有这个自由度。
%     这是一个有意的简化，代价会在最后 AC 校验里体现。
%
%   See also distflow_lin, distflow_eval, make_day_data, dro_24h.

%% ---------- 选项 ----------
opt = struct('gamma', 2000, 'voll', 500, 'curt_pen', 0.01, ...
             'vmax', 1.05, 'vmin', 0.90);
for k = 1:2:numel(varargin)
    name = varargin{k};
    if ~ischar(name) || ~isfield(opt, name)
        error('opt_dispatch_lindistflow: 未知选项 %s', num2str(name));
    end
    opt.(name) = varargin{k + 1};
end

nt = dat.nt;  nb = dat.nb;  npv = numel(dat.pv_bus);
v2min = opt.vmin^2;  v2max = opt.vmax^2;

%% ---------- 电压对注入的线性映射 ----------
Mv  = 2 * L.A_p * diag(L.r) * L.A_p';       % v = v_root - dvQ + Mv * p_inj
Q_br0 = -L.A_p' * (-dat.Qd);
dvQ   = 2 * L.A_p * (L.x .* Q_br0);         % 无功造成的固定压降 (nb x 1)

B_slack = zeros(nb,1); B_slack(1) = 1;
B_ess   = zeros(nb,1); B_ess(dat.ess_bus) = 1;
B_pv    = zeros(nb, npv);
for k = 1:npv, B_pv(dat.pv_bus(k), k) = 1; end

%% ---------- 变量布局 ----------
iD  = 1:nt;
iC  = nt + (1:nt);
iE  = 2*nt + (1:nt);
iG  = 3*nt + (1:nt);
iU  = 4*nt + (1:(npv*nt));
iLS = 4*nt + npv*nt + (1:(nb*nt));
iSM = 4*nt + npv*nt + nb*nt + (1:(nb*nt));
iSX = 4*nt + npv*nt + 2*nb*nt + (1:(nb*nt));
nvar = 4*nt + npv*nt + 3*nb*nt;

iUf  = @(k,t) iU( (t-1)*npv + k );
iLSf = @(j,t) iLS( (t-1)*nb + j );
iSMf = @(j,t) iSM( (t-1)*nb + j );
iSXf = @(j,t) iSX( (t-1)*nb + j );

%% ---------- 目标函数 ----------
f = zeros(nvar, 1);
f(iG)  = dat.price(:);
f(iLS) = opt.voll;
f(iSM) = opt.gamma;
f(iSX) = opt.gamma;
f(iU)  = opt.curt_pen;

%% ---------- 等式约束 ----------
Aeq = []; beq = [];
eta = dat.ess_eta;

% (1) SOC 动态（日循环：末时段 SOC 回到初值）
for t = 1:nt
    row = zeros(1, nvar);
    row(iE(t)) = 1;
    if t > 1, row(iE(t-1)) = -1; end
    row(iC(t)) = -eta;
    row(iD(t)) =  1/eta;
    Aeq = [Aeq; row];                                            %#ok<AGROW>
    beq = [beq; (t == 1) * dat.ess_e0];                          %#ok<AGROW>
end
row = zeros(1, nvar); row(iE(nt)) = 1;
Aeq = [Aeq; row];  beq = [beq; dat.ess_e0];

% (2) 功率平衡（无损 LinDistFlow：全网注入之和为 0）
%     物理形式： g + pv_delivered + d - c + shed = load
%               其中 pv_delivered = pmax - u
%     ->  g - u + d - c + shed = load - pmax
%     ★ 符号极易搞错。实测踩过：把 u 的符号写成 +、pmax 放到 RHS 正号、
%       充放电也反了，结果是"光伏反而抬高所需平衡机出力"，LP 把光伏全弃了
%       （curtail = 100% 可用量），而且照样 exitflag=1 报"求解成功"。
for t = 1:nt
    row = zeros(1, nvar);
    row(iG(t)) =  1;
    row(iD(t)) =  1;        % 放电 = 电源
    row(iC(t)) = -1;        % 充电 = 负荷
    rhs = sum(dat.Pd(:,t));
    for k = 1:npv
        row(iUf(k,t)) = -1;                         % 弃光 = 减少光伏出力
        rhs = rhs - pv_avail(t,k) * dat.pv_cap(k);
    end
    for j = 1:nb
        row(iLSf(j,t)) = 1;                         % 切负荷 = 减少需求
    end
    Aeq = [Aeq; row];                                            %#ok<AGROW>
    beq = [beq; rhs];                                            %#ok<AGROW>
end

%% ---------- 不等式约束（电压软约束）----------
% p_inj(t) = B_slack*g + B_ess*(d-c) + B_pv*(pmax - u) - (Pd - ls)
v_root = 1.0;
Aineq = []; bineq = [];
for t = 1:nt
    p_fix = -dat.Pd(:,t);
    for k = 1:npv
        p_fix = p_fix + B_pv(:,k) * pv_avail(t,k) * dat.pv_cap(k);
    end
    vfix = v_root - dvQ + Mv * p_fix;

    cg  = Mv * B_slack;
    ce  = Mv * B_ess;
    cpv = Mv * B_pv;

    for j = 1:nb
        % --- 下限: -(Mv*p_dec)_j - sm_j <= vfix_j - v2min
        row = zeros(1, nvar);
        row(iG(t)) = -cg(j);
        row(iD(t)) = -ce(j);
        row(iC(t)) =  ce(j);
        for k = 1:npv, row(iUf(k,t)) =  cpv(j,k); end
        for i = 1:nb, row(iLSf(i,t)) = -Mv(j,i); end   % 母线 i 切负荷抬高母线 j 电压
        row(iSMf(j,t)) = -1;
        Aineq = [Aineq; row];                                        %#ok<AGROW>
        bineq = [bineq; vfix(j) - v2min];                            %#ok<AGROW>

        % --- 上限: (Mv*p_dec)_j - sx_j <= v2max - vfix_j
        row = zeros(1, nvar);
        row(iG(t)) =  cg(j);
        row(iD(t)) =  ce(j);
        row(iC(t)) = -ce(j);
        for k = 1:npv, row(iUf(k,t)) = -cpv(j,k); end
        for i = 1:nb, row(iLSf(i,t)) =  Mv(j,i); end
        row(iSXf(j,t)) = -1;
        Aineq = [Aineq; row];                                        %#ok<AGROW>
        bineq = [bineq; v2max - vfix(j)];                            %#ok<AGROW>
    end
end

%% ---------- 变量上下界 ----------
lb = zeros(nvar, 1);
ub =  inf(nvar, 1);
ub(iD) = dat.ess_pmax;
ub(iC) = dat.ess_pmax;
ub(iE) = dat.ess_emax;   lb(iE) = dat.ess_emin;
ub(iG) = dat.gmax;
for t = 1:nt
    for k = 1:npv
        ub(iUf(k,t)) = pv_avail(t,k) * dat.pv_cap(k);
    end
    for j = 1:nb
        ub(iLSf(j,t)) = dat.Pd(j,t);        % 最多把该节点负荷切光（母线1无负荷 -> 自动为0）
    end
end

%% ---------- 求解 ----------
lpopts = optimoptions('linprog', 'Display', 'off', 'OptimalityTolerance', 1e-9);
[x, fval, exitflag, output] = linprog(f, Aineq, bineq, Aeq, beq, lb, ub, lpopts);

out = struct('exitflag', exitflag, 'fval', fval, 'output', output, 'nvar', nvar, ...
             'idx', struct('D',iD,'C',iC,'E',iE,'G',iG,'U',iU,'LS',iLS,'SM',iSM,'SX',iSX), ...
             'nt', nt, 'nb', nb, 'npv', npv);
if exitflag ~= 1
    warning('opt_dispatch_lindistflow: linprog exitflag = %d (%s)', exitflag, output.message);
end
end
