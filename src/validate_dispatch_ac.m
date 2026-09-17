function Tac = validate_dispatch_ac(mpc0, L, dat, x, out, varargin)
%VALIDATE_DISPATCH_AC  把 LinDistFlow LP 的解回灌 MATPOWER 全交流潮流逐时段校验
%
%   Tac = validate_dispatch_ac(mpc0, L, dat, x, out)
%
%   输入：
%     mpc0  基准算例（case33mg_der 不带光伏出力）
%     L     distflow_lin 返回的模型（用于算出简化模型预测的电压做对比）
%     dat   make_day_data 生成
%     x,out opt_dispatch_lindistflow 的解与索引
%
%   输出 Tac 表，逐时段给出：
%     vmin_ac / vmax_ac      交流潮流算出的真实电压
%     vmin_lin               简化模型预测的电压（同一时段）—— 用于量线性化误差
%     slack_ac / slack_lin   平衡机出力对比（差值 ≈ 网损，因为简化模型无损）
%     shed_ac                交流下"不切负荷会不会越限"的指示
%     load_pct               最大支路负载率
%
%   ★ 为什么要做这一步
%     LinDistFlow 是简化模型：忽略线损、忽略无功的电压耦合细节、平衡母线电压
%     固定为 1.0（真实交流 OPF 会把它顶到上限去降损）。
%     而且它的误差是**有方向的** —— 压降偏小、电压偏高，也就是**对欠压乐观**。
%     所以优化器认为合规的解，在真实交流潮流下可能已经越限。
%     这一步就是把这个 gap 量出来。
%
%   See also opt_dispatch_lindistflow, distflow_eval.

opt = struct('vmax', 1.05, 'vmin', 0.90, 'quiet', false);
for k = 1:2:numel(varargin)
    name = varargin{k};
    if ~isfield(opt, name), error('validate_dispatch_ac: 未知选项 %s', num2str(name)); end
    opt.(name) = varargin{k + 1};
end

define_constants;
mpopt = mpoption('verbose', 0, 'out.all', 0);

nt = dat.nt; nb = dat.nb; npv = out.npv;
ipv  = mpc0.isolar;
iess = mpc0.iess;
ix = out.idx;

Tac = table((1:nt)', nan(nt,1), nan(nt,1), nan(nt,1), nan(nt,1), nan(nt,1), ...
    nan(nt,1), nan(nt,1), nan(nt,1), false(nt,1), ...
    'VariableNames', {'hour','vmin_ac','vmax_ac','vmin_lin','vmax_lin', ...
                      'slack_ac_MW','slack_lin_MW','load_pct','shed_MW','converged'});

for t = 1:nt
    m = mpc0;
    shed = x(ix.LS((t-1)*nb + (1:nb)));
    m.bus(:, PD) = max(dat.Pd(:,t) - shed, 0);      % 应用切负荷
    m.bus(:, QD) = dat.Qd;

    % 光伏实际出力 = 可用上限 - 弃光
    pv_disp = dat.pv_avail_max(t,:)' - x(ix.U((t-1)*npv + (1:npv)));
    m.gen(ipv, PG) = max(pv_disp, 0);

    % 储能净出力（正=放电）
    m.gen(iess, PG) = x(ix.D(t)) - x(ix.C(t));

    r = runpf(m, mpopt);
    Tac.converged(t) = logical(r.success);
    if r.success
        Tac.vmin_ac(t) = min(r.bus(:, VM));
        Tac.vmax_ac(t) = max(r.bus(:, VM));
        Tac.slack_ac_MW(t) = sum(r.gen(1, PG));
        rate = m.branch(1, RATE_A);
        if rate > 0
            Tac.load_pct(t) = max(hypot(r.branch(:,PF), r.branch(:,PT))) / rate * 100;
        end
    end

    % 用**同样的注入**过一遍 LinDistFlow，得到简化模型预测的电压
    p_inj = zeros(nb, 1);
    p_inj(1) = x(ix.G(t));
    p_inj(dat.ess_bus) = p_inj(dat.ess_bus) + x(ix.D(t)) - x(ix.C(t));
    for k = 1:npv
        p_inj(dat.pv_bus(k)) = p_inj(dat.pv_bus(k)) + pv_disp(k);
    end
    p_inj = p_inj - (dat.Pd(:,t) - shed);
    [~, V_lin] = distflow_eval(L, p_inj, -dat.Qd, 1.0);

    Tac.vmin_lin(t)  = min(V_lin);
    Tac.vmax_lin(t)  = max(V_lin);
    Tac.slack_lin_MW(t) = x(ix.G(t));
    Tac.shed_MW(t) = sum(shed);
end

% ---- 摘要 ----
n_over  = sum(Tac.vmax_ac > opt.vmax + 1e-6);
n_under = sum(Tac.vmin_ac < opt.vmin - 1e-6);
dv = Tac.vmin_lin - Tac.vmin_ac;         % 简化模型 - 真实（正 = 简化模型偏乐观）

if opt.quiet, return; end

fprintf('\n=== AC validation of the LinDistFlow dispatch ===\n');
fprintf('converged            : %d / %d\n', sum(Tac.converged), nt);
fprintf('Vmin  AC / Lin       : %.4f / %.4f\n', min(Tac.vmin_ac), min(Tac.vmin_lin));
fprintf('Vmax  AC / Lin       : %.4f / %.4f\n', max(Tac.vmax_ac), max(Tac.vmax_lin));
fprintf('violations           : %d over Vmax, %d under Vmin  (limit %.2f / %.2f)\n', ...
    n_over, n_under, opt.vmax, opt.vmin);
fprintf('slack energy AC/Lin  : %.2f / %.2f MWh  (差 = %.2f MWh = 网损)\n', ...
    sum(Tac.slack_ac_MW), sum(Tac.slack_lin_MW), ...
    sum(Tac.slack_ac_MW) - sum(Tac.slack_lin_MW));
fprintf('max branch loading   : %.1f%%\n', max(Tac.load_pct));
fprintf('\n--- 线性化误差方向（关键）---\n');
fprintf('  mean(Vmin_lin - Vmin_ac) = %+.5f p.u.  (正 = 简化模型高估电压)\n', mean(dv));
fprintf('  max (Vmin_lin - Vmin_ac) = %+.5f p.u.\n', max(dv));
fprintf('  -> LinDistFlow 忽略线损，压降偏小、电压偏高，**对欠压乐观**。\n');
if n_under > 0
    fprintf('  -> 这个乐观方向的偏差确实造成了 %d 个时段在交流下欠压。\n', n_under);
end
fprintf('\n');
end
