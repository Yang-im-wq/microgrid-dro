%% demo_matpower.m —— MATPOWER 快速体验（跑一遍看它能干什么）
%
% 用法（二选一）：
%   MATLAB 里： cd D:\microgrid-dro\src  →  输入 demo_matpower
%   命令行   ： matlab -batch "cd('D:\microgrid-dro\src'); demo_matpower"
%
% 输出全部用 ASCII 标签，避免 Windows 控制台 GBK 把中文打乱。
%
% 单位提醒：case33mg 的 baseMVA = 1，所以 gen/bus 里的功率数值单位就是 MW。

function demo_matpower()

fprintf('\n');
fprintf('############################################################\n');
fprintf('#  MATPOWER 8.1 快速体验                                     #\n');
fprintf('############################################################\n\n');

define_constants;

% 关掉 MATPOWER 的详细打印，否则结果被 printpf 淹没
mpopt = mpoption('verbose', 0, 'out.all', 0);

v = mpver('all');
fprintf('[version] MATPOWER %s (%s)\n', v.Version, v.Date);
fprintf('[path]    %s\n\n', which('runpf'));

%% ==================== 1. 潮流 ====================
fprintf('============================================================\n');
fprintf('1. runpf -- power flow, IEEE 9 bus\n');
fprintf('============================================================\n');
t0 = tic;
r = runpf('case9', mpopt);
fprintf('converged=%d  iter=%d  time=%.3fs\n', r.success, r.iterations, toc(t0));
fprintf('load=%.1f MW  loss=%.2f MW  Vmin=%.4f pu @bus%d\n\n', ...
    sum(r.bus(:, PD)), sum(r.branch(:, PF) + r.branch(:, PT)), ...
    min(r.bus(:, VM)), find(r.bus(:, VM) == min(r.bus(:, VM)), 1));

%% ==================== 2. 最优潮流 ====================
fprintf('============================================================\n');
fprintf('2. runopf -- optimal power flow, IEEE 9 bus\n');
fprintf('============================================================\n');
o = runopf('case9', mpopt);
fprintf('converged=%d   cost=%.2f $/h\n', o.success, o.f);
fprintf('gen Pg (MW):');
fprintf(' %.2f', o.gen(:, PG));
fprintf('\n(official benchmark = 5296.69 $/h)\n\n');

%% ==================== 3. IEEE 33 节点配电网 ====================
fprintf('============================================================\n');
fprintf('3. IEEE 33-bus distribution (our microgrid base case)\n');
fprintf('============================================================\n');
for c = {'case33bw', 'case33mg'}
    nm = c{1};
    m = feval(nm);
    t0 = tic;
    r = runpf(nm, mpopt);
    fprintf('%-9s baseMVA=%-4g PF %.0f iter %.2fs | loss=%6.1f kW | Vmin=%.4f @bus%d\n', ...
        nm, m.baseMVA, r.iterations, toc(t0), ...
        sum(r.branch(:, PF) + r.branch(:, PT)) * 1000, ...
        min(r.bus(:, VM)), find(r.bus(:, VM) == min(r.bus(:, VM)), 1));
end
fprintf('(Vmin=0.9131 @bus18 is the IEEE 33-bus literature value)\n\n');

%% ==================== 4. 联络开关闭合 ====================
fprintf('============================================================\n');
fprintf('4. close 5 tie switches: radial feeder -> meshed microgrid\n');
fprintf('============================================================\n');
m = case33mg;
n_tie = sum(m.branch(:, BR_STATUS) == 0);
r_rad = runpf(m, mpopt);
L_rad = sum(r_rad.branch(:, PF) + r_rad.branch(:, PT));
fprintf('radial : loss=%6.1f kW   Vmin=%.4f pu\n', L_rad * 1000, min(r_rad.bus(:, VM)));
m.branch(:, BR_STATUS) = 1;
r_mesh = runpf(m, mpopt);
L_mesh = sum(r_mesh.branch(:, PF) + r_mesh.branch(:, PT));
fprintf('meshed : loss=%6.1f kW   Vmin=%.4f pu   (%d tie switches closed)\n', ...
    L_mesh * 1000, min(r_mesh.bus(:, VM)), n_tie);
fprintf('-> loss -%.1f%%   Vmin +%.4f pu\n\n', 100 * (1 - L_mesh / L_rad), ...
    min(r_mesh.bus(:, VM)) - min(r_rad.bus(:, VM)));

%% ==================== 5. 光伏渗透率扫描 ====================
fprintf('============================================================\n');
fprintf('5. add PV at bus18 (weakest end node), sweep penetration\n');
fprintf('============================================================\n');
m0 = case33mg;
load_MW = sum(m0.bus(:, PD));
ncol = size(m0.gen, 2);
pens = [0 0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0];
R = nan(numel(pens), 6);            % pv, slack[kW], loss[kW], vmin, vmax, vmax_bus
for k = 1:numel(pens)
    pen = pens(k);
    m = case33mg;
    pv = pen * load_MW;                              % MW
    if pen > 0
        m.gen(end + 1, 1:ncol) = zeros(1, ncol);
        g = size(m.gen, 1);
        m.gen(g, [GEN_BUS GEN_STATUS VG MBASE]) = [18 1 1 100];
        m.gen(g, PG) = pv;      % <-- 关键: PG 才是潮流用的实际注入, PMAX 只是 OPF 上界
        m.gen(g, [PMAX PMIN]) = [pv pv];
        m.gen(g, [QMAX QMIN]) = [0 0];               % unity power factor
        m.gencost(end + 1, 1:6) = [POLYNOMIAL 0 0 2 0 0];
        % bus 18 保持 PQ 类型: gen 在 PQ 母线上按"定功率注入"处理
        % (MATPOWER 的 makeSbus 对所有母线累加 gen 注入, 与 BUS_TYPE 无关)
        % -> 正是单位功率因数光伏的行为: 不参与调压, 电压让它自然顶上去
    end
    r = runpf(m, mpopt);
    if r.success
        [vmax, ib] = max(r.bus(:, VM));
        R(k, :) = [pv, r.gen(1, PG) * 1000, ...
            sum(r.branch(:, PF) + r.branch(:, PT)) * 1000, ...
            min(r.bus(:, VM)), vmax, ib];
    end
end
% 最小网损点（只在 PV>0 里找）
loss_tmp = R(:, 3);
loss_tmp(R(:, 1) == 0) = Inf;
[~, k_minloss] = min(loss_tmp);

fprintf('total load = %.3f MW   |   PV = unity-PF injection at bus 18 (end node)\n', load_MW);
fprintf('%-7s %-9s %-11s %-11s %-11s %-13s %s\n', ...
    'PV[%]', 'PV[MW]', 'Slack[kW]', 'loss[kW]', 'Vmin[pu]', 'Vmax[pu]', 'note');
fprintf('%s\n', repmat('-', 1, 74));
for k = 1:numel(pens)
    if isnan(R(k, 1))
        fprintf('%-7s not converged\n', sprintf('%.0f', pens(k) * 100));
        continue
    end
    note = '';
    if R(k, 5) > 1.1
        note = sprintf('OVER-VOLTAGE @bus%d', R(k, 6));
    elseif k == k_minloss
        note = 'min loss';
    end
    fprintf('%-7s %-9.3f %-11.1f %-11.1f %-11.4f %-13s %s\n', ...
        sprintf('%.0f', pens(k) * 100), R(k, 1), R(k, 2), R(k, 3), R(k, 4), ...
        sprintf('%.4f@%d', R(k, 5), R(k, 6)), note);
end
fprintf('\nphysics: PV loss first DROPS (local supply) then RISES (reverse flow),\n');
fprintf('         and at high penetration the end-node voltage crosses the 1.1 pu limit.\n\n');

fprintf('############################################################\n');
fprintf('#  try it yourself:                                        #\n');
fprintf('#    runpf(''case33mg'')        %% power flow                 #\n');
fprintf('#    runopf(''case33mg'')        %% optimal power flow         #\n');
fprintf('#    help mpoption       %% all tunable options            #\n');
fprintf('#    test_matpower       %% run official test suite (slow)  #\n');
fprintf('############################################################\n');
end
