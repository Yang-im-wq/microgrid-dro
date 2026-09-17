function test_p5_distflow()
%TEST_P5_DISTFLOW  P5 验收：LinDistFlow 精度对标 MATPOWER 全交流潮流
%
%   跑法：matlab -batch "cd('D:\microgrid-dro\src'); test_p5_distflow"
%
%   ★ 这一步是 P5 的地基
%     后面所有 DRO / 随机 / 确定性的优化都在 LinDistFlow 这个**简化模型**上做。
%     如果它跟真实交流潮流差得远，优化出来的结论就没有意义。
%     所以先用 MATPOWER 全交流潮流当"真值"，量出线性化误差到底多大。
%
%   做法：跑 AC 潮流拿到收敛解 → 从解里取出各节点真实注入 → 把**同样的注入**
%         喂给 LinDistFlow → 比较两者算出的电压。
%         两者之差**纯粹**是线性化误差（注入完全相同，不是模型差异）。
%
%   输出全 ASCII，避免 Windows 控制台 GBK 乱码。

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'src'), fullfile(root, 'cases'));
define_constants;
mpopt = mpoption('verbose', 0, 'out.all', 0);

fprintf('\n############################################################\n');
fprintf('#  P5 : LinDistFlow accuracy vs full AC power flow          #\n');
fprintf('############################################################\n\n');

L = distflow_lin(case33mg_der('pv_penetration', 0, 'storage', false));
fprintf('network: %d buses, %d branches, radial-connected = %d\n\n', ...
    L.nb, L.nbr, L.connected);

% 工况：覆盖"轻载/基态/重载/高光伏"几个电压偏离最大的场景
cfgs = {
    'base                      ', 0,    false, 1.00, 0;
    'light load    (x0.5)      ', 0,    false, 0.50, 0;
    'heavy load    (x1.4)      ', 0,    false, 1.40, 0;
    'PV 80%  (min loss)        ', 0.8,  false, 1.00, 0;
    'PV 150%                   ', 1.5,  false, 1.00, 0;
    'PV 200% + heavy load      ', 2.0,  false, 1.20, 0;
    'PV 100% + storage charge  ', 1.0,  true,  1.00, -1.0;
};

fprintf('%-26s | %8s | %8s | %8s | %8s\n', ...
    'case', 'max|dV|', 'mean|dV|', 'Vmin_AC', 'Vmin_Lin');
fprintf('%s\n', repmat('-', 1, 74));

worst = 0;
for k = 1:size(cfgs, 1)
    m = case33mg_der('pv_penetration', cfgs{k,2}, 'storage', cfgs{k,3}, ...
                     'load_scale', cfgs{k,4});
    if cfgs{k,5} ~= 0 && isfield(m, 'iess') && ~isempty(m.iess)
        m.gen(m.iess, PG) = cfgs{k,5};       % 储能定功率（充电为负）
    end
    r = runpf(m, mpopt);
    if ~r.success
        fprintf('%-26s | NOT CONVERGED\n', cfgs{k,1});
        continue
    end

    % ---- 从 AC 解里取真实注入（发电 - 负荷）----
    nb = size(r.bus, 1);
    p_gen = accumarray(r.gen(:, GEN_BUS), r.gen(:, PG), [nb 1]);
    q_gen = accumarray(r.gen(:, GEN_BUS), r.gen(:, QG), [nb 1]);
    p_inj = p_gen - r.bus(:, PD);
    q_inj = q_gen - r.bus(:, QD);

    % ---- 同样的注入喂给 LinDistFlow ----
    v_root = r.bus(L.root, VM) ^ 2;
    [~, V_lin] = distflow_eval(L, p_inj, q_inj, v_root);
    V_ac = r.bus(:, VM);

    d = abs(V_lin - V_ac);
    fprintf('%-26s | %8.5f | %8.5f | %8.4f | %8.4f\n', ...
        cfgs{k,1}, max(d), mean(d), min(V_ac), min(V_lin));
    worst = max(worst, max(d));
end

fprintf('\n--- summary ---\n');
fprintf('worst voltage error over all cases = %.5f p.u. (%.3f%%)\n', worst, worst * 100);
if worst < 0.005
    fprintf('verdict: EXCELLENT  (< 0.5%%)\n');
elseif worst < 0.010
    fprintf('verdict: GOOD       (< 1.0%%)\n');
elseif worst < 0.020
    fprintf('verdict: ACCEPTABLE (< 2.0%%) -- LinDistFlow 的典型精度，可用于优化\n');
else
    fprintf('verdict: POOR       -- 结论不可信，需要换模型\n');
end

% ★ 关键：误差的方向
%   LinDistFlow 忽略了线损，所以它算出的压降**偏小**、电压**偏高**。
%   也就是它对"欠压"是**乐观的** —— 优化器可能以为电压合规，实际交流潮流下已越限。
%   这是必须靠 AC 校验兜住的方向性偏差。
fprintf('\nnote: LinDistFlow 忽略线损 -> 压降偏小、电压偏高（对欠压乐观）。\n');
fprintf('      该方向性偏差必须靠最后一步 AC 校验兜住。\n');

% ---- 支路潮流：LinDistFlow 忽略线损，所以会比 AC 小一点 ----
m = case33mg_der('pv_penetration', 1.0, 'storage', false);
r = runpf(m, mpopt);
nb = size(r.bus,1);
p_inj = accumarray(r.gen(:,GEN_BUS), r.gen(:,PG), [nb 1]) - r.bus(:,PD);
P_br_lin = -L.A_p' * p_inj;          % 注意负号：见 distflow_eval 的符号约定说明
P_br_ac  = r.branch(L.ci, PF);
loss_ac  = sum(r.branch(:, PF) + r.branch(:, PT));
fprintf('\nbranch flow check:\n');
fprintf('  LinDistFlow 忽略线损，所以它算的支路潮流应**略小于** AC：\n');
fprintf('    max |dP| = %.4f MW   (AC 总网损 = %.4f MW，量级吻合即正常)\n', ...
    max(abs(P_br_lin - P_br_ac)), loss_ac);
fprintf('    max P_lin = %.4f MW ; max P_ac = %.4f MW\n', ...
    max(P_br_lin), max(P_br_ac));
fprintf('\n');
end
