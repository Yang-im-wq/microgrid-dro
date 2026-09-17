function plot_overview()
%PLOT_OVERVIEW  项目全景图：把 P0→P5 每一阶段的关键结果拼成一张仪表盘
%
%   跑法：matlab -batch "cd('D:\microgrid-dro\src'); plot_overview"
%   输出：results/figures/OVERVIEW.png
%
%   八个面板，对应五个阶段的每一步：
%     (a) P0 数据桥     —— 光伏概率预测 → 标幺日内曲线
%     (b) P1 微电网算例 —— 33 节点拓扑 + DER 位置 + 潮流方向
%     (c) P2 渗透率扫描 —— 网损与电压，标出"好窗口"
%     (d) P2 网络重构   —— 分支交换收敛到经典最优 139.55 kW
%     (e) P3 经济调度   —— 调度堆叠 + 弃光
%     (f) P4 MOST+AC    —— 24h 时序调度与其 AC 校验揪出的越限
%     (g) P5 DRO        —— 确定性 / 随机 / 分布鲁棒 三方对比
%     (h) 汇总          —— 各约束的触发点 + Vmax 限值敏感性
%
%   除 (b)(f 的电压部分) 需跑潮流外，(f) 的调度与其余全部读已落盘 CSV。

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'cases'));
define_constants;
mpopt = mpoption('verbose', 0, 'out.all', 0);

f = figure('Position', [30 30 1780 900], 'Color', 'w', 'Visible', 'off');
RD = @(nm) readtable(fullfile(root, 'results', 'data', nm), 'Delimiter', ',');

%% ================= (a) P0: 光伏标幺日内曲线 =================
subplot(2,4,1);
Tp = readtable(fullfile(root, 'data', 'pv_daily_profile.csv'), 'Delimiter', ',');
QL = {'q05','q10','q20','q30','q40','q50','q60','q70','q80','q90','q95'};
M = Tp{:, QL};
hh = Tp.hour + Tp.minute/60;
fill([hh; flipud(hh)], [M(:,1); flipud(M(:,11))], [0.12 0.47 0.71], ...
    'FaceAlpha', 0.18, 'EdgeColor', 'none'); hold on;
fill([hh; flipud(hh)], [M(:,3); flipud(M(:,9))], [0.12 0.47 0.71], ...
    'FaceAlpha', 0.30, 'EdgeColor', 'none');
plot(hh, M(:,6), '-', 'Color', [0.12 0.47 0.71], 'LineWidth', 2.0);
yline(1, '--', 'installed capacity', 'Color', [0.84 0.15 0.16], 'FontSize', 7);
xlabel('hour'); ylabel('normalized output [p.u.]');
title('(a) P0 data bridge', 'FontSize', 10);
grid on; box on; xlim([0 24]); ylim([0 1.05]);
text(1, 0.88, sprintf('q50 peak %.2f\nband q05-q95', max(M(:,6))), 'FontSize', 7, 'Color', [0.35 0.35 0.35]);

%% ================= (b) P1: 网络拓扑 =================
subplot(2,4,2);
m = case33mg_der('pv_penetration', 1.5);
r = runpf(m, mpopt);
draw_topo(m, r, layout33());
title('(b) P1 microgrid case', 'FontSize', 10);

%% ================= (c) P2: 渗透率扫描 =================
subplot(2,4,3);
S5 = RD('sweep_pv_vmax105.csv');
S1 = RD('sweep_pv_vmax110.csv');
xx = S5.pen * 100;
yyaxis left
plot(xx, S5.loss_kW, '-', 'LineWidth', 2.0, 'Color', [0.12 0.47 0.71]);
ylabel('Loss [kW]', 'FontSize', 8);
yyaxis right
plot(xx, S5.vmax_pu, '-', 'LineWidth', 1.8, 'Color', [0.84 0.15 0.16]); hold on;
plot(S1.pen*100, S1.vmax_pu, '--', 'LineWidth', 1.2, 'Color', [0.95 0.55 0.55]);
yline(1.05, ':', 'Color', [0.84 0.15 0.16], 'FontSize', 7);
yline(1.10, ':', 'Color', [0.95 0.55 0.55], 'FontSize', 7);
ylabel('Bus Vmax [p.u.]', 'FontSize', 8);
xlabel('PV penetration [% of load]', 'FontSize', 8);
title('(c) P2 penetration sweep', 'FontSize', 10);
grid on; box on; xlim([0 round(max(xx))]);

%% ================= (d) P2: 重构收敛 =================
subplot(2,4,4);
Rr = RD('reconfig_radial_pv000_vmax105.csv');
plot(Rr.step, Rr.loss_kW, 'o-', 'LineWidth', 1.8, 'MarkerSize', 5, ...
    'Color', [0.30 0.55 0.20], 'MarkerFaceColor', [0.30 0.55 0.20]);
yline(Rr.loss_kW(1), '--', sprintf('base %.1f', Rr.loss_kW(1)), 'Color', [0.5 0.5 0.5], 'FontSize', 7);
text(Rr.step(end), Rr.loss_kW(end)+5, sprintf('%.2f kW\nclassic optimum', Rr.loss_kW(end)), ...
    'FontSize', 7, 'HorizontalAlignment', 'right');
xlabel('branch-exchange iteration', 'FontSize', 8); ylabel('Loss [kW]', 'FontSize', 8);
title('(d) P2 reconfiguration', 'FontSize', 10);
grid on; box on;

%% ================= (e) P3: 经济调度 =================
subplot(2,4,5);
E3 = RD('ed_sweep_pen200_vmax105.csv');
area(E3.avail, [E3.slack_MW, E3.pv_disp_MW], 'LineWidth', 0.8); hold on;
plot(E3.avail, E3.pv_avail_MW, 'k--', 'LineWidth', 1.4);
fill([E3.avail; flipud(E3.avail)], [E3.pv_avail_MW; flipud(E3.pv_disp_MW + E3.slack_MW)], ...
    [0.85 0.33 0.10], 'FaceAlpha', 0.28, 'EdgeColor', 'none');
xlabel('PV availability [p.u.]', 'FontSize', 8); ylabel('Power [MW]', 'FontSize', 8);
title('(e) P3 AC OPF dispatch', 'FontSize', 10);
legend({'slack','PV','available','curtailed'}, 'Location', 'northwest', 'FontSize', 6);
grid on; box on; xlim([0 1]);

%% ================= (f) P4: MOST 24h + AC 校验 =================
subplot(2,4,6);
M4 = RD('most_24h_pen150_vmax105.csv');
hr = M4.hour;
plot(hr, M4.vmax_pu, '-o', 'LineWidth', 2, 'MarkerSize', 3, 'Color', [0.84 0.15 0.16]); hold on;
plot(hr, M4.vmin_pu, '-o', 'LineWidth', 2, 'MarkerSize', 3, 'Color', [0.12 0.47 0.71]);
yline(1.05, '--', 'Vmax 1.05', 'Color', [0.84 0.15 0.16], 'FontSize', 7);
yline(0.90, '--', 'Vmin 0.90', 'Color', [0.12 0.47 0.71], 'FontSize', 7);
bad = hr(M4.vmin_pu < 0.90 - 1e-6);
if ~isempty(bad)
    plot(bad, M4.vmin_pu(M4.vmin_pu < 0.90 - 1e-6), 'v', 'MarkerSize', 7, ...
        'MarkerFaceColor', [0.95 0.75 0.10], 'MarkerEdgeColor', 'k');
end
xlabel('hour', 'FontSize', 8); ylabel('Bus voltage [p.u.]', 'FontSize', 8);
title(sprintf('(f) P4 MOST + AC check: %d violations', numel(bad)), 'FontSize', 10);
grid on; box on; xlim([1 24]);

%% ================= (g) P5: DET vs SAA vs DRO =================
subplot(2,4,7);
D5 = RD('dro_compare_pen150_eps050_g50000hi.csv');
Mx = [D5.E_cost, D5.CVaR95, D5.worst_cost];
b = bar(Mx, 'grouped'); hold on;
set(gca, 'XTickLabel', D5.plan, 'FontSize', 8);
ylabel('Cost [$/day]', 'FontSize', 8);
title('(g) P5 DET vs SAA vs DRO', 'FontSize', 10);
legend({'E[cost]','CVaR95','worst'}, 'Location', 'northwest', 'FontSize', 6);
grid on; box on;
yl = [min(Mx(:))*0.97, max(Mx(:))*1.06];
ylim(yl);
text(2, yl(2) - 0.02*(yl(2)-yl(1)), sprintf('spread across plans < 0.15%%  ->  DRO buys little here'), ...
    'FontSize', 6.5, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', 'Color', [0.3 0.3 0.3]);

%% ================= (h) 汇总 =================
subplot(2,4,8);
ev = [56.5, 80.0, 102.7, 171.8, 244.9, 237.7];
nm = {'Vmin \geq 0.95 cured', 'min loss', 'slack reverses', ...
      'Vmax 1.05 exceeded', 'Vmax 1.10 exceeded', 'thermal congestion'};
col = [0.30 0.55 0.20; 0.30 0.55 0.20; 0.45 0.45 0.45; ...
       0.84 0.15 0.16; 0.95 0.60 0.60; 0.55 0.35 0.70];
hold on;
for k = 1:numel(ev)
    barh(k, ev(k), 0.62, 'FaceColor', col(k,:), 'EdgeColor', 'none');
    text(ev(k)+4, k, sprintf('%.0f%%', ev(k)), 'FontSize', 7, 'VerticalAlignment', 'middle');
end
fill([56.5 171.8 171.8 56.5], [0.35 0.35 numel(ev)+0.65 numel(ev)+0.65], ...
    [0.30 0.55 0.20], 'FaceAlpha', 0.10, 'EdgeColor', 'none');
set(gca, 'YTick', 1:numel(ev), 'YTickLabel', nm, 'YDir', 'reverse', 'FontSize', 7);
xlabel('PV penetration [% of load]', 'FontSize', 8);
xlim([0 275]);
title('(h) hosting-capacity limits', 'FontSize', 10);
grid on; box on;
text(114, 0.7, 'good window', 'FontSize', 7, 'Color', [0.2 0.45 0.2], 'HorizontalAlignment', 'center');

sgtitle({'IEEE 33-bus microgrid optimization platform : P0 \rightarrow P5 walkthrough', ...
    'MATPOWER 8.1 + PV probabilistic forecast  |  Vmax = 1.05 (ANSI C84.1 Range A)  |  every layer AC-validated'}, ...
    'FontSize', 12);

out = fullfile(root, 'results', 'figures', 'OVERVIEW.png');
exportgraphics(f, out, 'Resolution', 130);
close(f);
fprintf('saved: %s\n', out);
end


%% ---------------- helpers ----------------
function draw_topo(m, r, XY)
define_constants;
hold on;
Smax = max(hypot(r.branch(:,PF), r.branch(:,PT)));
for k = find(m.branch(:, BR_STATUS) == 0)'
    plot(XY([m.branch(k,F_BUS) m.branch(k,T_BUS)],1), ...
         XY([m.branch(k,F_BUS) m.branch(k,T_BUS)],2), ...
         '--', 'Color', [0.80 0.80 0.80], 'LineWidth', 0.6);
end
for k = find(m.branch(:, BR_STATUS) == 1)'
    S = hypot(r.branch(k,PF), r.branch(k,PT));
    if r.branch(k, PF) >= 0, c = [0.15 0.35 0.75]; else, c = [0.82 0.15 0.12]; end
    plot(XY([m.branch(k,F_BUS) m.branch(k,T_BUS)],1), ...
         XY([m.branch(k,F_BUS) m.branch(k,T_BUS)],2), '-', ...
         'Color', c, 'LineWidth', 0.5 + 3.0*S/max(Smax,eps));
end
vm = r.bus(:, VM);
cmap = jet(256);
lo = 0.90; hi = 1.05;
ci = max(1, min(256, round((vm-lo)/(hi-lo)*255)+1));
scatter(XY(:,1), XY(:,2), 32, cmap(ci,:), 'filled', ...
    'MarkerEdgeColor', [0.25 0.25 0.25], 'LineWidth', 0.3);
for b = [18 22 25 33]
    plot(XY(b,1), XY(b,2), 's', 'MarkerSize', 8, 'LineWidth', 1.6, ...
        'MarkerEdgeColor', [0.05 0.55 0.10], 'MarkerFaceColor', 'none');
end
plot(XY(18,1), XY(18,2), 'o', 'MarkerSize', 5, 'LineWidth', 1.4, ...
    'MarkerEdgeColor', [0.05 0.55 0.75], 'MarkerFaceColor', 'none');
plot(XY(1,1), XY(1,2), 'p', 'MarkerSize', 14, 'LineWidth', 1.1, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.95 0.85 0.20]);
axis equal; axis off; xlim([0 20]); ylim([-9.5 1.5]); colormap(cmap);
end


function XY = layout33()
XY = nan(33,2);
for i = 1:18, XY(i,:) = [i, 0]; end
for k = 1:4, XY(18+k,:) = [2,   -k]; end
for k = 1:3, XY(22+k,:) = [4.2, -k]; end
for k = 1:8, XY(25+k,:) = [15,  -k]; end
end
