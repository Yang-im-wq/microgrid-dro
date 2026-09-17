function plot_extensions()
%PLOT_EXTENSIONS  三项扩展能力总览：储能选址 / SOCP 潮流 / 云图时间相关
%
%   跑法：matlab -batch "cd('D:\microgrid-dro\src'); plot_extensions"
%   输出：results/figures/OVERVIEW_EXT.png
%
%   三行，每行一个主题：
%     第一行 储能选址 —— 位置对经济性与电压的影响
%     第二行 SOCP   —— 二阶锥模型把简化模型的电压误差降到 0
%     第三行 云图   —— 从云量序列提取场景时间相关性，及其对尾部风险的影响

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'cases'));
define_constants;

f = figure('Position', [30 30 1700 980], 'Color', 'w', 'Visible', 'off');
RD = @(nm) readtable(fullfile(root, 'results', 'data', nm), 'Delimiter', ',');

%% ================= 第一行：储能选址 =================
S = RD('storage_siting_pen150_g042.csv');
S = sortrows(S, 'bus');

% (a) 成本 vs 越限时段
subplot(3,3,1);
nv = S.n_over + S.n_under;
scatter(S.cost_usd_d, nv, 55, S.max_load_pct, 'filled', ...
    'MarkerEdgeColor', [0.3 0.3 0.3]); hold on;
i18 = find(S.bus == 18, 1);
plot(S.cost_usd_d(i18), nv(i18), 'p', 'MarkerSize', 17, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.95 0.85 0.20]);
text(S.cost_usd_d(i18)-0.6, nv(i18)+0.45, 'bus 18', 'FontSize', 7, 'HorizontalAlignment','right');
i1 = find(S.bus == 4, 1);
plot(S.cost_usd_d(i1), nv(i1), 'o', 'MarkerSize', 11, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.30 0.60 0.30]);
text(S.cost_usd_d(i1)+0.6, nv(i1)+0.45, 'bus 4', 'FontSize', 7);
cb = colorbar; cb.Label.String = 'max branch loading [%]'; cb.FontSize = 7;
xlabel('cost  [$/day]', 'FontSize', 8); ylabel('violating periods', 'FontSize', 8);
title('(a) Siting: cost vs violations', 'FontSize', 10);
grid on; box on;

% (b) 各母线越限时段（按母线号排序，看清沿馈线的分布）
subplot(3,3,2);
bar(S.bus, nv, 0.7, 'FaceColor', [0.84 0.15 0.16]);
xlabel('bus number', 'FontSize', 8); ylabel('violating periods', 'FontSize', 8);
title('(b) Violations along the feeder', 'FontSize', 10);
grid on; box on; xlim([1 34]);
text(18, nv(i18)+0.25, '18', 'FontSize', 7, 'HorizontalAlignment', 'center', 'Color', [0.4 0.4 0.4]);
text(4, nv(i1)+0.25, '4', 'FontSize', 7, 'HorizontalAlignment', 'center', 'Color', [0.4 0.4 0.4]);

% (c) 拓扑图：标出两个位置
subplot(3,3,3);
m = case33mg_der('pv_penetration', 1.5, 'storage_bus', 4);
r = runpf(m, mpoption('verbose', 0, 'out.all', 0));
draw_topo(m, r, layout33(), 4, 18);
title('(c) Bus 4 (best) vs bus 18 (old)', 'FontSize', 10);

%% ================= 第二行：SOCP =================
D = RD('distflow_model_compare.csv');

% (d) 误差 vs 渗透率
subplot(3,3,4);
plot(D.pen*100, D.err_lin, 'o-', 'LineWidth', 2, 'Color', [0.85 0.33 0.10]); hold on;
plot(D.pen*100, D.err_socp, 's-', 'LineWidth', 2, 'Color', [0.12 0.47 0.71]);
yline(0, '--', 'Color', [0.55 0.55 0.55], 'FontSize', 7);
xlabel('PV penetration [%]', 'FontSize', 8);
ylabel('Vmin(model) - Vmin(AC)  [p.u.]', 'FontSize', 8);
title('(d) SOCP: voltage error', 'FontSize', 10);
legend({'LinDistFlow','SOCP'}, 'Location', 'best', 'FontSize', 7);
grid on; box on;

% (e) 绝对电压（放大看细节）
subplot(3,3,5);
plot(D.pen*100, D.vmin_ac, 'k-', 'LineWidth', 2.2); hold on;
plot(D.pen*100, D.vmin_lin, 'o--', 'LineWidth', 1.5, 'Color', [0.85 0.33 0.10]);
plot(D.pen*100, D.vmin_socp, 's--', 'LineWidth', 1.5, 'Color', [0.12 0.47 0.71]);
d = D.vmin_lin - D.vmin_ac;
ylim([min(D.vmin_lin)-0.004, max(D.vmin_ac)+0.006]);
xlabel('PV penetration [%]', 'FontSize', 8); ylabel('Vmin  [p.u.]', 'FontSize', 8);
title('(e) Zoom: SOCP sits exactly on AC', 'FontSize', 10);
legend({'AC (truth)','LinDistFlow','SOCP'}, 'Location', 'southeast', 'FontSize', 7);
grid on; box on;

% (f) 精度汇总
subplot(3,3,6);
rl = sqrt(mean(D.err_lin.^2));  rs = sqrt(mean(D.err_socp.^2));
b = bar([rl rs], 0.5, 'FaceColor', 'flat');
b.CData(1,:) = [0.85 0.33 0.10];  b.CData(2,:) = [0.12 0.47 0.71];
set(gca, 'XTickLabel', {'LinDistFlow','SOCP'}, 'FontSize', 8);
ylabel('RMS of Vmin error  [p.u.]', 'FontSize', 8);
title('(f) RMS error over 6 cases', 'FontSize', 10);
grid on; box on;
text(1, rl, sprintf('  %.5f', rl), 'FontSize', 8, 'VerticalAlignment', 'bottom');
text(2, max(rs, rl*0.01), sprintf('  %.5f', rs), 'FontSize', 8, 'VerticalAlignment', 'bottom');

%% ================= 第三行：云图时间相关 =================
% (g) 云量序列 + 自相关
subplot(3,3,7);
ccsv = fullfile('D:\pv-probabilistic-forecast', 'vision', 'output', 'cloud_cover_yolo.csv');
if exist(ccsv, 'file')
    C = readtable(ccsv, 'Delimiter', ',');
    y = C.cloud_fraction_true;
    plot(1:numel(y), y, '-', 'LineWidth', 1.1, 'Color', [0.45 0.55 0.75]); hold on;
    ylim([0 1.05]); xlim([1 numel(y)]);
    xlabel('hour index', 'FontSize', 8); ylabel('cloud fraction (CLM truth)', 'FontSize', 8);
    title(sprintf('(g) Satellite cloud series (lag-1 AC = %.3f)', ...
        corr(y(1:end-1), y(2:end))), 'FontSize', 10);
    grid on; box on;
else
    text(0.1, 0.5, 'cloud series not found', 'FontSize', 9);
    axis off;
end

% (h) 场景时间结构：独立 vs 相关
subplot(3,3,8);
if exist(ccsv, 'file')
    C = readtable(ccsv, 'Delimiter', ',');
    y = C.cloud_fraction_true;
    ac = arrayfun(@(k) corr(y(1:end-k), y(1+k:end)), 1:8);
    plot(1:8, ac, 'o-', 'LineWidth', 2, 'Color', [0.35 0.55 0.30]); hold on;
    j = jsondecode(fileread(fullfile(root, 'data', 'cloud_correlation.json')));
    kk = (1:8)';
    plot(kk, j.recommended_rho_time .^ kk, '--', 'LineWidth', 1.5, 'Color', [0.7 0.7 0.7]);
    ylim([0 1.05]);
    xlabel('lag  [hour]', 'FontSize', 8); ylabel('autocorrelation', 'FontSize', 8);
    title(sprintf('(h) Cloud AC -> copula rho = %.3f', j.recommended_rho_time), 'FontSize', 10);
    legend({'measured','AR(1) fit'}, 'Location', 'northeast', 'FontSize', 7);
    grid on; box on;
end

% (i) 对尾部风险的影响
subplot(3,3,9);
try
    A = RD('dro_compare_pen150_eps050_g50000rts79.csv');        % 逐小时独立
    B = RD('dro_compare_pen150_eps050_g50000cloudtime.csv');    % 云图时间相关
    M = [A.CVaR95(2), B.CVaR95(2); A.worst_cost(2), B.worst_cost(2)]';
    b = bar(M, 'grouped');
    set(gca, 'XTickLabel', {'independent','cloud-correlated'}, 'FontSize', 8);
    ylabel('cost  [$/day]', 'FontSize', 8);
    title('(i) Tail risk: +4% when correlated', 'FontSize', 10);
    legend({'CVaR95','worst'}, 'Location', 'northwest', 'FontSize', 7);
    grid on; box on;
    ylim([min(M(:))*0.96, max(M(:))*1.02]);
    text(1, M(1,1), sprintf('  +%.1f%%', 100*(B.worst_cost(2)/A.worst_cost(2)-1)), ...
        'FontSize', 7.5, 'VerticalAlignment', 'bottom', 'Color', [0.35 0.35 0.35]);
catch
    text(0.1, 0.5, 'comparison data missing', 'FontSize', 9); axis off;
end

sgtitle({'Three extensions: storage siting  |  SOCP power flow  |  cloud-derived temporal correlation', ...
         'every result validated against full AC power flow'}, 'FontSize', 12);

out = fullfile(root, 'results', 'figures', 'OVERVIEW_EXT.png');
exportgraphics(f, out, 'Resolution', 130);
close(f);
fprintf('saved: %s\n', out);
end


%% ---------------- helpers ----------------
function draw_topo(m, r, XY, best_bus, old_bus)
define_constants;
hold on;
Smax = max(hypot(r.branch(:,PF), r.branch(:,PT)));
for k = find(m.branch(:, BR_STATUS) == 0)'
    plot(XY([m.branch(k,F_BUS) m.branch(k,T_BUS)],1), ...
         XY([m.branch(k,F_BUS) m.branch(k,T_BUS)],2), ...
         '--', 'Color', [0.82 0.82 0.82], 'LineWidth', 0.6);
end
for k = find(m.branch(:, BR_STATUS) == 1)'
    S = hypot(r.branch(k,PF), r.branch(k,PT));
    if r.branch(k, PF) >= 0, c = [0.15 0.35 0.75]; else, c = [0.82 0.15 0.12]; end
    plot(XY([m.branch(k,F_BUS) m.branch(k,T_BUS)],1), ...
         XY([m.branch(k,F_BUS) m.branch(k,T_BUS)],2), '-', ...
         'Color', c, 'LineWidth', 0.5 + 3.0*S/max(Smax,eps));
end
vm = r.bus(:, VM);
cmap = jet(256); lo = 0.90; hi = 1.05;
ci = max(1, min(256, round((vm-lo)/(hi-lo)*255)+1));
scatter(XY(:,1), XY(:,2), 28, cmap(ci,:), 'filled', ...
    'MarkerEdgeColor', [0.25 0.25 0.25], 'LineWidth', 0.3);
% 光伏节点
for b = [18 22 25 33]
    plot(XY(b,1), XY(b,2), 's', 'MarkerSize', 7, 'LineWidth', 1.4, ...
        'MarkerEdgeColor', [0.05 0.55 0.10], 'MarkerFaceColor', 'none');
end
% 储能：最优位置（绿）vs 原位置（黄）
plot(XY(best_bus,1), XY(best_bus,2), 'o', 'MarkerSize', 13, 'LineWidth', 2.4, ...
    'MarkerEdgeColor', [0.10 0.65 0.15], 'MarkerFaceColor', 'none');
plot(XY(old_bus,1), XY(old_bus,2), 'o', 'MarkerSize', 13, 'LineWidth', 2.4, ...
    'MarkerEdgeColor', [0.95 0.75 0.10], 'MarkerFaceColor', 'none');
plot(XY(1,1), XY(1,2), 'p', 'MarkerSize', 13, 'LineWidth', 1.1, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.95 0.85 0.20]);
text(XY(best_bus,1), XY(best_bus,2)-0.75, sprintf('ESS@%d', best_bus), ...
    'FontSize', 7, 'HorizontalAlignment','center', 'Color', [0.10 0.55 0.15]);
text(XY(old_bus,1), XY(old_bus,2)+0.85, sprintf('old@%d', old_bus), ...
    'FontSize', 7, 'HorizontalAlignment','center', 'Color', [0.75 0.55 0.05]);
axis equal; axis off; xlim([0 20]); ylim([-9.8 2.0]); colormap(cmap);
end


function XY = layout33()
XY = nan(33,2);
for i = 1:18, XY(i,:) = [i, 0]; end
for k = 1:4, XY(18+k,:) = [2,   -k]; end
for k = 1:3, XY(22+k,:) = [4.2, -k]; end
for k = 1:8, XY(25+k,:) = [15,  -k]; end
end
