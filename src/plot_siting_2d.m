function plot_siting_2d()
%PLOT_SITING_2D  储能选址 + 定容二维扫描的可视化
%
%   跑法：matlab -batch "cd('D:\microgrid-dro\src'); plot_siting_2d"
%   输出：results/figures/storage_siting_2d.png
%
%   读 results/data/storage_siting_pen150_g042.csv（由 optimize_storage_siting 生成）

root = fileparts(fileparts(mfilename('fullpath')));
T = readtable(fullfile(root, 'results', 'data', 'storage_siting_pen150_g042.csv'), ...
              'Delimiter', ',');
T.nviol = T.n_over + T.n_under;
caps = unique(T.cap_MW);
buses = unique(T.bus);
nc = numel(caps);  nb = numel(buses);

% 拼成 母线 x 容量 的矩阵
V = nan(nb, nc);   C = nan(nb, nc);   L = nan(nb, nc);
for i = 1:height(T)
    r = find(buses == T.bus(i));  c = find(caps == T.cap_MW(i));
    V(r,c) = T.nviol(i);  C(r,c) = T.cost_usd_d(i);  L(r,c) = T.max_load_pct(i);
end

f = figure('Position', [60 60 1520 470], 'Color', 'w', 'Visible', 'off');

% (a) 越限时段数热力图
subplot(1,3,1);
imagesc(V');
colormap(gca, flipud(hot));
cb = colorbar; cb.Label.String = 'violating periods'; cb.FontSize = 8;
set(gca, 'XTick', 1:4:nb, 'XTickLabel', buses(1:4:end), ...
         'YTick', 1:nc, 'YTickLabel', arrayfun(@(c) sprintf('%.1f MW', c), caps, 'UniformOutput', false), ...
         'FontSize', 8);
xlabel('storage bus', 'FontSize', 9); ylabel('storage capacity', 'FontSize', 9);
title('(a) Voltage violations (darker = worse)', 'FontSize', 11);
% 标出**综合最优**（按成本 + 1e4×越限 + 500×切负荷 的分数，与排名的口径一致）。
% 注意不能只按"越限最少"标 —— 那样会选中 bus 6/0.5MW（同样 0 越限但成本更高），
% 与排名表不一致。
T.score = T.cost_usd_d + 1e4 * T.nviol + 500 * T.shed_MWh;
[~, ibestT] = min(T.score);
rb = find(buses == T.bus(ibestT));  cb2 = find(caps == T.cap_MW(ibestT));
subplot(1,3,1);
hold on;
rectangle('Position', [rb-0.5, cb2-0.5, 1, 1], 'EdgeColor', [0.1 0.6 0.1], 'LineWidth', 2.5);
text(rb, cb2-0.75, sprintf('best: bus %d / %.1f MW', buses(rb), caps(cb2)), ...
    'FontSize', 7.5, 'HorizontalAlignment', 'center', 'Color', [0.1 0.5 0.1]);

% (b) 成本热力图
subplot(1,3,2);
imagesc(C');
colormap(gca, flipud(summer));
cb = colorbar; cb.Label.String = 'cost  [$/day]'; cb.FontSize = 8;
set(gca, 'XTick', 1:4:nb, 'XTickLabel', buses(1:4:end), ...
         'YTick', 1:nc, 'YTickLabel', arrayfun(@(c) sprintf('%.1f MW', c), caps, 'UniformOutput', false), ...
         'FontSize', 8);
xlabel('storage bus', 'FontSize', 9); ylabel('storage capacity', 'FontSize', 9);
title('(b) Cost (uniform in bus — only 0.4% spread)', 'FontSize', 11);

% (c) 三条代表母线的容量-越限曲线
subplot(1,3,3);
pick = [4 7 18];
cols = [0.10 0.60 0.15; 0.20 0.45 0.75; 0.85 0.55 0.10];
hold on;
for k = 1:numel(pick)
    r = find(buses == pick(k));
    if isempty(r), continue; end
    plot(caps, V(r,:), 'o-', 'LineWidth', 2.2, 'MarkerSize', 8, ...
        'Color', cols(k,:), 'MarkerFaceColor', cols(k,:));
end
xlabel('storage capacity  [MW]', 'FontSize', 9);
ylabel('violating periods', 'FontSize', 9);
title('(c) Capacity helps only at the right bus', 'FontSize', 11);
legend(arrayfun(@(b) sprintf('bus %d', b), pick, 'UniformOutput', false), ...
    'Location', 'northeast', 'FontSize', 8);
grid on; box on; xlim([0.3 2.2]); ylim([-0.5 10]);

sgtitle({'Storage joint siting + sizing: 32 buses \times 3 capacities = 96 combos', ...
    'Siting decides WHETHER the problem is solved; sizing decides IF IT IS ENOUGH'}, ...
    'FontSize', 12);

out = fullfile(root, 'results', 'figures', 'storage_siting_2d.png');
exportgraphics(f, out, 'Resolution', 140);
close(f);
fprintf('saved: %s\n', out);
end
