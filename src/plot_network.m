function plot_network(varargin)
%PLOT_NETWORK  IEEE 33 节点馈线拓扑图：DER 位置 + 潮流方向
%
%   plot_network()
%   plot_network('pv_penetration', 1.5)
%
%   画两张图并排：
%     (a) 基态（无光伏）：潮流从变电站单向流出，末端电压被拉低
%     (b) 高光伏：末端光伏反送，潮流方向翻转
%
%   母线颜色 = 电压（jet 色标）；支路粗细 = 有功潮流大小；
%   蓝色 = 正向（变电站 -> 馈线），红色 = 反向（馈线 -> 变电站）。
%   绿色方块 = 光伏，青色圆点 = 储能，红色虚线 = 常开联络开关。

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'cases'));
define_constants;
mpopt = mpoption('verbose', 0, 'out.all', 0);

opt = struct('pv_penetration', 1.5, 'save', true);
for k = 1:2:numel(varargin)
    if ~isfield(opt, varargin{k}), error('plot_network: 未知选项'); end
    opt.(varargin{k}) = varargin{k+1};
end

f = figure('Position', [60 40 1500 700], 'Color', 'w', 'Visible', 'off');

% 先跑两次潮流，取**全局**电压范围 —— 两图共用同一色标才能横向比较
M = cell(1,2); R = cell(1,2);
M{1} = case33mg_der('pv_penetration', 0);            R{1} = runpf(M{1}, mpopt);
M{2} = case33mg_der('pv_penetration', opt.pv_penetration); R{2} = runpf(M{2}, mpopt);
vlo = min(0.90, min([R{1}.bus(:,VM); R{2}.bus(:,VM)]));
vhi = max(1.05, max([R{1}.bus(:,VM); R{2}.bus(:,VM)]));

subplot(1,2,1);
draw_one(M{1}, R{1}, vlo, vhi, '(a) base case : no PV (one-way flow)');
subplot(1,2,2);
h = draw_one(M{2}, R{2}, vlo, vhi, ...
    sprintf('(b) PV = %.0f%% : reverse flow at the ends', opt.pv_penetration*100));

% 单一图例，放在整张图下方（两面板共用同一套符号与色标）
lg = legend(h, {'substation (slack)', 'PV (unity PF)', 'storage', ...
    'tie switch (open)', 'flow: forward', 'flow: REVERSE'}, ...
    'Orientation', 'horizontal', 'FontSize', 9);
lg.Units = 'normalized';
lg.Position = [0.16, 0.02, 0.55, 0.05];

cb = colorbar('Position', [0.93, 0.30, 0.012, 0.42]);
cb.Label.String = 'Bus voltage [p.u.]  (shared scale)';
clim([vlo vhi]);

sgtitle('IEEE 33-bus feeder : topology, DER siting and power flow direction', 'FontSize', 13);

outdir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results', 'figures');
if ~exist(outdir, 'dir'), mkdir(outdir); end
out = fullfile(outdir, sprintf('network_pv%03.0f.png', opt.pv_penetration*100));
exportgraphics(f, out, 'Resolution', 140);
close(f);
fprintf('figure saved: %s\n', out);
end


function h = draw_one(m, r, vlo, vhi, ttl)
define_constants;
XY = layout33();
nb = size(m.bus, 1);
vm = r.bus(:, VM);

hold on;
Smax = max(hypot(r.branch(:,PF), r.branch(:,PT)));
% 常开联络开关：灰色虚线（不承载潮流）
for k = find(m.branch(:, BR_STATUS) == 0)'
    plot(XY([m.branch(k,F_BUS) m.branch(k,T_BUS)],1), ...
         XY([m.branch(k,F_BUS) m.branch(k,T_BUS)],2), ...
         '--', 'Color', [0.78 0.78 0.78], 'LineWidth', 0.8);
end
% 闭合支路：粗细 = 潮流大小，颜色 = 方向
for k = find(m.branch(:, BR_STATUS) == 1)'
    S = hypot(r.branch(k,PF), r.branch(k,PT));
    lw = 0.5 + 4.0 * S / max(Smax, eps);
    if r.branch(k, PF) >= 0
        c = [0.15 0.35 0.75];
    else
        c = [0.82 0.15 0.12];
    end
    plot(XY([m.branch(k,F_BUS) m.branch(k,T_BUS)],1), ...
         XY([m.branch(k,F_BUS) m.branch(k,T_BUS)],2), '-', 'Color', c, 'LineWidth', lw);
end

% 母线（共用色标）
cmap = jet(256);
ci = max(1, min(256, round((vm - vlo) / (vhi - vlo) * 255) + 1));
scatter(XY(:,1), XY(:,2), 95, cmap(ci,:), 'filled', ...
    'MarkerEdgeColor', [0.2 0.2 0.2], 'LineWidth', 0.6);

% 母线编号
for b = 1:nb
    if b <= 18 || mod(b,2)==0
        text(XY(b,1), XY(b,2)+0.30, num2str(b), 'FontSize', 7, ...
            'HorizontalAlignment', 'center', 'Color', [0.25 0.25 0.25]);
    end
end

% DER 标记（绿方 = 光伏，青圈 = 储能）
for b = [18 22 25 33]
    plot(XY(b,1), XY(b,2), 's', 'MarkerSize', 13, 'LineWidth', 2.2, ...
        'MarkerEdgeColor', [0.05 0.55 0.10], 'MarkerFaceColor', 'none');
end
plot(XY(18,1), XY(18,2), 'o', 'MarkerSize', 9, 'LineWidth', 2.0, ...
    'MarkerEdgeColor', [0.05 0.55 0.75], 'MarkerFaceColor', 'none');
% 变电站
plot(XY(1,1), XY(1,2), 'p', 'MarkerSize', 22, 'LineWidth', 1.5, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.95 0.85 0.20]);

axis equal; axis off;
xlim([0 22]); ylim([-9.5 2.2]);
title(ttl, 'FontSize', 11);
colormap(cmap);

% 返回句柄（供整图图例使用）
h = [plot(nan, nan, 'p', 'MarkerSize', 14, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', [0.95 0.85 0.20]), ...
     plot(nan, nan, 's', 'MarkerSize', 11, 'LineWidth', 2.2, ...
        'MarkerEdgeColor', [0.05 0.55 0.10], 'MarkerFaceColor', 'none'), ...
     plot(nan, nan, 'o', 'MarkerSize', 8, 'LineWidth', 2.0, ...
        'MarkerEdgeColor', [0.05 0.55 0.75], 'MarkerFaceColor', 'none'), ...
     plot(nan, nan, '--', 'Color', [0.78 0.78 0.78], 'LineWidth', 1.2), ...
     plot(nan, nan, '-', 'Color', [0.15 0.35 0.75], 'LineWidth', 2.5), ...
     plot(nan, nan, '-', 'Color', [0.82 0.15 0.12], 'LineWidth', 2.5)];
end


function XY = layout33()
% IEEE 33 节点的手工布局：主干 1-18 水平，三条支线向下垂
XY = nan(33, 2);
for i = 1:18, XY(i,:) = [i, 0]; end              % 主干 1..18
for k = 1:4, XY(18+k,:) = [2,   -k]; end         % 支线 A: 19..22
for k = 1:3, XY(22+k,:) = [4.2, -k]; end         % 支线 B: 23..25
for k = 1:8, XY(25+k,:) = [15,  -k]; end         % 支线 C: 26..33
end
