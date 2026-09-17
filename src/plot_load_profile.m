function plot_load_profile()
%PLOT_LOAD_PROFILE  负荷曲线对比：原合成曲线 vs IEEE RTS-79 标准曲线
%
%   跑法：matlab -batch "cd('D:\microgrid-dro\src'); plot_load_profile"
%   输出：results/figures/load_profile_compare.png
%
%   重点看"峰均比"（峰值 ÷ 均值）—— 它决定了峰值会不会顶到设备上限。

root = fileparts(fileparts(mfilename('fullpath')));
outdir = fullfile(root, 'results', 'figures');
if ~exist(outdir, 'dir'), mkdir(outdir); end

% 原合成曲线（build_most_data.default_load_shape 里那条，此处照抄以便对比）
synth = [0.62 0.58 0.55 0.54 0.55 0.62 0.72 0.82 0.89 0.92 0.92 0.90 ...
         0.87 0.86 0.88 0.92 1.00 1.12 1.20 1.22 1.13 1.00 0.84 0.71]';
synth = synth / mean(synth);

V = readtable(fullfile(root, 'data', 'load_profile.variants.csv'), 'Delimiter', ',');
cur = readtable(fullfile(root, 'data', 'load_profile.csv'), 'Delimiter', ',');
hr = (1:24)';

f = figure('Position', [80 80 1420 460], 'Color', 'w', 'Visible', 'off');

%% (a) 新旧对比 —— 重点是峰均比
subplot(1,3,1);
plot(hr, synth, '-o', 'LineWidth', 2.2, 'MarkerSize', 5, ...
    'Color', [0.85 0.33 0.10], 'MarkerFaceColor', [0.85 0.33 0.10]); hold on;
plot(hr, cur.load_pu, '-s', 'LineWidth', 2.2, 'MarkerSize', 5, ...
    'Color', [0.12 0.47 0.71], 'MarkerFaceColor', [0.12 0.47 0.71]);
yline(1, '--', 'daily mean = 1.0', 'Color', [0.5 0.5 0.5], 'FontSize', 8);
xlabel('hour'); ylabel('load shape factor  [p.u.]');
title('(a) Old synthetic  vs  IEEE RTS-79', 'FontSize', 11);
legend({sprintf('old synthetic  (peak/mean = %.2f)', max(synth)/mean(synth)), ...
        sprintf('RTS-79 winter Wd  (peak/mean = %.2f)', max(cur.load_pu)/mean(cur.load_pu))}, ...
    'Location', 'northwest', 'FontSize', 8);
grid on; box on; xlim([1 24]); ylim([0.5 1.45]);

%% (b) RTS-79 的 6 种变体
subplot(1,3,2);
names = {'winter_wd','winter_we','summer_wd','summer_we','springfall_wd','springfall_we'};
lbl   = {'Winter Wd','Winter We','Summer Wd','Summer We','Spr/Fall Wd','Spr/Fall We'};
cols  = lines(numel(names));
hold on;
for k = 1:numel(names)
    plot(hr, V.(names{k}), '-', 'LineWidth', 1.6, 'Color', cols(k,:), 'DisplayName', lbl{k});
end
yline(1, '--', 'HandleVisibility', 'off', 'Color', [0.6 0.6 0.6], 'FontSize', 7);
xlabel('hour'); ylabel('load shape factor  [p.u.]');
title('(b) IEEE RTS-79 all six variants', 'FontSize', 11);
legend('Location', 'northwest', 'FontSize', 7, 'NumColumns', 2);
grid on; box on; xlim([1 24]);

%% (c) 峰均比 —— 决定峰值会不会顶到设备上限
subplot(1,3,3);
allv = [max(synth)/mean(synth); cellfun(@(n) max(V.(n))/mean(V.(n)), names)'];
alln = [{'old synthetic'}, lbl];
b = barh(allv, 0.65, 'FaceColor', 'flat');
for k = 1:numel(allv)
    if k == 1, b.CData(k,:) = [0.85 0.33 0.10]; else, b.CData(k,:) = cols(k-1,:); end
    text(allv(k)+0.012, k, sprintf('%.2f', allv(k)), 'FontSize', 8, ...
        'VerticalAlignment', 'middle');
end
set(gca, 'YTick', 1:numel(allv), 'YTickLabel', alln, 'YDir', 'reverse', 'FontSize', 8);
xlabel('peak / mean ratio');
title('(c) Peak-to-mean ratio', 'FontSize', 11);
xlim([0 max(allv)*1.18]);
grid on; box on;
text(0.03, numel(allv)+0.4, '峰均比越高，峰值越容易顶到设备上限（变电站/线路）', ...
    'FontSize', 7, 'Color', [0.35 0.35 0.35]);

sgtitle('Load profile: synthetic (old) \rightarrow IEEE RTS-79 standard (new)', 'FontSize', 12);

out = fullfile(outdir, 'load_profile_compare.png');
exportgraphics(f, out, 'Resolution', 140);
close(f);
fprintf('saved: %s\n', out);

fprintf('\n--- 峰均比 ---\n');
fprintf('  %-16s %.3f\n', 'old synthetic', max(synth)/mean(synth));
for k = 1:numel(names)
    fprintf('  %-16s %.3f\n', lbl{k}, max(V.(names{k}))/mean(V.(names{k})));
end
end
