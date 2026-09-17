function T = sweep_pv_penetration(varargin)
%SWEEP_PV_PENETRATION  光伏渗透率扫描：网损 / 电压 / 支路负载
%
%   T = sweep_pv_penetration()
%   T = sweep_pv_penetration('vmax', 1.10)
%   T = sweep_pv_penetration('pen', 0:0.1:2.0, 'plot', false)
%
%   选项：
%     'pen'   渗透率向量（装机/总负荷）   (默认 0:0.05:2.5)
%     'vmax'  电压上限 p.u.               (默认 1.05)
%     'vmin'  电压下限 p.u.               (默认 0.90)
%     'save'  是否把数据存成 CSV          (默认 true)
%     'plot'  是否出图                    (默认 true)
%     'thermal' 是否加支路热限值          (默认 true)
%
%   输出：
%     T  table，一行一个渗透率；同时（save=true 时）写入
%        results/data/sweep_pv_vmax{上限}.csv  —— 文件名带 vmax 标记，
%        这样「1.05 vs 1.10」的对比可以直接读两份 CSV 叠一起画。
%
%   ★ 电压限值的选取
%     上界用 1.05（ANSI C84.1 Range A 的配电网实际限值），而不是 case33mg 默认的
%     1.1（那是输电系统的默认值）。实测本馈线过 1.05 在 172% 渗透率，过 1.10 要 245%
%     —— 差 40% 以上，直接决定「光伏承载力」这个结论。
%     下界保持 0.9：本馈线基态就有 21/33 个节点低于 0.95（Vmin=0.9038@bus18），
%     用 0.95 会让基态本身不可行。低电压是这条馈线的固有问题，而光伏恰好能改善它。

%% ---------- 选项 ----------
opt = struct('pen', 0:0.05:2.5, 'vmax', 1.05, 'vmin', 0.90, ...
             'save', true, 'plot', true, 'thermal', true);
for k = 1:2:numel(varargin)
    name = varargin{k};
    if ~ischar(name) || ~isfield(opt, name)
        error('sweep_pv_penetration: 未知选项 ''%s''', num2str(name));
    end
    opt.(name) = varargin{k + 1};
end

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'cases'));
define_constants;
mpopt = mpoption('verbose', 0, 'out.all', 0);

pens = opt.pen(:);
n = numel(pens);
total_load = sum(case33mg().bus(:, PD));
base_gen_rows = size(case33mg().gen, 1);

pen          = nan(n,1);
pv_total_MW  = nan(n,1);
pv_actual_MW = nan(n,1);
loss_kW      = nan(n,1);
loss_pct     = nan(n,1);
slack_kW     = nan(n,1);
vmax_pu      = nan(n,1);
vmax_bus     = nan(n,1);
vmin_pu      = nan(n,1);
vmin_bus     = nan(n,1);
maxS_MVA     = nan(n,1);
max_load_pct = nan(n,1);
n_over_vmax  = nan(n,1);
n_under_vmin = nan(n,1);
n_congested  = nan(n,1);
converged    = false(n,1);

for k = 1:n
    m = case33mg_der('pv_penetration', pens(k), 'thermal', opt.thermal, ...
                     'vmax', opt.vmax, 'vmin', opt.vmin);
    r = runpf(m, mpopt);
    converged(k) = logical(r.success);
    if ~r.success, continue; end

    S = hypot(r.branch(:, PF), r.branch(:, PT));
    rate = m.branch(1, RATE_A);
    if rate <= 0, rate = Inf; end

    pen(k)          = pens(k);
    pv_total_MW(k)  = pens(k) * total_load;
    ipv             = base_gen_rows + (1:4);
    if size(r.gen, 1) >= max(ipv)
        pv_actual_MW(k) = sum(r.gen(ipv, PG));
    else
        pv_actual_MW(k) = 0;
    end
    loss_kW(k)      = sum(r.branch(:, PF) + r.branch(:, PT)) * 1000;
    loss_pct(k)     = loss_kW(k) / (total_load * 1000) * 100;
    slack_kW(k)     = r.gen(1, PG) * 1000;
    [vmax_pu(k), vmax_bus(k)] = max(r.bus(:, VM));
    [vmin_pu(k), vmin_bus(k)] = min(r.bus(:, VM));
    maxS_MVA(k)     = max(S);
    max_load_pct(k) = maxS_MVA(k) / rate * 100;
    n_over_vmax(k)  = sum(r.bus(:, VM) > opt.vmax + 1e-9);
    n_under_vmin(k) = sum(r.bus(:, VM) < opt.vmin - 1e-9);
    n_congested(k)  = sum(S > rate * (1 + 1e-9));
end

T = table(pen, pv_total_MW, pv_actual_MW, loss_kW, loss_pct, slack_kW, ...
          vmax_pu, vmax_bus, vmin_pu, vmin_bus, maxS_MVA, max_load_pct, ...
          n_over_vmax, n_under_vmin, n_congested, converged);

%% ---------- 存数据 ----------
tag = sprintf('vmax%02.0f', opt.vmax * 100);
if opt.save
    outdir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results', 'data');
    if ~exist(outdir, 'dir'), mkdir(outdir); end
    out = fullfile(outdir, sprintf('sweep_pv_%s.csv', tag));
    writetable(T, out);
    fprintf('data saved: %s\n', out);
    % 元数据另存，便于对比时确认两边的设置一致
    meta = struct('vmax', opt.vmax, 'vmin', opt.vmin, 'thermal', opt.thermal, ...
                  'total_load_MW', total_load, 'pen_min', min(pens), ...
                  'pen_max', max(pens), 'n_points', n, ...
                  'generated', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fid = fopen(fullfile(outdir, sprintf('sweep_pv_%s.meta.json', tag)), 'w');
    fprintf(fid, '%s', jsonencode(meta));
    fclose(fid);
end

%% ---------- 出图 ----------
kk = converged;
% 临界渗透率先算出来 —— 摘要段要用，不能只在出图分支里算
kmin  = find(loss_kW == min(loss_kW(kk)), 1);
pen_v = crit_pen(pen(kk), vmax_pu(kk), opt.vmax, +1);
pen_c = crit_pen(pen(kk), max_load_pct(kk), 100, +1);

if opt.plot && any(kk)
    f = figure('Position', [80 80 1450 430], 'Color', 'w', 'Visible', 'off');

    subplot(1,3,1);
    plot(pen(kk)*100, loss_kW(kk), '-', 'LineWidth', 2, 'Color', [0.12 0.47 0.71]); hold on;
    plot(pen(kmin)*100, loss_kW(kmin), 'o', 'MarkerSize', 9, ...
        'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerEdgeColor', 'none');
    text(pen(kmin)*100 + 6, loss_kW(kmin) + 18, ...
        sprintf('min %.0f kW @ %.0f%%', loss_kW(kmin), pen(kmin)*100), 'FontSize', 9);
    yline(loss_kW(1), '--', sprintf('base %.0f kW', loss_kW(1)), ...
        'Color', [0.5 0.5 0.5], 'FontSize', 8);
    xlabel('PV penetration  [% of total load]');
    ylabel('Network loss  [kW]');
    title('(a) Loss: drops, then rises');
    grid on; box on; xlim([0 max(pen)*100]);

    subplot(1,3,2);
    plot(pen(kk)*100, vmax_pu(kk), '-', 'LineWidth', 2, 'Color', [0.84 0.15 0.16]); hold on;
    plot(pen(kk)*100, vmin_pu(kk), '-', 'LineWidth', 2, 'Color', [0.12 0.47 0.71]);
    yline(opt.vmax, '--', sprintf('Vmax %.2f', opt.vmax), 'Color', [0.84 0.15 0.16], ...
        'FontSize', 8, 'LabelHorizontalAlignment', 'left');
    yline(opt.vmin, '--', sprintf('Vmin %.2f', opt.vmin), 'Color', [0.12 0.47 0.71], ...
        'FontSize', 8, 'LabelHorizontalAlignment', 'left');
    if ~isnan(pen_v)
        xline(pen_v*100, '-.', sprintf('%.0f%%', pen_v*100), 'Color', [0.84 0.15 0.16], ...
            'FontSize', 9, 'LabelOrientation', 'horizontal');
    end
    xlabel('PV penetration  [% of total load]');
    ylabel('Bus voltage  [p.u.]');
    title(sprintf('(b) Vmax %.2f binds at %.0f%%', opt.vmax, pen_v*100));
    grid on; box on; xlim([0 max(pen)*100]);
    legend({'Vmax', 'Vmin'}, 'Location', 'southeast', 'FontSize', 9);

    subplot(1,3,3);
    plot(pen(kk)*100, max_load_pct(kk), '-', 'LineWidth', 2, 'Color', [0.30 0.55 0.20]); hold on;
    yline(100, '--r', 'thermal rating', 'FontSize', 8, 'LabelHorizontalAlignment', 'left');
    if ~isnan(pen_c)
        xline(pen_c*100, '-.', sprintf('%.0f%%', pen_c*100), 'Color', [0.30 0.55 0.20], ...
            'FontSize', 9, 'LabelOrientation', 'horizontal');
    end
    xlabel('PV penetration  [% of total load]');
    ylabel('Max branch loading  [%]');
    title('(c) Thermal: not binding at realistic PV');
    grid on; box on; xlim([0 max(pen)*100]);
    ylim([0, max(max_load_pct(kk)) * 1.3]);

    sgtitle(sprintf('case33mg\\_der : PV penetration sweep  (Vmax=%.2f, Vmin=%.2f)', ...
        opt.vmax, opt.vmin), 'FontSize', 13);

    outdir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results', 'figures');
    if ~exist(outdir, 'dir'), mkdir(outdir); end
    outfig = fullfile(outdir, sprintf('sweep_pv_%s.png', tag));
    exportgraphics(f, outfig, 'Resolution', 140);
    close(f);
    fprintf('figure saved: %s\n', outfig);
end

%% ---------- 摘要 ----------
fprintf('\n--- summary (Vmax=%.2f, Vmin=%.2f) ---\n', opt.vmax, opt.vmin);
fprintf('base loss (0%% PV)        = %.1f kW\n', loss_kW(1));
fprintf('min loss                  = %.1f kW @ %.0f%%\n', min(loss_kW(kk)), ...
    pen(find(loss_kW == min(loss_kW(kk)), 1))*100);
fprintf('slack sign reversal       = %s\n', pstr(crit_pen(pen(kk), slack_kW(kk), 0, -1)));
fprintf('Vmax %.2f first exceeded    = %s\n', opt.vmax, pstr(pen_v));
% Vmin 的 0.90 下限基态就没越（基态 0.9038），所以报 ANSI 的 0.95 更有意义：
% 它是"光伏把馈线电压治好了"的标志点。
fprintf('Vmin >= 0.95 (ANSI) at      = %s\n', ...
    pstr(crit_pen(pen(kk), vmin_pu(kk), 0.95, +1)));
fprintf('thermal congestion        = %s\n', pstr(pen_c));
fprintf('max loading over sweep    = %.1f%%\n', max(max_load_pct(kk)));
fprintf('\n');
end


function p = crit_pen(pens, y, thr, dir)
if nargin < 4, dir = 1; end
p = NaN;
s = y - thr;
for k = 2:numel(s)
    if ~isfinite(s(k-1)) || ~isfinite(s(k)), continue; end
    rising  = (s(k-1) < 0) && (s(k) >= 0);
    falling = (s(k-1) >= 0) && (s(k) < 0);
    if (dir > 0 && rising) || (dir < 0 && falling)
        t = s(k-1) / (s(k-1) - s(k));
        p = pens(k-1) + t * (pens(k) - pens(k-1));
        return
    end
end
end


function s = pstr(p)
if isnan(p)
    s = 'not reached within sweep';
else
    s = sprintf('%.1f%% penetration', p * 100);
end
end
