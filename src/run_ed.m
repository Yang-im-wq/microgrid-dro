function T = run_ed(varargin)
%RUN_ED  确定性经济调度（单时段 AC OPF）
%
%   T = run_ed()
%   T = run_ed('pv_penetration', 2.0, 'avail', 0:0.1:1)
%   T = run_ed('reconfig', true)
%
%   选项：
%     'pv_penetration' 光伏装机渗透率          (默认 2.0)
%     'avail'          光伏可用系数向量 0~1    (默认 0:0.05:1)
%     'vmax'/'vmin'    电压限值                (默认 1.05 / 0.90)
%     'storage'        是否含储能              (默认 false，见下)
%     'reconfig'       是否对比网络配置        (默认 false)
%     'save'/'plot'                            (默认 true / true)
%
%   ★ 为什么单时段 OPF 必须关掉储能
%     单时段 OPF 没有 SOC 约束（SOC 是跨时段的耦合），储能又通常是零成本，
%     于是它会**任意充放**。实测 pen=100% 时储能自由放电 0.697 MW，凭空把
%     光伏挤出 0.697 MW，制造出根本不存在的"弃光"。
%     储能的真正建模在多时段调度 —— 见 P4 的 MOST。
%     所以本函数默认 'storage', false。
%
%   ★ 关于"弃光"的判定
%     本算例的平衡机（bus1，变电站）PMIN = 0，即**不允许向主网反送**。
%     所以光伏出力被硬性限制在「负荷 + 网损」以内，超出的必须弃掉。
%     这是一个真实的并网约束（很多分布式电源接入协议禁止或限制反送），
%     但它会主导结果 —— 想放开就把平衡机 PMIN 改负（见 'allow_export'）。
%
%   See also case33mg_der, runopf, sweep_pv_penetration, reconfigure.

%% ---------- 选项 ----------
opt = struct('pv_penetration', 2.0, 'avail', 0:0.05:1, ...
             'vmax', 1.05, 'vmin', 0.90, 'storage', false, ...
             'reconfig', false, 'allow_export', false, 'save', true, 'plot', true);
for k = 1:2:numel(varargin)
    name = varargin{k};
    if ~ischar(name) || ~isfield(opt, name)
        error('run_ed: 未知选项 ''%s''', num2str(name));
    end
    opt.(name) = varargin{k + 1};
end

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'cases'));
define_constants;
mpopt = mpoption('verbose', 0, 'out.all', 0);

avail = opt.avail(:);
n = numel(avail);
n_pv = 4;
base_rows = size(case33mg().gen, 1);
ipv = base_rows + (1:n_pv);

A          = nan(n,1);
pv_avail   = nan(n,1);
pv_disp    = nan(n,1);
pv_curt    = nan(n,1);
slack_MW   = nan(n,1);
loss_MW    = nan(n,1);
cost       = nan(n,1);
vmax_pu    = nan(n,1);
vmin_pu    = nan(n,1);
lmp_max    = nan(n,1);
lmp_min    = nan(n,1);
n_over     = nan(n,1);
ok         = false(n,1);

for k = 1:n
    m = case33mg_der('pv_penetration', opt.pv_penetration, 'storage', opt.storage, ...
                     'vmax', opt.vmax, 'vmin', opt.vmin);
    % 光伏可用出力上限 = 装机 × 可用系数（OPF 在 [0, 该上限] 内自由削减）
    m.gen(ipv, PMAX) = m.gen(ipv, PMAX) * avail(k);
    m.gen(ipv, PG)   = m.gen(ipv, PMAX);       % 潮流暖启动
    if opt.allow_export
        m.gen(1, PMIN) = -m.gen(1, PMAX);      % 允许反送
    end

    o = runopf(m, mpopt);
    ok(k) = logical(o.success);
    if ~o.success, continue; end

    A(k)        = avail(k);
    pv_avail(k) = sum(m.gen(ipv, PMAX));
    pv_disp(k)  = sum(o.gen(ipv, PG));
    pv_curt(k)  = max(pv_avail(k) - pv_disp(k), 0);
    slack_MW(k) = o.gen(1, PG);
    loss_MW(k)  = sum(o.branch(:, PF) + o.branch(:, PT));
    cost(k)     = o.f;
    vmax_pu(k)  = max(o.bus(:, VM));
    vmin_pu(k)  = min(o.bus(:, VM));
    lmp_max(k)  = max(o.bus(:, LAM_P));
    lmp_min(k)  = min(o.bus(:, LAM_P));
    n_over(k)   = sum(o.bus(:, VM) > opt.vmax + 1e-6);
end

T = table(A, pv_avail, pv_disp, pv_curt, slack_MW, loss_MW, cost, ...
          vmax_pu, vmin_pu, lmp_max, lmp_min, n_over, ok, ...
    'VariableNames', {'avail','pv_avail_MW','pv_disp_MW','pv_curt_MW', ...
                      'slack_MW','loss_MW','cost_usd_h','vmax_pu','vmin_pu', ...
                      'lmp_max','lmp_min','n_over_vmax','converged'});

%% ---------- 摘要 ----------
% 注意：弃光量在低可用率下有 ~1e-5 MW 量级的数值噪声（求解器容差），
% 用 1e-6 当阈值会误报"一开始就弃光"。取 1 kW 作为有意义的阈值。
TOL = 1e-3;
kk = ok;
k_curt_start = find(pv_curt(kk) > TOL, 1);
fprintf('\n=== deterministic ED (single period AC OPF) ===\n');
fprintf('PV penetration = %.0f%%  (installed %.3f MW, load %.3f MW)\n', ...
    opt.pv_penetration*100, opt.pv_penetration * sum(case33mg().bus(:,PD)), ...
    sum(case33mg().bus(:,PD)));
fprintf('storage = %d, allow_export = %d, Vmax = %.2f\n\n', ...
    opt.storage, opt.allow_export, opt.vmax);
if isempty(k_curt_start)
    fprintf('curtailment never occurs over the swept availability range\n');
else
    fprintf('curtailment starts at avail = %.2f (PV available %.3f MW)\n', ...
        A(k_curt_start), pv_avail(k_curt_start));
    fprintf('PV dispatch saturates at %.3f MW\n', pv_disp(find(kk,1,'last')));
end
k_slack0 = find(slack_MW(kk) < TOL, 1);
if isempty(k_slack0)
    fprintf('slack never drops to zero (min %.4f MW)\n', min(slack_MW(kk)));
else
    fprintf('slack drops to ~0 at avail = %.2f\n', A(k_slack0));
end
fprintf('max curtailment = %.3f MW (%.1f%% of available), Vmax then = %.4f\n', ...
    max(pv_curt(kk)), 100 * max(pv_curt(kk)) / max(pv_avail(kk)), ...
    vmax_pu(find(pv_curt == max(pv_curt(kk)), 1)));

%% ---------- 存数据 ----------
tag = sprintf('pen%03.0f_vmax%03.0f%s%s', opt.pv_penetration*100, opt.vmax*100, ...
    tern(opt.storage,'_stor',''), tern(opt.allow_export,'_exp',''));
if opt.save
    outdir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results', 'data');
    if ~exist(outdir, 'dir'), mkdir(outdir); end
    out = fullfile(outdir, sprintf('ed_sweep_%s.csv', tag));
    writetable(T, out);
    fprintf('data saved: %s\n', out);
end

%% ---------- 出图 ----------
if opt.plot && any(kk)
    f = figure('Position', [80 80 1450 430], 'Color', 'w', 'Visible', 'off');

    % (a) 调度堆叠图：slack + 光伏实发，并画出可用光伏与弃光
    subplot(1,3,1);
    xx = A(kk);
    area(xx, [slack_MW(kk), pv_disp(kk)], 'LineWidth', 1);
    hold on;
    plot(xx, pv_avail(kk), 'k--', 'LineWidth', 1.6);
    gapx = [xx; flipud(xx)];
    gapy = [pv_avail(kk); flipud(pv_disp(kk) + slack_MW(kk))];
    fill(gapx, gapy, [0.85 0.33 0.10], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    xlabel('PV availability  [fraction of installed]');
    ylabel('Power  [MW]');
    title(sprintf('(a) Dispatch stack (PV pen %.0f%%)', opt.pv_penetration*100));
    legend({'slack (grid import)', 'PV dispatched', 'PV available', 'curtailed'}, ...
        'Location', 'northwest', 'FontSize', 8);
    grid on; box on; xlim([0 1]);

    % (b) 成本 与 弃光
    subplot(1,3,2);
    yyaxis left
    plot(xx, cost(kk), '-', 'LineWidth', 2, 'Color', [0.12 0.47 0.71]);
    ylabel('Cost  [$/h]');
    yyaxis right
    plot(xx, pv_curt(kk), '-', 'LineWidth', 2, 'Color', [0.85 0.33 0.10]);
    ylabel('Curtailed PV  [MW]');
    xlabel('PV availability  [fraction of installed]');
    title('(b) Cost vs curtailment');
    grid on; box on; xlim([0 1]);
    legend({'cost', 'curtailment'}, 'Location', 'east', 'FontSize', 8);

    % (c) 电压
    subplot(1,3,3);
    plot(xx, vmax_pu(kk), '-', 'LineWidth', 2, 'Color', [0.84 0.15 0.16]); hold on;
    plot(xx, vmin_pu(kk), '-', 'LineWidth', 2, 'Color', [0.12 0.47 0.71]);
    yline(opt.vmax, '--', sprintf('Vmax %.2f', opt.vmax), 'Color', [0.84 0.15 0.16], ...
        'FontSize', 8, 'LabelHorizontalAlignment', 'left');
    yline(opt.vmin, ':', sprintf('Vmin %.2f', opt.vmin), 'Color', [0.12 0.47 0.71], ...
        'FontSize', 8, 'LabelHorizontalAlignment', 'left');
    xlabel('PV availability  [fraction of installed]');
    ylabel('Bus voltage  [p.u.]');
    title('(c) Voltage (OPF respects limits)');
    grid on; box on; xlim([0 1]);
    legend({'Vmax', 'Vmin'}, 'Location', 'southeast', 'FontSize', 8);

    sgtitle(sprintf('Deterministic economic dispatch : AC OPF, Vmax=%.2f, storage=%d', ...
        opt.vmax, opt.storage), 'FontSize', 13);

    outdir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results', 'figures');
    if ~exist(outdir, 'dir'), mkdir(outdir); end
    outfig = fullfile(outdir, sprintf('ed_sweep_%s.png', tag));
    exportgraphics(f, outfig, 'Resolution', 140);
    close(f);
    fprintf('figure saved: %s\n', outfig);
end
end


function s = tern(c, a, b)
if c, s = a; else, s = b; end
end
