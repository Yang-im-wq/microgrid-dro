function [mdo, T, Tac] = run_most_24h(varargin)
%RUN_MOST_24H  24 小时 MOST 时序调度（光伏 + 储能）+ MATPOWER AC 校验
%
%   [mdo, T, Tac] = run_most_24h()
%   [...] = run_most_24h('pv_penetration', 1.0)
%
%   这是本项目里 MOST 第一次真正发挥作用 —— 储能的 SOC 是跨时段耦合，
%   单时段 OPF 表达不了（P3 里我们不得不把储能关掉）。
%
%   ★ 为什么必须有 AC 校验这一步
%     MOST 是**纯 DC** 的（见 most/README.md：implementation is limited to DC
%     power flow modeling）。它看不见电压，也看不见线路容量以外的任何物理约束。
%     而本算例真正的约束恰恰是电压（P2 实测：Vmax 1.05 在 172% 渗透率越限，
%     而热容量在现实渗透率下根本碰不到）。
%     所以流程是两段式：MOST 定调度曲线 → 逐时段回灌 MATPOWER 跑全交流潮流，
%     检查电压是否越限。**MOST 的解只有在 AC 校验下才成立。**
%
%   选项：
%     'pv_penetration'  光伏装机渗透率   (默认 1.5)
%     'vmax'/'vmin'     电压限值         (默认 1.05 / 0.90)
%     'cyclic'          储能日循环       (默认 true)
%     'save'/'plot'                      (默认 true / true)
%
%   See also build_most_data, most, loadmd, case33mg_der.

%% ---------- 选项 ----------
opt = struct('pv_penetration', 1.5, 'vmax', 1.05, 'vmin', 0.90, ...
             'cyclic', true, 'save', true, 'plot', true);
for k = 1:2:numel(varargin)
    name = varargin{k};
    if ~ischar(name) || ~isfield(opt, name)
        error('run_most_24h: 未知选项 ''%s''', num2str(name));
    end
    opt.(name) = varargin{k + 1};
end

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'src'), fullfile(root, 'cases'));
define_constants;

%% ---------- 1. 装配 MOST 输入 ----------
[mpc, xgd, sd, profiles, nt, meta] = build_most_data( ...
    'pv_penetration', opt.pv_penetration, 'vmax', opt.vmax, 'vmin', opt.vmin);

%% ---------- 2. MOST 求解 ----------
mpopt = mpoption('verbose', 0, 'out.all', 0);
mpopt = mpoption(mpopt, 'model', 'DC');
mpopt = mpoption(mpopt, 'most.dc_model', 1);          % DC 网络（MOST 仅支持 DC）
mpopt = mpoption(mpopt, 'most.solver', 'MIPS');
if opt.cyclic
    mpopt = mpoption(mpopt, 'most.storage.cyclic', 1); % 日循环：末时段 SOC = 初时段
end

mdi = loadmd(mpc, nt, xgd, sd, [], profiles);
mdo = most(mdi, mpopt);

if ~mdo.results.success          % 注意是 mdo.results.success，不是 mdo.success
    error('run_most_24h: MOST 未收敛');
end

EPg  = mdo.results.ExpectedDispatch;        % ng x nt  期望出力 (MW)
SoC  = mdo.Storage.ExpectedStorageState;    % ns x nt  SOC (MWh)
iess = mpc.iess;
ipv  = mpc.isolar;

%% ---------- 3. AC 校验（逐时段回灌全交流潮流）----------
% 注意：MOST 的 DC 解把储能当作"净注入"，但储能本身不产生电压支撑；
% 回灌时用的是 DC 解的 PG，电压由 AC 潮流重新算出来。
Tac = table((1:nt)', nan(nt,1), nan(nt,1), false(nt,1), nan(nt,1), ...
    'VariableNames', {'hour','vmin_pu','vmax_pu','converged','max_load_pct'});
for t = 1:nt
    m = mdo.flow(t, 1, 1).mpc;              % 该时段的完整算例（含 DC 解出的出力）
    r = runpf(m, mpoption('verbose', 0, 'out.all', 0));
    Tac.converged(t) = logical(r.success);
    if r.success
        Tac.vmin_pu(t) = min(r.bus(:, VM));
        Tac.vmax_pu(t) = max(r.bus(:, VM));
        rate = m.branch(1, RATE_A);
        if rate > 0
            Tac.max_load_pct(t) = max(hypot(r.branch(:,PF), r.branch(:,PT))) / rate * 100;
        end
    end
end
n_over = sum(Tac.vmax_pu > opt.vmax + 1e-6);
n_under = sum(Tac.vmin_pu < opt.vmin - 1e-6);

%% ---------- 4. 结果表 ----------
T = table((1:nt)', profiles(meta.i_pv_profile).values(:,1,1), ...
    profiles(meta.i_load_profile).values(:,1,1), ...
    sum(EPg(ipv,:), 1)', EPg(1,:)', EPg(iess,:)', SoC(1,:)', ...
    Tac.vmin_pu, Tac.vmax_pu, Tac.max_load_pct, ...
    'VariableNames', {'hour','pv_avail_pu','load_pu','pv_MW','slack_MW', ...
                      'storage_MW','soc_MWh','vmin_pu','vmax_pu','max_load_pct'});

%% ---------- 5. 摘要 ----------
fprintf('\n=== MOST 24h scheduling (PV + storage) ===\n');
fprintf('PV penetration = %.0f%%   storage = %.1f MWh / %.1f MW   cyclic = %d\n', ...
    opt.pv_penetration*100, meta.storage_MWh, meta.storage_MW, opt.cyclic);
fprintf('PV forecast = REAL (q50 of LSTM quantile) ; load = %s\n', meta.load_source);
fprintf('\nobjective = %.2f $/day\n', mdo.results.f);
fprintf('slack energy = %.3f MWh/day\n', sum(EPg(1,:)));
fprintf('PV dispatched = %.3f MWh/day  (available %.3f MWh/day)\n', ...
    sum(EPg(ipv,:), 'all'), sum(profiles(meta.i_pv_profile).values(:,1,1)) * meta.pv_cap_total_MW);
fprintf('storage throughput = %.3f MWh, SOC range [%.2f, %.2f] MWh\n', ...
    sum(abs(EPg(iess,:))), min(SoC), max(SoC));
fprintf('\n--- AC validation (feed MOST dispatch into full AC power flow) ---\n');
fprintf('converged           : %d / %d periods\n', sum(Tac.converged), nt);
fprintf('Vmax over all periods: %.4f  (limit %.2f, %d periods exceed)\n', ...
    max(Tac.vmax_pu), opt.vmax, n_over);
fprintf('Vmin over all periods: %.4f  (limit %.2f, %d periods violate)\n', ...
    min(Tac.vmin_pu), opt.vmin, n_under);
fprintf('max branch loading  : %.1f%%\n', max(Tac.max_load_pct));
if n_over > 0 || n_under > 0
    fprintf(['\n>>> MOST 给出的是 DC 最优解，但回灌全交流潮流后有 %d 个时段越上限、\n' ...
             '>>> %d 个时段越下限 —— 这就是"DC 看不见电压"的直接后果。\n' ...
             '>>> MOST 的解不能直接采用，必须经过 AC 校验/修正。<<<\n'], n_over, n_under);
else
    fprintf('\n>>> 全部时段电压合规：MOST 的 DC 解在 AC 校验下成立。<<<\n');
end

%% ---------- 6. 存数据 ----------
tag = sprintf('pen%03.0f_vmax%03.0f', opt.pv_penetration*100, opt.vmax*100);
if opt.save
    outdir = fullfile(root, 'results', 'data');
    if ~exist(outdir, 'dir'), mkdir(outdir); end
    writetable(T, fullfile(outdir, sprintf('most_24h_%s.csv', tag)));
    fid = fopen(fullfile(outdir, sprintf('most_24h_%s.meta.json', tag)), 'w');
    fprintf(fid, '%s', jsonencode(meta));
    fclose(fid);
    fprintf('\ndata saved: results/data/most_24h_%s.csv\n', tag);
end

%% ---------- 7. 出图 ----------
if opt.plot
    f = figure('Position', [70 70 1500 460], 'Color', 'w', 'Visible', 'off');
    hr = (1:nt)';

    subplot(1,3,1);
    % 储能画净注入：正 = 放电（供负荷），负 = 充电（消耗）。
    % 早先写成 max(-Pg,0) 并标成 "discharge" 是符号搞反了 —— 那其实是充电量。
    stack = [EPg(1,:)', sum(EPg(ipv,:),1)', EPg(iess,:)'];
    area(hr, stack, 'LineWidth', 0.9); hold on;
    plot(hr, profiles(meta.i_pv_profile).values(:,1,1) * meta.pv_cap_total_MW, 'k--', 'LineWidth', 1.6);
    plot(hr, profiles(meta.i_load_profile).values(:,1,1) * meta.load_mean_MW, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2);
    xlabel('hour'); ylabel('Power [MW]');
    title('(a) MOST dispatch (DC)', 'FontSize', 11);
    legend({'slack import','PV','storage net (+dis / -chg)','PV available','load'}, ...
        'Location', 'northwest', 'FontSize', 7);
    grid on; box on; xlim([1 nt]);

    subplot(1,3,2);
    plot(hr, SoC, '-o', 'LineWidth', 2, 'Color', [0.05 0.55 0.75], 'MarkerSize', 4);
    hold on;
    yline(min(sd.MinStorageLevel), ':', 'SOC min', 'FontSize', 8, 'LabelHorizontalAlignment','left');
    yline(max(sd.MaxStorageLevel), ':', 'SOC max', 'FontSize', 8, 'LabelHorizontalAlignment','left');
    xlabel('hour'); ylabel('State of charge [MWh]');
    title(sprintf('(b) Storage SOC  (%.1f MWh / %.1f MW)', meta.storage_MWh, meta.storage_MW), ...
        'FontSize', 11);
    grid on; box on; xlim([1 nt]);

    subplot(1,3,3);
    plot(hr, Tac.vmax_pu, '-o', 'LineWidth', 2, 'Color', [0.84 0.15 0.16], 'MarkerSize', 4); hold on;
    plot(hr, Tac.vmin_pu, '-o', 'LineWidth', 2, 'Color', [0.12 0.47 0.71], 'MarkerSize', 4);
    yline(opt.vmax, '--', sprintf('Vmax %.2f', opt.vmax), 'Color', [0.84 0.15 0.16], 'FontSize', 8);
    yline(opt.vmin, ':', sprintf('Vmin %.2f', opt.vmin), 'Color', [0.12 0.47 0.71], 'FontSize', 8);
    xlabel('hour'); ylabel('Bus voltage [p.u.]');
    title(sprintf('(c) AC validation: %d over Vmax, %d under Vmin', n_over, n_under), ...
        'FontSize', 11);
    legend({'Vmax','Vmin'}, 'Location', 'southeast', 'FontSize', 8);
    grid on; box on; xlim([1 nt]);

    sgtitle(sprintf(['MOST 24h scheduling + AC validation   ' ...
        '(PV %.0f%%, Vmax %.2f, storage %.1f MWh)'], ...
        opt.pv_penetration*100, opt.vmax, meta.storage_MWh), 'FontSize', 12);

    outdir = fullfile(root, 'results', 'figures');
    if ~exist(outdir, 'dir'), mkdir(outdir); end
    outfig = fullfile(outdir, sprintf('most_24h_%s.png', tag));
    exportgraphics(f, outfig, 'Resolution', 140);
    close(f);
    fprintf('figure saved: %s\n', outfig);
end
end
