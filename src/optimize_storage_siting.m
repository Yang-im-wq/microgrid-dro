function T = optimize_storage_siting(varargin)
%OPTIMIZE_STORAGE_SITING  储能选址优化：扫描候选母线，量化接入位置的影响
%
%   T = optimize_storage_siting()
%   T = optimize_storage_siting('bus_list', 2:33, 'save', true)
%
%   背景：项目早期实测发现储能接在 bus 18（馈线最弱末端）时，
%   光是在那里充 1 MW 就会把电压压到 0.832 —— **同样的储能放在不同位置，
%   效果可能天差地别**。本函数就是把这个差异量化出来。
%
%   对每个候选母线 b：
%     1. 把储能改接到 b（通过 make_day_data 的 storage_bus 选项）
%     2. 用 LinDistFlow 调度模型求日前计划（单场景，q50 预测）
%     3. 把计划**回灌 MATPOWER 全交流潮流**逐时段校验
%     4. 记录：成本、切负荷、最低/最高电压、越限时段数、最大支路负载
%
%   选项：
%     'bus_list'   候选母线（默认 2:33，即全部带负荷的母线）
%     'pen'        光伏渗透率（默认 1.5）
%     'gmax'       变电站容量上限 MW（默认 4.2）
%     'save'/'plot'
%
%   ★ 指标怎么读
%     - 成本低 + 切负荷少   -> 经济性好
%     - 越限时段数少        -> 对电压友好
%     - 三者常常不一致，**这正是选址问题的核心张力**
%
%   See also make_day_data, opt_dispatch_lindistflow, validate_dispatch_ac.

opt = struct('bus_list', 2:33, 'pen', 1.5, 'gmax', 4.2, ...
             'vmax', 1.05, 'vmin', 0.90, 'save', true, 'plot', true);
for k = 1:2:numel(varargin)
    name = varargin{k};
    if ~ischar(name) || ~isfield(opt, name)
        error('optimize_storage_siting: 未知选项 ''%s''', num2str(name));
    end
    opt.(name) = varargin{k + 1};
end

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'src'), fullfile(root, 'cases'));
define_constants;
mpopt = mpoption('verbose', 0, 'out.all', 0);

buses = opt.bus_list(:);
nb_cand = numel(buses);
fprintf('\n############################################################\n');
fprintf('#  储能选址优化：扫描 %d 个候选母线                          \n', nb_cand);
fprintf('############################################################\n\n');

L = distflow_lin(case33mg_der('pv_penetration', 0, 'storage', false));

B     = nan(nb_cand,1);
cost  = nan(nb_cand,1);
shed  = nan(nb_cand,1);
curt  = nan(nb_cand,1);
vminA = nan(nb_cand,1);
vmaxA = nan(nb_cand,1);
nover = nan(nb_cand,1);
nunder= nan(nb_cand,1);
loadp = nan(nb_cand,1);
ok    = false(nb_cand,1);

fprintf('%-5s | %10s | %8s | %8s | %8s | %6s | %6s | %8s\n', ...
    'bus', 'cost$/d', 'shedMWh', 'VminAC', 'VmaxAC', '#over', '#under', 'maxLoad%');
fprintf('%s\n', repmat('-', 1, 76));

for i = 1:nb_cand
    b = buses(i);
    try
        dat = make_day_data('pv_penetration', opt.pen, 'gmax', opt.gmax, ...
                            'storage_bus', b, 'vmax', opt.vmax, 'vmin', opt.vmin);
    catch ME
        fprintf('%-5d | 算例构造失败: %s\n', b, ME.message);
        continue
    end
    avail = dat.pv_avail_max ./ max(dat.pv_cap);

    [x, o] = opt_dispatch_lindistflow(L, dat, avail, 'vmax', opt.vmax, 'vmin', opt.vmin);
    if o.exitflag ~= 1
        fprintf('%-5d | 优化不可行 (exitflag=%d)\n', b, o.exitflag);
        continue
    end

    % 回灌全交流潮流校验
    mpc0 = case33mg_der('pv_penetration', opt.pen, 'storage', true, ...
                        'storage_bus', b, 'ramp', 'pmax');
    Tac = validate_dispatch_ac(mpc0, L, dat, x, o, 'vmax', opt.vmax, ...
                               'vmin', opt.vmin, 'quiet', true);

    B(i)     = b;
    cost(i)  = o.fval;
    shed(i)  = sum(x(o.idx.LS));
    curt(i)  = sum(x(o.idx.U));
    vminA(i) = min(Tac.vmin_ac);
    vmaxA(i) = max(Tac.vmax_ac);
    nover(i) = sum(Tac.vmax_ac > opt.vmax + 1e-6);
    nunder(i)= sum(Tac.vmin_ac < opt.vmin - 1e-6);
    loadp(i) = max(Tac.load_pct);
    ok(i)    = true;

    fprintf('%-5d | %10.1f | %8.3f | %8.4f | %8.4f | %6d | %6d | %8.1f\n', ...
        b, cost(i), shed(i), vminA(i), vmaxA(i), nover(i), nunder(i), loadp(i));
end

T = table(B, cost, shed, curt, vminA, vmaxA, nover, nunder, loadp, ok, ...
    'VariableNames', {'bus','cost_usd_d','shed_MWh','curt_MWh','vmin_ac', ...
                      'vmax_ac','n_over','n_under','max_load_pct','ok'});
T = T(T.ok, :);

%% ---------- 排名 ----------
% 综合分：越限时段数是硬伤（权重高），然后是成本与切负荷
score = T.cost_usd_d + 1e4 * (T.n_over + T.n_under) + 500 * T.shed_MWh;
T.score = score;
T = sortrows(T, 'score');
T.rank = (1:height(T))';

fprintf('\n--- 排名（综合分 = 成本 + 1e4×越限时段 + 500×切负荷）---\n');
fprintf('%-5s %-5s | %10s | %8s | %6s | %6s | %9s\n', ...
    'rank','bus','cost$/d','shedMWh','#over','#under','maxLoad%');
fprintf('%s\n', repmat('-', 1, 64));
for i = 1:min(8, height(T))
    fprintf('%-5d %-5d | %10.1f | %8.3f | %6d | %6d | %9.1f\n', ...
        T.rank(i), T.bus(i), T.cost_usd_d(i), T.shed_MWh(i), ...
        T.n_over(i), T.n_under(i), T.max_load_pct(i));
end

fprintf('\n--- 关键对比 ---\n');
b_def = find(T.bus == 18, 1);
if ~isempty(b_def)
    fprintf('  当前默认位置 bus 18 : 成本 %.1f $/d，越限 %d 时段\n', ...
        T.cost_usd_d(b_def), T.n_over(b_def)+T.n_under(b_def));
end
fprintf('  最优位置     bus %-3d : 成本 %.1f $/d，越限 %d 时段\n', ...
    T.bus(1), T.cost_usd_d(1), T.n_over(1)+T.n_under(1));
fprintf('  成本最优     bus %-3d : 成本 %.1f $/d\n', ...
    T.bus(find(T.cost_usd_d==min(T.cost_usd_d),1)), min(T.cost_usd_d));
fprintf('  电压最优     bus %-3d : 越限 %d 时段\n', ...
    T.bus(find(T.n_over+T.n_under==min(T.n_over+T.n_under),1)), ...
    min(T.n_over+T.n_under));

%% ---------- 存数据 ----------
if opt.save
    outdir = fullfile(root, 'results', 'data');
    if ~exist(outdir, 'dir'), mkdir(outdir); end
    out = fullfile(outdir, sprintf('storage_siting_pen%03.0f_g%03.0f.csv', ...
        opt.pen*100, opt.gmax*10));
    writetable(T, out);
    fprintf('\ndata saved: %s\n', out);
end

%% ---------- 出图 ----------
if opt.plot
    f = figure('Position', [70 70 1480 460], 'Color', 'w', 'Visible', 'off');

    subplot(1,3,1);
    scatter(T.cost_usd_d, T.n_over + T.n_under, 60, T.max_load_pct, 'filled', ...
        'MarkerEdgeColor', [0.3 0.3 0.3]); hold on;
    cb = colorbar; cb.Label.String = 'max branch loading [%]';
    i18 = find(T.bus == 18, 1);  i1 = 1;
    if ~isempty(i18)
        plot(T.cost_usd_d(i18), T.n_over(i18)+T.n_under(i18), 'p', 'MarkerSize', 18, ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.95 0.85 0.20]);
        text(T.cost_usd_d(i18)+1, T.n_over(i18)+T.n_under(i18)+0.3, 'bus 18 (current)', 'FontSize', 7);
    end
    plot(T.cost_usd_d(i1), T.n_over(i1)+T.n_under(i1), 'o', 'MarkerSize', 12, ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.30 0.60 0.30]);
    text(T.cost_usd_d(i1)+1, T.n_over(i1)+T.n_under(i1)-0.5, sprintf('bus %d (best)', T.bus(i1)), 'FontSize', 7);
    xlabel('cost  [$/day]'); ylabel('violating periods  [count]');
    title('(a) Cost vs voltage violations', 'FontSize', 11);
    grid on; box on;

    subplot(1,3,2);
    % 用"越限时段数"而不是成本：各母线的**成本几乎一样**（差 0.4%），
    % 柱子看不出区别；真正有差异的是电压越限时段数。
    nv = T.n_over + T.n_under;
    nshow = min(16, height(T));
    barh(1:nshow, nv(1:nshow), 0.6, 'FaceColor', [0.84 0.15 0.16]);
    set(gca, 'YTick', 1:nshow, ...
        'YTickLabel', arrayfun(@(b) sprintf('bus %d', b), T.bus(1:nshow), 'UniformOutput', false), ...
        'YDir', 'reverse', 'FontSize', 7);
    xlabel('violating periods  [count]');
    title('(b) Voltage violations by bus (best 16)', 'FontSize', 11);
    grid on; box on;
    xlim([0 max([nv(1:nshow); 1])*1.15]);

    subplot(1,3,3);
    hold on;
    scatter(T.max_load_pct, T.shed_MWh, 60, T.n_over + T.n_under, 'filled', ...
        'MarkerEdgeColor', [0.3 0.3 0.3]);
    cb2 = colorbar; cb2.Label.String = 'violating periods';
    if ~isempty(i18)
        plot(T.max_load_pct(i18), T.shed_MWh(i18), 'p', 'MarkerSize', 18, ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.95 0.85 0.20]);
    end
    xlabel('max branch loading  [%]'); ylabel('load shed  [MWh]');
    title('(c) Thermal stress vs load shedding', 'FontSize', 11);
    grid on; box on;

    sgtitle(sprintf('Storage siting sweep over %d buses   (PV %.0f%%, gmax %.1f MW)', ...
        height(T), opt.pen*100, opt.gmax), 'FontSize', 12);

    outdir = fullfile(root, 'results', 'figures');
    if ~exist(outdir, 'dir'), mkdir(outdir); end
    outfig = fullfile(outdir, sprintf('storage_siting_pen%03.0f.png', opt.pen*100));
    exportgraphics(f, outfig, 'Resolution', 140);
    close(f);
    fprintf('figure saved: %s\n', outfig);
end
end
