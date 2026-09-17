function R = dro_24h(varargin)
%DRO_24H  P5 主实验：确定性 vs 随机规划 vs 分布鲁棒（Wasserstein/L1 模糊集）
%
%   R = dro_24h()
%   R = dro_24h('eps', 0.05, 'ntrain', 12, 'ntest', 200)
%
%   流程：
%     1. 从 11 分位数的概率预测采样出训练场景 / 测试场景（互不重叠）
%     2. 三种策略各求一个日前储能调度：
%          DET  只用中位数 q50（确定性）
%          SAA  最小化训练场景的平均成本（随机规划 / 样本平均近似）
%          DRO  平均成本 + ε·(最坏 − 最好)（分布鲁棒）
%     3. 把三个调度放到**同一批没参与优化的测试场景**上评估：
%          期望成本 / 期望切负荷 / 95% CVaR / 越限概率 / 最坏成本
%
%   预期结论（也是 DRO 的标准卖点）：
%     DRO 的"样本内"成本略高于 SAA（保守的代价），
%     但在**样本外**的尾部指标（CVaR、最坏成本、切负荷、越限概率）上更好。
%
%   选项：
%     'eps'      鲁棒半径          (默认 0.02)
%     'ntrain'   训练场景数        (默认 12)
%     'ntest'    测试场景数        (默认 200)
%     'pen'      光伏渗透率        (默认 1.5)
%     'seed_train' / 'seed_test'   随机种子
%
%   See also solve_multi_scenario, eval_oss, gen_pv_scenarios.

opt = struct('eps', 0.02, 'ntrain', 12, 'ntest', 200, 'pen', 1.5, ...
             'seed_train', 11, 'seed_test', 99, 'save', true, 'plot', true, ...
             'gamma', 2000, 'voll', 500, 'tag', '');
for k = 1:2:numel(varargin)
    name = varargin{k};
    if ~ischar(name) || ~isfield(opt, name)
        error('dro_24h: 未知选项 ''%s''', num2str(name));
    end
    opt.(name) = varargin{k + 1};
end
spopt = {'gamma', opt.gamma, 'voll', opt.voll};
evopt = {'gamma', opt.gamma, 'voll', opt.voll};

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'src'), fullfile(root, 'cases'));
define_constants;

fprintf('\n############################################################\n');
fprintf('#  P5 : DET  vs  SAA  vs  DRO                               #\n');
fprintf('############################################################\n\n');

%% ---------- 1. 模型与数据 ----------
mpc0 = case33mg_der('pv_penetration', opt.pen, 'ramp', 'pmax');
L    = distflow_lin(mpc0);
dat  = make_day_data('pv_penetration', opt.pen);
nt = dat.nt;  npv = numel(dat.pv_bus);
fprintf('network %d buses / %d branches ; PV %d nodes, %.3f MW total ; load %.2f MWh/day\n', ...
    L.nb, L.nbr, npv, sum(dat.pv_cap), sum(dat.Pd(:)));
fprintf('substation limit %.2f MW ; storage %.2f MWh / %.2f MW\n\n', ...
    dat.gmax, dat.ess_emax, dat.ess_pmax);

%% ---------- 2. 场景 ----------
[Xi_tr, qlev, Q] = gen_pv_scenarios(root, nt, npv, opt.ntrain, opt.seed_train);
Xi_te            = gen_pv_scenarios(root, nt, npv, opt.ntest,  opt.seed_test);
% 确定性场景 = 中位数 q50
i50 = find(qlev == 0.5, 1);
Xi_det = repmat(Q(:,i50), [1 npv 1]);

fprintf('scenarios: train %d, test %d (independent draws)\n', opt.ntrain, opt.ntest);
fprintf('PV daily energy per scenario: q05 %.1f / q50 %.1f / q95 %.1f MWh\n\n', ...
    sum(Q(:,1))*sum(dat.pv_cap), sum(Q(:,i50))*sum(dat.pv_cap), sum(Q(:,end))*sum(dat.pv_cap));

%% ---------- 3. 三种策略 ----------
names = {'DET', 'SAA', 'DRO'};
X = cell(3,1); info = cell(3,1);
t0 = tic;
[x1, i1] = solve_multi_scenario(L, dat, Xi_det, 'det', spopt{:});
fprintf('%-4s solved  obj = %10.2f   (%.1fs)\n', 'DET', i1.obj, toc(t0));
t0 = tic;
[x2, i2] = solve_multi_scenario(L, dat, Xi_tr, 'saa', spopt{:});
fprintf('%-4s solved  obj = %10.2f   (%.1fs)\n', 'SAA', i2.obj, toc(t0));
t0 = tic;
[x3, i3] = solve_multi_scenario(L, dat, Xi_tr, 'dro', 'eps', opt.eps, spopt{:});
fprintf('%-4s solved  obj = %10.2f   (eps = %.3f)  (%.1fs)\n\n', 'DRO', i3.obj, opt.eps, toc(t0));
X = {x1, x2, x3};  info = {i1, i2, i3};

%% ---------- 4. 样本外评估 ----------
fprintf('--- out-of-sample evaluation (%d unseen scenarios) ---\n', opt.ntest);
fprintf('%-5s | %11s | %9s | %11s | %9s | %8s | %11s\n', ...
    'plan', 'E[cost]$/d', 'E[shed]MWh', 'CVaR95$/d', 'worst$/d', 'viol%', 'E[Vmin]');
fprintf('%s\n', repmat('-', 1, 80));
E = cell(3,1);
for k = 1:3
    E{k} = eval_oss(L, dat, info{k}.x_first, Xi_te, evopt{:});
    fprintf('%-5s | %11.1f | %9.3f | %11.1f | %9.1f | %7.1f%% | %11.4f\n', ...
        names{k}, E{k}.cost_exp, E{k}.shed_exp, E{k}.cost_cvar95, ...
        E{k}.worst_cost, 100*E{k}.viol_prob, E{k}.vmin_exp);
end

%% ---------- 5. 汇总 ----------
R = struct('names', {names}, 'E', {E}, 'info', {info}, 'Xi_tr', Xi_tr, ...
           'Xi_te', Xi_te, 'dat', dat, 'eps', opt.eps, 'pen', opt.pen);

fprintf('\n--- 相对 SAA 的差异（正 = 更差）---\n');
fprintf('%-5s | %11s | %11s | %11s | %11s\n', ...
    'plan', 'dE[cost]%', 'dCVaR95%', 'dE[shed]%', 'dWorst%');
fprintf('%s\n', repmat('-', 1, 62));
for k = [1 3]
    fprintf('%-5s | %+11.2f | %+11.2f | %+11.2f | %+11.2f\n', names{k}, ...
        100*(E{k}.cost_exp   - E{2}.cost_exp)  / E{2}.cost_exp, ...
        100*(E{k}.cost_cvar95- E{2}.cost_cvar95)/ E{2}.cost_cvar95, ...
        100*(E{k}.shed_exp   - E{2}.shed_exp)  / max(E{2}.shed_exp, 1e-9), ...
        100*(E{k}.worst_cost - E{2}.worst_cost)/ E{2}.worst_cost);
end

%% ---------- 6. 存数据 ----------
tag = sprintf('pen%03.0f_eps%03.0f_g%05.0f%s', opt.pen*100, opt.eps*1000, opt.gamma, opt.tag);
if opt.save
    outdir = fullfile(root, 'results', 'data');
    if ~exist(outdir, 'dir'), mkdir(outdir); end
    T = table(names', ...
        cellfun(@(e) e.cost_exp,    E), cellfun(@(e) e.shed_exp, E), ...
        cellfun(@(e) e.cost_cvar95, E), cellfun(@(e) e.worst_cost, E), ...
        cellfun(@(e) e.viol_prob,   E), cellfun(@(e) e.vmin_exp, E), ...
        'VariableNames', {'plan','E_cost','E_shed_MWh','CVaR95','worst_cost', ...
                          'viol_prob','E_Vmin'});
    writetable(T, fullfile(outdir, sprintf('dro_compare_%s.csv', tag)));
    fprintf('\ndata saved: results/data/dro_compare_%s.csv\n', tag);
end

%% ---------- 7. 出图 ----------
if opt.plot
    f = figure('Position', [70 70 1480 460], 'Color', 'w', 'Visible', 'off');
    cols = [0.30 0.45 0.70; 0.30 0.60 0.30; 0.80 0.35 0.15];

    subplot(1,3,1);
    hold on;
    for k = 1:3
        histogram(E{k}.cost, 30, 'FaceColor', cols(k,:), 'FaceAlpha', 0.45, ...
            'EdgeColor', 'none', 'DisplayName', names{k});
    end
    xlabel('Out-of-sample cost  [$/day]'); ylabel('scenarios');
    title('(a) Cost distribution (unseen scenarios)', 'FontSize', 11);
    legend('Location', 'northwest', 'FontSize', 8); grid on; box on;

    subplot(1,3,2);
    M = [cellfun(@(e) e.cost_exp, E), cellfun(@(e) e.cost_cvar95, E), ...
         cellfun(@(e) e.worst_cost, E)];
    bar(M, 'grouped'); hold on;
    set(gca, 'XTickLabel', names, 'FontSize', 9);
    ylabel('Cost  [$/day]');
    title('(b) Expected / CVaR95 / worst', 'FontSize', 11);
    legend({'E[cost]','CVaR95','worst'}, 'Location', 'northwest', 'FontSize', 8);
    grid on; box on;

    subplot(1,3,3);
    hold on;
    for k = 1:3
        cdfplot_manual(E{k}.cost, cols(k,:), names{k});
    end
    xlabel('Out-of-sample cost  [$/day]'); ylabel('CDF');
    title('(c) Tail comparison', 'FontSize', 11);
    legend('Location', 'southeast', 'FontSize', 8); grid on; box on;

    sgtitle(sprintf(['DET vs SAA vs DRO   (PV %.0f%%, eps = %.3f, ' ...
        'storage %.1f MWh, %d test scenarios)'], ...
        opt.pen*100, opt.eps, dat.ess_emax, opt.ntest), 'FontSize', 12);

    outdir = fullfile(root, 'results', 'figures');
    if ~exist(outdir, 'dir'), mkdir(outdir); end
    outfig = fullfile(outdir, sprintf('dro_compare_%s.png', tag));
    exportgraphics(f, outfig, 'Resolution', 140);
    close(f);
    fprintf('figure saved: %s\n', outfig);
end
end


function cdfplot_manual(v, c, nm)
vs = sort(v(:));
y = (1:numel(vs))' / numel(vs);
plot(vs, y, '-', 'LineWidth', 2, 'Color', c, 'DisplayName', nm);
end
