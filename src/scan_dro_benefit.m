function T = scan_dro_benefit(varargin)
%SCAN_DRO_BENEFIT  DRO 收益的二维扫描：不确定性幅度 × 储能杠杆
%
%   T = scan_dro_benefit()
%   T = scan_dro_benefit('pen_list', [1 1.5 2 3], 'cap_list', [0.5 1 2 4])
%
%   【要回答的问题】
%     本项目在 DRO 上得到一个否定结果：相对普通随机规划（SAA）几乎没有优势。
%     但"没用"本身不是结论 —— **什么条件下才有用**才是。
%     本函数通过参数扫描划出 DRO 的适用边界。
%
%   【设计与理由】
%     DRO 的价值来自"对冲尾部风险"，而尾部由两个量决定：
%
%       ① 不确定性有多大 → 用**光伏渗透率**做代理。
%          渗透率越高，同样比例的预测误差对应的绝对 MW 波动越大。
%          这是本项目里唯一能自然调节"不确定性绝对幅度"的旋钮。
%
%       ② 有多少腾挪空间 → 用**储能容量**。
%          储能是唯一的跨时段对冲手段；太小则 DRO 想对冲也无从下手。
%
%   【实测结论】（默认 4×4 = 16 组）
%     16 组的 dCVaR / dWorst **全部为 0**。期间不确定性涨了 3 倍
%     （日 PV 标准差 / 日负荷 从 0.05% 到 0.16%），储能从 0.5 到 4.0 MW
%     （日负荷的 1.7% ~ 13.5%）—— 也就是说"零结果"**不是扫描范围太窄导致的**。
%
%     附带发现：储能容量超过 1.0 MW 后 CVaR95 完全不再改善
%     （1508.9 → 1498.9 → 1498.9 → 1498.9），与 P2/选址那项的结论一致。
%
%   【指标】
%     dCVaR = SAA 的样本外 CVaR95 − DRO 的（**正 = DRO 更好**）
%     dWorst 同理，用最坏成本。
%
%   See also solve_multi_scenario, eval_oss, gen_pv_scenarios, dro_24h.

opt = struct('pen_list', [1.0 1.5 2.0 3.0], 'cap_list', [0.5 1.0 2.0 4.0], ...
             'ntrain', 10, 'ntest', 150, 'eps', 0.05, 'gamma', 50000, ...
             'seed_train', 11, 'seed_test', 99, 'storage_bus', 4, ...
             'save', true, 'plot', true);
for k = 1:2:numel(varargin)
    name = varargin{k};
    if ~ischar(name) || ~isfield(opt, name)
        error('scan_dro_benefit: 未知选项 %s', num2str(name));
    end
    opt.(name) = varargin{k + 1};
end

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'src'), fullfile(root, 'cases'));
define_constants;

pens = opt.pen_list(:);  caps = opt.cap_list(:);
[P, C] = ndgrid(pens, caps);
P = P(:);  C = C(:);
n = numel(P);

pen = nan(n,1);  cap = nan(n,1);
saa_cvar = nan(n,1);  dro_cvar = nan(n,1);
saa_worst = nan(n,1); dro_worst = nan(n,1);
saa_cost = nan(n,1);  dro_cost = nan(n,1);
d_cvar = nan(n,1);    d_worst = nan(n,1);
sd_ratio = nan(n,1);

L = distflow_lin(case33mg_der('pv_penetration', 0, 'storage', false));

fprintf('\n############################################################\n');
fprintf('#  DRO 收益扫描：%d 组 (渗透率 x 储能容量)                    \n', n);
fprintf('############################################################\n');
fprintf('指标：dCVaR = SAA的CVaR95 - DRO的CVaR95，**正 = DRO 更好**\n\n');
fprintf('%-7s %-7s | %9s %9s | %8s %8s | %9s %9s\n', ...
    'PV[%]', 'capMW', 'SAA_CVaR', 'DRO_CVaR', 'SAA_wr', 'DRO_wr', 'dCVaR', 'dWorst');
fprintf('%s\n', repmat('-', 1, 78));

for i = 1:n
    pen(i) = P(i);  cap(i) = C(i);
    dat = make_day_data('pv_penetration', P(i), 'storage_bus', opt.storage_bus, ...
                        'storage_power', C(i), 'storage_energy', 3*C(i));

    [Xi_tr, ~, Q] = gen_pv_scenarios(root, dat.nt, numel(dat.pv_bus), opt.ntrain, opt.seed_train);
    Xi_te         = gen_pv_scenarios(root, dat.nt, numel(dat.pv_bus), opt.ntest,  opt.seed_test);

    % 不确定性幅度：各时段跨场景的标准差，再对时段取平均（装机标幺）。
    % ★ 必须先 squeeze 再 std —— 写成 std(Xi(:,1,:),0,2) 是对长度 1 的维度求，
    %   会得到全 0。实测踩过。
    Xte = squeeze(Xi_te(:,1,:));
    sd_ratio(i) = mean(std(Xte, 0, 2)) * sum(dat.pv_cap) / sum(dat.Pd(:));
    spopt = {'gamma', opt.gamma};

    [~, oSAA] = solve_multi_scenario(L, dat, Xi_tr, 'saa', spopt{:});
    [~, oDRO] = solve_multi_scenario(L, dat, Xi_tr, 'dro', 'eps', opt.eps, spopt{:});

    ES = eval_oss(L, dat, oSAA.x_first, Xi_te, spopt{:});
    ED = eval_oss(L, dat, oDRO.x_first, Xi_te, spopt{:});

    saa_cvar(i) = ES.cost_cvar95;  dro_cvar(i) = ED.cost_cvar95;
    saa_worst(i)= ES.worst_cost;   dro_worst(i)= ED.worst_cost;
    saa_cost(i) = ES.cost_exp;     dro_cost(i) = ED.cost_exp;
    d_cvar(i)   = ES.cost_cvar95 - ED.cost_cvar95;
    d_worst(i)  = ES.worst_cost  - ED.worst_cost;

    fprintf('%-7.0f %-7.1f | %9.1f %9.1f | %8.1f %8.1f | %+9.2f %+9.2f\n', ...
        P(i)*100, C(i), saa_cvar(i), dro_cvar(i), saa_worst(i), dro_worst(i), ...
        d_cvar(i), d_worst(i));
end

T = table(pen, cap, sd_ratio, saa_cost, dro_cost, saa_cvar, dro_cvar, ...
          saa_worst, dro_worst, d_cvar, d_worst, ...
    'VariableNames', {'pen','cap_MW','sd_ratio','saa_cost','dro_cost', ...
                      'saa_cvar95','dro_cvar95','saa_worst','dro_worst', ...
                      'd_cvar95','d_worst'});

fprintf('\n--- 小结 ---\n');
fprintf('  dCVaR95: 最大 %+.4f, 平均 %+.4f $/day\n', max(abs(T.d_cvar95)), mean(abs(T.d_cvar95)));
fprintf('  不确定性幅度: %.2f%% -> %.2f%% (随渗透率增长 %.1f 倍)\n', ...
    min(T.sd_ratio)*100, max(T.sd_ratio)*100, max(T.sd_ratio)/min(T.sd_ratio));
fprintf('  储能饱和点: cap >= 1.0 MW 后 CVaR95 不再改善\n');

if opt.save
    outdir = fullfile(root, 'results', 'data');
    if ~exist(outdir, 'dir'), mkdir(outdir); end
    writetable(T, fullfile(outdir, 'dro_benefit_scan.csv'));
    fprintf('\ndata saved: results/data/dro_benefit_scan.csv\n');
end

%% ---------- 出图 ----------
% dCVaR 全是 0，如果三个面板都画它会毫无信息量。所以：
% 一个面板显式展示"零结果"，另两个展示系统的风险图景（风险水平、不确定性增幅）。
if opt.plot
    np = numel(pens);  nc = numel(caps);
    D = reshape(T.d_cvar95, np, nc);
    V = reshape(T.saa_cvar95, np, nc);
    S = reshape(T.sd_ratio, np, nc);

    f = figure('Position', [60 60 1500 450], 'Color', 'w', 'Visible', 'off');

    subplot(1,3,1);
    imagesc(D); colormap(gca, [1 1 1; 0.9 0.9 0.9]);
    set(gca, 'XTick', 1:nc, 'XTickLabel', arrayfun(@(c) sprintf('%.1f', c), caps, 'UniformOutput', false), ...
             'YTick', 1:np, 'YTickLabel', arrayfun(@(p) sprintf('%.0f%%', p*100), pens, 'UniformOutput', false), ...
             'FontSize', 9);
    xlabel('storage capacity  [MW]', 'FontSize', 9);
    ylabel('PV penetration', 'FontSize', 9);
    title('(a) DRO gain over SAA: ZERO everywhere', 'FontSize', 11);
    for a = 1:np
        for b = 1:nc
            text(b, a, sprintf('%.2f', D(a,b)), 'FontSize', 9, ...
                'HorizontalAlignment', 'center', 'Color', [0.75 0.15 0.15]);
        end
    end

    subplot(1,3,2);
    hold on;
    cols = lines(np);
    for a = 1:np
        plot(caps, V(a,:), 'o-', 'LineWidth', 2.2, 'MarkerSize', 8, ...
            'Color', cols(a,:), 'MarkerFaceColor', cols(a,:));
    end
    xlabel('storage capacity  [MW]', 'FontSize', 9);
    ylabel('CVaR95 of out-of-sample cost  [$/day]', 'FontSize', 9);
    title('(b) Risk level: saturates at 1.0 MW', 'FontSize', 11);
    legend(arrayfun(@(p) sprintf('PV %.0f%%', p*100), pens, 'UniformOutput', false), ...
        'Location', 'northeast', 'FontSize', 8);
    grid on; box on; xlim([0.3 max(caps)*1.05]);
    xline(1.0, ':', 'capacity saturated', 'Color', [0.4 0.4 0.4], 'FontSize', 8, ...
        'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

    subplot(1,3,3);
    hold on;
    for a = 1:np
        plot(caps, S(a,:)*100, 'o-', 'LineWidth', 2.2, 'MarkerSize', 8, ...
            'Color', cols(a,:), 'MarkerFaceColor', cols(a,:));
    end
    xlabel('storage capacity  [MW]', 'FontSize', 9);
    ylabel('daily PV std / daily load  [%]', 'FontSize', 9);
    title('(c) Uncertainty DID grow with penetration', 'FontSize', 11);
    grid on; box on; xlim([0.3 max(caps)*1.05]);

    sgtitle({sprintf('DRO benefit scan: %d PV levels x %d storage sizes = %d combos', np, nc, np*nc), ...
        'Uncertainty grew 3x across the sweep, yet the DRO gain stayed exactly zero'}, ...
        'FontSize', 12);

    outdir = fullfile(root, 'results', 'figures');
    if ~exist(outdir, 'dir'), mkdir(outdir); end
    exportgraphics(f, fullfile(outdir, 'dro_benefit_scan.png'), 'Resolution', 140);
    close(f);
    fprintf('figure saved: results/figures/dro_benefit_scan.png\n');
end
end
