function [Trad, Tmesh] = reconfigure(varargin)
%RECONFIGURE  IEEE 33 节点网络重构 / 微电网环网运行分析
%
%   [Trad, Tmesh] = reconfigure()
%   [Trad, Tmesh] = reconfigure('pv_penetration', 0.8)
%
%   选项：
%     'pv_penetration' 光伏渗透率            (默认 0，即经典重构基准)
%     'vmax' / 'vmin'  电压限值              (默认 1.05 / 0.90)
%     'save'           存 CSV                (默认 true)
%     'plot'           出图                  (默认 true)
%
%   做两件**不同**的事（容易混淆，务必分清）：
%
%   【A】径向重构 (radial reconfiguration)
%       经典配电网重构：**合上一个联络开关 + 断开环路上一个分段开关**，
%       保持网络仍是辐射状（33 节点 / 32 条闭合支路 / 无环）。
%       枚举方法：对每个联络开关，合上它形成唯一一个基本回路，
%       再穷举"断开该回路上哪条支路"，即可得到全部单步可达的辐射配置。
%       这是 Baran-Wu / Civanlar 一脉的标准做法。
%
%   【B】微电网环网运行 (meshed microgrid operation)
%       合上任意子集的联络开关，**允许环网**（2^5 = 32 种组合）。
%       此时网络并非辐射状，但对微电网而言这是合法的运行方式
%       —— 多电源互为备用、降低网损。case33mg 的 "mg" 正是指这个。
%
%   注意：MATPOWER 的前推回代算法 (pf.alg='PQSUM') **遇环直接报错**，
%   所以环网态必须用牛顿法 (默认 'NR')。

%% ---------- 选项 ----------
opt = struct('pv_penetration', 0, 'vmax', 1.05, 'vmin', 0.90, ...
             'save', true, 'plot', true);
for k = 1:2:numel(varargin)
    name = varargin{k};
    if ~ischar(name) || ~isfield(opt, name)
        error('reconfigure: 未知选项 ''%s''', num2str(name));
    end
    opt.(name) = varargin{k + 1};
end

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'cases'));
define_constants;
mpopt = mpoption('verbose', 0, 'out.all', 0);

% 基准算例（可带光伏）
mpc0 = case33mg_der('pv_penetration', opt.pv_penetration, ...
                    'vmax', opt.vmax, 'vmin', opt.vmin);
nb = size(mpc0.bus, 1);

tie_idx  = find(mpc0.branch(:, BR_STATUS) == 0);    % 常开联络开关（5 条）
seci_idx = find(mpc0.branch(:, BR_STATUS) == 1);    % 常闭分段开关（32 条）
fprintf('network: %d buses, %d branches (%d closed, %d tie open)\n', ...
    nb, size(mpc0.branch,1), numel(seci_idx), numel(tie_idx));
fprintf('PV penetration = %.0f%%, Vmax = %.2f, Vmin = %.2f\n\n', ...
    opt.pv_penetration*100, opt.vmax, opt.vmin);

%% ==================== A. 径向重构（分支交换局部搜索） ====================
% 单步换开关只能覆盖"基态一步可达"的邻域，够不到经典最优（139.55 kW 需要同时
% 换掉多个开关）。所以用标准的 branch-exchange 局部搜索：反复"合上一条开着的
% 支路 -> 断开其回路上一条闭着的支路"，每步取改进最大的一次，直到无法改进。
fprintf('=== A. radial reconfiguration (branch-exchange local search) ===\n');

nbr = size(mpc0.branch, 1);
open_set   = tie_idx(:);
closed_set = setdiff((1:nbr)', open_set);

r0 = runpf_with_open(mpc0, open_set, mpopt);
cur_loss = sum(r0.branch(:,PF) + r0.branch(:,PT)) * 1000;
fprintf('  base: open %s -> loss = %.2f kW\n', mat2str(sort(open_set)'), cur_loss);

hist_open = {open_set};
hist_loss = cur_loss;
n_eval = 0;

for iter = 1:30
    best_L = Inf; best_open = [];
    for i = open_set'                       % 合上支路 i -> 形成唯一回路
        a = mpc0.branch(i, F_BUS);
        b = mpc0.branch(i, T_BUS);
        path = tree_path(mpc0.branch(:, F_BUS), mpc0.branch(:, T_BUS), ...
                         closed_set, nb, a, b);
        for j = path(:)'                    % 断开回路上的支路 j -> 重新辐射
            new_open = [setdiff(open_set, i); j];
            m2 = mpc0;
            m2.branch(:, BR_STATUS) = 1;
            m2.branch(new_open, BR_STATUS) = 0;
            if ~is_radial_connected(m2.branch(:,F_BUS), m2.branch(:,T_BUS), ...
                                    m2.branch(:,BR_STATUS), nb)
                continue
            end
            r2 = runpf(m2, mpopt);
            n_eval = n_eval + 1;
            if ~r2.success, continue; end
            L = sum(r2.branch(:,PF) + r2.branch(:,PT)) * 1000;
            if L < best_L
                best_L = L;
                best_open = new_open;
            end
        end
    end
    if isempty(best_open) || best_L >= cur_loss - 1e-9
        break                                % 无法改进 -> 收敛
    end
    open_set   = best_open;
    closed_set = setdiff((1:nbr)', open_set);
    cur_loss   = best_L;
    hist_open{end+1} = open_set;             %#ok<AGROW>
    hist_loss(end+1) = cur_loss;             %#ok<AGROW>
    fprintf('  iter %2d: loss -> %7.2f kW  (open: %s)\n', ...
        iter, cur_loss, mat2str(sort(open_set)'));
end
fprintf('  converged after %d improving steps (%d power flows evaluated)\n', ...
    numel(hist_loss) - 1, n_eval);
fprintf('  best radial config = %.2f kW  (saving %.1f%% vs base)\n', ...
    cur_loss, (1 - cur_loss/hist_loss(1)) * 100);

% 最优配置的详细指标
m_best = mpc0;
m_best.branch(:, BR_STATUS) = 1;
m_best.branch(open_set, BR_STATUS) = 0;
r_best = runpf(m_best, mpopt);
fprintf('    open branches: %s\n', mat2str(sort(open_set)'));
fprintf('    Vmin = %.4f @bus%d,  Vmax = %.4f,  max loading = %.1f%%\n', ...
    min(r_best.bus(:,VM)), find(r_best.bus(:,VM)==min(r_best.bus(:,VM)),1), ...
    max(r_best.bus(:,VM)), ...
    max(hypot(r_best.branch(:,PF), r_best.branch(:,PT))) / m_best.branch(1,RATE_A) * 100);

Trad = table((0:numel(hist_loss)-1)', hist_loss(:), ...
    cellfun(@(x) strjoin(string(sort(x)'), '|'), hist_open(:), 'UniformOutput', false), ...
    'VariableNames', {'step','loss_kW','open_branches'});
% 注意：支路列表用 '|' 分隔而不是空格 —— 空格分隔的单元格会让 MATLAB 的
% readtable 把分隔符误探测成空格，把 3 列读成 5 列（实测踩过）。
% 读取方也建议显式写 readtable(p, 'Delimiter', ',')。

%% ==================== B. 微电网环网运行 ====================
fprintf('\n=== B. meshed microgrid operation (32 tie-switch subsets) ===\n');
nt = numel(tie_idx);
Mmask = (0:2^nt-1)';
Mloss = nan(2^nt,1); Mvmax = nan(2^nt,1); Mvmin = nan(2^nt,1);
Mload = nan(2^nt,1); Mok = false(2^nt,1); Mntie = zeros(2^nt,1);
for k = 1:2^nt
    bits = bitget(k-1, 1:nt) > 0;
    m = mpc0;
    m.branch(tie_idx(bits), BR_STATUS) = 1;
    Mntie(k) = sum(bits);
    r = runpf(m, mpopt);
    Mok(k) = logical(r.success);
    if ~r.success, continue; end
    Mloss(k) = sum(r.branch(:,PF) + r.branch(:,PT)) * 1000;
    Mvmax(k) = max(r.bus(:,VM));
    Mvmin(k) = min(r.bus(:,VM));
    Smax = max(hypot(r.branch(:,PF), r.branch(:,PT)));
    rate = m.branch(1, RATE_A);
    if rate > 0, Mload(k) = Smax / rate * 100; else, Mload(k) = NaN; end
end
Mok(1) = true;
Mloss(1) = hist_loss(1); Mvmax(1) = max(r0.bus(:,VM)); Mvmin(1) = min(r0.bus(:,VM));
Mload(1) = max(hypot(r0.branch(:,PF), r0.branch(:,PT))) / mpc0.branch(1,RATE_A) * 100;

[best_m, km] = min(Mloss);
fprintf('  no tie closed (radial) = %.2f kW\n', Mloss(1));
fprintf('  best meshed            = %.2f kW  (saving %.1f%%, %d ties closed)\n', ...
    best_m, (1 - best_m/Mloss(1)) * 100, Mntie(km));
fprintf('  all 5 ties closed      = %.2f kW (saving %.1f%%)\n', ...
    Mloss(end), (1 - Mloss(end)/Mloss(1)) * 100);
fprintf('  Vmin: radial %.4f -> all-meshed %.4f\n', Mvmin(1), Mvmin(end));

Tmesh = table(Mmask, Mntie, Mok, Mloss, Mvmin, Mvmax, Mload, ...
    'VariableNames', {'mask','n_tie_closed','converged', ...
                      'loss_kW','vmin_pu','vmax_pu','max_load_pct'});

%% ==================== 存数据 ====================
tag = sprintf('pv%03.0f_vmax%03.0f', opt.pv_penetration*100, opt.vmax*100);
if opt.save
    outdir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results', 'data');
    if ~exist(outdir, 'dir'), mkdir(outdir); end
    writetable(Trad,  fullfile(outdir, sprintf('reconfig_radial_%s.csv', tag)));
    writetable(Tmesh, fullfile(outdir, sprintf('reconfig_meshed_%s.csv', tag)));
    fprintf('\ndata saved: results/data/reconfig_{radial,meshed}_%s.csv\n', tag);
end

%% ==================== 出图 ====================
if opt.plot
    f = figure('Position', [80 80 1250 440], 'Color', 'w', 'Visible', 'off');

    subplot(1,2,1);
    steps = (0:numel(hist_loss)-1)';
    plot(steps, hist_loss, 'o-', 'LineWidth', 1.6, 'MarkerSize', 5, ...
        'Color', [0.12 0.47 0.71], 'MarkerFaceColor', [0.12 0.47 0.71]);
    hold on;
    plot(steps(end), hist_loss(end), 'o', 'MarkerSize', 10, ...
        'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerEdgeColor', 'none');
    yline(hist_loss(1), '--', sprintf('base %.1f kW', hist_loss(1)), ...
        'Color', [0.5 0.5 0.5], 'FontSize', 8);
    xlabel('branch-exchange iteration');
    ylabel('Network loss  [kW]');
    title(sprintf('(a) Radial reconfiguration: %.1f -> %.1f kW (-%.1f%%)', ...
        hist_loss(1), hist_loss(end), (1 - hist_loss(end)/hist_loss(1))*100));
    grid on; box on;

    subplot(1,2,2);
    yyaxis left
    plot(Mntie, Mloss, 'o', 'MarkerSize', 7, 'MarkerFaceColor', [0.30 0.55 0.20], ...
        'MarkerEdgeColor', 'none');
    ylabel('Network loss  [kW]');
    yyaxis right
    plot(Mntie, Mvmin, 's', 'MarkerSize', 5, 'MarkerFaceColor', [0.84 0.15 0.16], ...
        'MarkerEdgeColor', 'none');
    ylabel('Vmin  [p.u.]');
    xlabel('number of tie switches closed (0 = radial)');
    title('(b) Meshed microgrid operation');
    grid on; box on;
    legend({'loss', 'Vmin'}, 'Location', 'northeast', 'FontSize', 9);

    sgtitle(sprintf('IEEE 33-bus : reconfiguration vs meshed operation  (PV = %.0f%%)', ...
        opt.pv_penetration*100), 'FontSize', 13);

    outdir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results', 'figures');
    if ~exist(outdir, 'dir'), mkdir(outdir); end
    outfig = fullfile(outdir, sprintf('reconfigure_%s.png', tag));
    exportgraphics(f, outfig, 'Resolution', 140);
    close(f);
    fprintf('figure saved: %s\n', outfig);
end
end


%% ==================== helpers ====================
function r = runpf_with_open(mpc0, open_set, mpopt)
% 把 open_set 里的支路置为断开，其余闭合，跑一次潮流
define_constants;
m = mpc0;
m.branch(:, BR_STATUS) = 1;
m.branch(open_set, BR_STATUS) = 0;
r = runpf(m, mpopt);
end


function path = tree_path(ef, et, eidx, nb, a, b)
% 在由 eidx 指定的边构成的树上，BFS 找 a->b 的路径，返回边在分支表中的行号（按路径顺序）
% ef/et 是**全部**支路的首末母线，eidx 是构成树的支路子集
adj = cell(nb, 1);
for k = eidx(:)'
    adj{ef(k)}(end+1) = k;
    adj{et(k)}(end+1) = k;
end
pred_edge = zeros(nb, 1);       % 到达该节点所经由的支路
visited   = false(nb, 1);
q = a; visited(a) = true;
while ~isempty(q)
    v = q(1); q(1) = [];
    if v == b, break; end
    for k = adj{v}
        if ef(k) == v, w = et(k); else, w = ef(k); end
        if ~visited(w)
            visited(w) = true; pred_edge(w) = k; q(end+1) = w;   %#ok<AGROW>
        end
    end
end
if ~visited(b)
    path = [];
    return
end
path = [];
v = b;
while v ~= a
    k = pred_edge(v);
    if k == 0, path = []; return; end
    path = [k path];                                         %#ok<AGROW>
    if ef(k) == v, v = et(k); else, v = ef(k); end
end
end


function ok = is_radial_connected(ef, et, status, nb)
% 辐射 + 连通：闭合支路数 == nb-1 且全图连通
ci = find(status);
if numel(ci) ~= nb - 1
    ok = false; return
end
adj = cell(nb, 1);
for k = ci(:)'
    adj{ef(k)}(end+1) = et(k);
    adj{et(k)}(end+1) = ef(k);
end
visited = false(nb, 1);
q = 1; visited(1) = true; cnt = 1;
while ~isempty(q)
    v = q(1); q(1) = [];
    for w = adj{v}
        if ~visited(w)
            visited(w) = true; cnt = cnt + 1; q(end+1) = w;      %#ok<AGROW>
        end
    end
end
ok = (cnt == nb);
end
