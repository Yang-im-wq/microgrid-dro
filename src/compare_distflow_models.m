function T = compare_distflow_models(varargin)
%COMPARE_DISTFLOW_MODELS  潮流模型精度对比：LinDistFlow vs SOCP vs 全交流潮流
%
%   T = compare_distflow_models()
%   T = compare_distflow_models('pen_list', [0 0.8 1.5 2.5])
%
%   方法（**固定同一组注入**，只比电压预测）：
%     1. 对每个光伏渗透率跑一次 MATPOWER 全交流潮流 → 得到真实电压与真实注入
%     2. 把**完全相同的注入**分别喂给两个简化模型：
%          LinDistFlow —— 解析式（路径矩阵）
%          SOCP        —— 二阶锥模型的"潮流模式"（固定注入、极小化 Σl 让锥绷紧）
%     3. 比较三者的电压
%
%   ★ 为什么这样比
%     注入相同，差异就**纯粹来自电压模型本身**，不掺杂"优化目标不同"的干扰。
%
%     SOCP 的 DistFlow 方程与真实交流潮流**结构一致**（都含线损项 (r²+x²)l），
%     只是把 l = |I|² 松弛成二阶锥。Farivar & Low (2013) 证明对辐射网该松弛是紧的，
%     所以 SOCP 的电压应该几乎精确。
%     LinDistFlow 丢掉了线损项，压降系统性偏小、电压偏高 —— **对欠压乐观**。
%     这张表就是把"偏多少"量出来。
%
%   选项：'pen_list' 渗透率列表（默认 [0 0.5 1.0 1.5 2.0 2.5]）
%         'hour' 评估时段（默认 19，晚高峰）
%
%   See also opt_dispatch_socp, distflow_eval, distflow_lin.

opt = struct('pen_list', (0:0.5:2.5)', 'hour', 19, ...
             'vmax', 1.05, 'vmin', 0.90, 'save', true, 'plot', true);
for k = 1:2:numel(varargin)
    name = varargin{k};
    if ~ischar(name) || ~isfield(opt, name)
        error('compare_distflow_models: 未知选项 ''%s''', num2str(name));
    end
    opt.(name) = varargin{k + 1};
end

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'src'), fullfile(root, 'cases'));
define_constants;
mpopt = mpoption('verbose', 0, 'out.all', 0);

pens = opt.pen_list(:);  n = numel(pens);  h = opt.hour;

fprintf('\n############################################################\n');
fprintf('#  潮流模型精度对比：LinDistFlow  vs  SOCP  vs  全交流潮流    #\n');
fprintf('############################################################\n');
fprintf('评估时段 = %d 时（晚高峰）；三个模型用**完全相同的注入**\n\n', h);
fprintf('%-8s | %9s | %9s | %9s | %10s | %10s\n', ...
    'PV[%]', 'Vmin_AC', 'Vmin_Lin', 'Vmin_SOCP', 'err_Lin', 'err_SOCP');
fprintf('%s\n', repmat('-', 1, 68));

pen = nan(n,1);  vmin_ac = nan(n,1);  vmin_lin = nan(n,1);  vmin_socp = nan(n,1);
vmax_ac = nan(n,1);  vmax_lin = nan(n,1);  vmax_socp = nan(n,1);
err_lin = nan(n,1);  err_socp = nan(n,1);  ncones = nan(n,1);

for i = 1:n
    pen(i) = pens(i);
    dat = make_day_data('pv_penetration', pens(i));
    L   = distflow_lin(case33mg_der('pv_penetration', 0, 'storage', false));
    nb  = L.nb;  npv = numel(dat.pv_bus);

    % ★ 交流潮流必须用**与该时段相同的负荷**。
    %   最初我直接 runpf(m) 用的是基准负荷，而 LinDistFlow 那边用的是 h 时的负荷，
    %   两者根本不是同一个运行点 —— 比较出来的误差没有意义（还会符号相反）。
    %   实测踩过：当时 LinDistFlow 的误差是负的，与之前 P5 的正号结论矛盾，就是这条线索。
    m = case33mg_der('pv_penetration', pens(i), 'storage', false, 'ramp', 'pmax');
    m.bus(:, PD) = dat.Pd(:, h);
    m.bus(:, QD) = dat.Qd;
    r = runpf(m, mpopt);
    if ~r.success, continue; end

    % ---- 从交流解得真实注入 ----
    if isfield(m, 'isolar')
        pv_disp = m.gen(m.isolar, PMAX);        % 各光伏满发（未弃光）
    else
        pv_disp = [];                           % 渗透率为 0 时没有光伏
    end
    g_ac    = r.gen(1, PG);                     % 变电站真实出力

    % ---- (a) LinDistFlow 解析电压 ----
    p_inj = zeros(nb,1);
    p_inj(1) = g_ac;
    for k = 1:npv
        p_inj(dat.pv_bus(k)) = p_inj(dat.pv_bus(k)) + pv_disp(k);
    end
    p_inj = p_inj - dat.Pd(:,h);
    [~, V_lin] = distflow_eval(L, p_inj, -dat.Qd, 1.0);

    % ---- (b) SOCP 潮流模式：固定同样的注入 ----
    av = pv_disp' ./ dat.pv_cap;
    fix = struct('g', g_ac, 'u', zeros(npv,1), 'ls', zeros(nb,1));
    [~, oS] = opt_dispatch_socp(L, dat, av, 'hour', h, 'fix', fix, ...
                                'vmax', opt.vmax, 'vmin', opt.vmin);
    V_socp = oS.V;

    vmin_ac(i)=min(r.bus(:,VM));  vmin_lin(i)=min(V_lin);  vmin_socp(i)=min(V_socp);
    vmax_ac(i)=max(r.bus(:,VM));  vmax_lin(i)=max(V_lin);  vmax_socp(i)=max(V_socp);
    err_lin(i)  = vmin_lin(i)  - vmin_ac(i);
    err_socp(i) = vmin_socp(i) - vmin_ac(i);
    ncones(i)   = oS.exitflag;

    fprintf('%-8.0f | %9.4f | %9.4f | %9.4f | %+10.5f | %+10.5f\n', ...
        pens(i)*100, vmin_ac(i), vmin_lin(i), vmin_socp(i), err_lin(i), err_socp(i));
end

T = table(pen, vmin_ac, vmin_lin, vmin_socp, err_lin, err_socp, ...
          vmax_ac, vmax_lin, vmax_socp, ncones, ...
    'VariableNames', {'pen','vmin_ac','vmin_lin','vmin_socp','err_lin','err_socp', ...
                      'vmax_ac','vmax_lin','vmax_socp','socp_exitflag'});

%% ---------- 汇总 ----------
k = ~isnan(err_lin) & ~isnan(err_socp);
fprintf('\n--- 电压预测误差（模型 − 真实交流，负号 = 低估）---\n');
fprintf('  LinDistFlow : 平均 %+.5f p.u.  最大 %+.5f p.u.  (RMS %.5f)\n', ...
    mean(err_lin(k)), max(abs(err_lin(k))), sqrt(mean(err_lin(k).^2)));
fprintf('  SOCP        : 平均 %+.5f p.u.  最大 %+.5f p.u.  (RMS %.5f)\n', ...
    mean(err_socp(k)), max(abs(err_socp(k))), sqrt(mean(err_socp(k).^2)));
r = sqrt(mean(err_lin(k).^2)) / max(sqrt(mean(err_socp(k).^2)), eps);
fprintf('\n  -> SOCP 的电压误差是 LinDistFlow 的 1/%.1f\n', r);
fprintf('  -> LinDistFlow 误差方向：%s（平均 %+.5f），即对欠压**乐观**\n', ...
    tern(mean(err_lin(k))>0, '偏高', '偏低'), mean(err_lin(k)));

if opt.save
    outdir = fullfile(root, 'results', 'data');
    if ~exist(outdir, 'dir'), mkdir(outdir); end
    out = fullfile(outdir, 'distflow_model_compare.csv');
    writetable(T, out);
    fprintf('\ndata saved: %s\n', out);
end

if opt.plot
    f = figure('Position', [80 80 1360 440], 'Color', 'w', 'Visible', 'off');

    subplot(1,3,1);
    plot(pen*100, err_lin, 'o-', 'LineWidth', 2, 'Color', [0.85 0.33 0.10]); hold on;
    plot(pen*100, err_socp, 's-', 'LineWidth', 2, 'Color', [0.12 0.47 0.71]);
    yline(0, '--', 'Color', [0.5 0.5 0.5], 'FontSize', 8);
    xlabel('PV penetration [%]'); ylabel('Vmin(model) - Vmin(AC)  [p.u.]');
    title('(a) Voltage prediction error', 'FontSize', 11);
    legend({'LinDistFlow','SOCP'}, 'Location', 'best', 'FontSize', 8);
    grid on; box on;

    subplot(1,3,2);
    plot(pen*100, vmin_ac, 'k-', 'LineWidth', 2); hold on;
    plot(pen*100, vmin_lin, 'o--', 'LineWidth', 1.6, 'Color', [0.85 0.33 0.10]);
    plot(pen*100, vmin_socp, 's--', 'LineWidth', 1.6, 'Color', [0.12 0.47 0.71]);
    yline(opt.vmin, ':', 'Vmin limit', 'Color', [0.5 0.5 0.5], 'FontSize', 8);
    xlabel('PV penetration [%]'); ylabel('Vmin  [p.u.]');
    title('(b) Absolute voltages', 'FontSize', 11);
    legend({'AC (truth)','LinDistFlow','SOCP'}, 'Location', 'best', 'FontSize', 8);
    grid on; box on;

    subplot(1,3,3);
    b = bar([sqrt(mean(err_lin(k).^2)), sqrt(mean(err_socp(k).^2))], 0.5, ...
        'FaceColor', 'flat');
    b.CData(1,:) = [0.85 0.33 0.10];  b.CData(2,:) = [0.12 0.47 0.71];
    set(gca, 'XTickLabel', {'LinDistFlow','SOCP'}, 'FontSize', 9);
    ylabel('RMS of Vmin error  [p.u.]');
    title('(c) Accuracy summary', 'FontSize', 11);
    grid on; box on;
    text(1, sqrt(mean(err_lin(k).^2)), sprintf('  %.5f', sqrt(mean(err_lin(k).^2))), ...
        'FontSize', 8, 'VerticalAlignment', 'bottom');
    text(2, sqrt(mean(err_socp(k).^2)), sprintf('  %.5f', sqrt(mean(err_socp(k).^2))), ...
        'FontSize', 8, 'VerticalAlignment', 'bottom');

    sgtitle('LinDistFlow  vs  SOCP DistFlow  —  same injections, compared against AC', ...
        'FontSize', 12);
    outdir = fullfile(root, 'results', 'figures');
    if ~exist(outdir, 'dir'), mkdir(outdir); end
    outfig = fullfile(outdir, 'distflow_model_compare.png');
    exportgraphics(f, outfig, 'Resolution', 140);
    close(f);
    fprintf('figure saved: %s\n', outfig);
end
end


function s = tern(c, a, b)
if c, s = a; else, s = b; end
end
