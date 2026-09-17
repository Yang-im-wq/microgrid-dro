function T = test_vroot(varargin)
%TEST_VROOT  变电站电压：固定 1.0 vs 放开为决策变量
%
%   T = test_vroot()
%
%   背景：单时段/多时段模型一直把变电站电压平方固定为 1.0。
%   但真实变电站有有载调压变压器（OLTC），可以调压 —— 交流最优潮流就会把它
%   顶到上限去降损（实测 P3 里交流 OPF 的 Vmax 一直贴在 1.05）。
%   放开它是一个**真实存在的运行自由度**，代价是压缩电压上限的裕度。
%
%   本脚本对若干渗透率比较两种设定的：
%     目标值、变电站电压、越限时段数、以及**交流潮流校验后的真实最低电压**
%   —— 最后一列是关键：简化模型说电压会怎样，得让交流潮流说了算。

opt = struct('pen_list', [0.5 1.0 1.5 2.0 2.5]', 'storage_bus', 4);
for k = 1:2:numel(varargin)
    name = varargin{k};
    if ~ischar(name) || ~isfield(opt, name)
        error('test_vroot: 未知选项 %s', num2str(name));
    end
    opt.(name) = varargin{k + 1};
end

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'src'), fullfile(root, 'cases'));
define_constants;
mpopt = mpoption('verbose', 0, 'out.all', 0);

pens = opt.pen_list(:);  n = numel(pens);
pen   = nan(n,1);
objF  = nan(n,1);  objR  = nan(n,1);      % fixed / released
vrF   = nan(n,1);  vrR   = nan(n,1);      % 变电站电压
vminAC_F = nan(n,1);  vminAC_R = nan(n,1); % 交流校验后的真实最低电压
loadF = nan(n,1);  loadR = nan(n,1);      % 最大支路负载率

fprintf('\n############################################################\n');
fprintf('#  变电站电压：固定 1.0  vs  放开为决策变量                    #\n');
fprintf('############################################################\n');
fprintf('%-8s | %10s %10s | %8s %8s | %10s %10s\n', ...
    'PV[%]', 'obj固定', 'obj放开', 'Vr固定', 'Vr放开', 'AC-Vmin固定', 'AC-Vmin放开');
fprintf('%s\n', repmat('-', 1, 74));

for i = 1:n
    pen(i) = pens(i);
    L   = distflow_lin(case33mg_der('pv_penetration', 0, 'storage', false));
    dat = make_day_data('pv_penetration', pens(i), 'storage_bus', opt.storage_bus);
    av  = dat.pv_avail_max ./ max(dat.pv_cap);

    [xF, oF] = opt_dispatch_lindistflow(L, dat, av, 'vroot', 'fixed');
    [xR, oR] = opt_dispatch_lindistflow(L, dat, av, 'vroot', 'free');
    if oF.exitflag ~= 1 || oR.exitflag ~= 1
        fprintf('%-8.0f | 求解失败 (fixed=%d, free=%d)\n', pens(i)*100, oF.exitflag, oR.exitflag);
        continue
    end
    objF(i) = oF.fval;   objR(i) = oR.fval;
    vrF(i)  = sqrt(xF(oF.idx.iVR));   vrR(i) = sqrt(xR(oR.idx.iVR));

    % 回灌交流潮流，看真实最低电压与最大支路负载
    [vminAC_F(i), loadF(i)] = ac_check(dat, av, xF, oF, pens(i), opt.storage_bus);
    [vminAC_R(i), loadR(i)] = ac_check(dat, av, xR, oR, pens(i), opt.storage_bus);

    fprintf('%-8.0f | %10.1f %10.1f | %8.4f %8.4f | %10.4f %10.4f\n', ...
        pens(i)*100, objF(i), objR(i), vrF(i), vrR(i), vminAC_F(i), vminAC_R(i));
end

T = table(pen, objF, objR, vrF, vrR, vminAC_F, vminAC_R, loadF, loadR, ...
    'VariableNames', {'pen','obj_fixed','obj_released','vroot_fixed','vroot_released', ...
                      'vmin_ac_fixed','vmin_ac_released','load_fixed','load_released'});

k = ~isnan(objF) & ~isnan(objR);
fprintf('\n--- 小结 ---\n');
fprintf('  放开变电站电压后目标值变化: 平均 %+.2f $/d\n', mean(objR(k)-objF(k)));
fprintf('  变电站电压: 固定 %.4f -> 放开后平均 %.4f（范围 %.4f ~ %.4f）\n', ...
    mean(vrF(k)), mean(vrR(k)), min(vrR(k)), max(vrR(k)));
fprintf('  交流校验后的真实最低电压: 固定 %.4f / 放开 %.4f\n', ...
    min(vminAC_F(k)), min(vminAC_R(k)));

outdir = fullfile(root, 'results', 'data');
if ~exist(outdir, 'dir'), mkdir(outdir); end
writetable(T, fullfile(outdir, 'vroot_compare.csv'));
fprintf('\ndata saved: results/data/vroot_compare.csv\n');
end


function [vmin_ac, max_load] = ac_check(dat, av, x, o, pen, storage_bus)
define_constants;
mpopt = mpoption('verbose', 0, 'out.all', 0);
nt = dat.nt; nb = dat.nb; npv = numel(dat.pv_bus);
m0 = case33mg_der('pv_penetration', pen, 'storage', true, ...
                  'storage_bus', storage_bus, 'ramp', 'pmax');
m0.bus(:, QD) = dat.Qd;
vmin_ac = inf;  max_load = 0;
for t = 1:nt
    m = m0;
    m.bus(:, PD) = max(dat.Pd(:,t) - x(o.idx.LS((t-1)*nb + (1:nb))), 0);
    for k2 = 1:npv
        m.gen(m.isolar(k2), PG) = max(av(t,k2)*dat.pv_cap(k2) ...
            - x(o.idx.U((t-1)*npv + k2)), 0);
    end
    m.gen(m.iess, PG) = x(o.idx.D(t)) - x(o.idx.C(t));
    r = runpf(m, mpopt);
    if ~r.success, continue; end
    vmin_ac = min(vmin_ac, min(r.bus(:, VM)));
    max_load = max(max_load, ...
        max(hypot(r.branch(:,PF), r.branch(:,PT))) / m.branch(1, RATE_A) * 100);
end
end
