function test_p1_regression()
%TEST_P1_REGRESSION  P1 验收：case33mg_der 的正确性回归测试
%
%   跑法：matlab -batch "cd('D:\microgrid-dro\src'); test_p1_regression"
%
%   核心闸门是 TEST 1：把 DER 全部关掉后，case33mg_der 必须与标准 case33mg
%   逐字段一致。这保证我们"叠加" DER 时没有偷偷改坏标准算例 —— 一旦改坏，
%   后面所有结果都没法跟文献对标了。
%
%   输出全 ASCII，避免 Windows 控制台 GBK 乱码。

define_constants;
mpopt = mpoption('verbose', 0, 'out.all', 0);
n_pass = 0; n_fail = 0;
tol = 0;      % 回归要求逐比特一致

% case33mg_der 在 ../cases/ 下，脚本从任何地方跑都能找到
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'cases'));

fprintf('\n');
fprintf('############################################################\n');
fprintf('#  P1 regression test : case33mg_der vs case33mg           #\n');
fprintf('############################################################\n\n');

%% ==================== TEST 1: 回归闸门 ====================
fprintf('--- TEST 1: regression gate (all DER off == case33mg) ---\n');
base = case33mg;
bare = case33mg_der('pv_penetration', 0, 'storage', false, 'thermal', false);

checks = { ...
    'bus',     base.bus,     bare.bus; ...
    'branch',  base.branch,  bare.branch; ...
    'gen',     base.gen,     bare.gen; ...
    'gencost', base.gencost, bare.gencost; ...
    'baseMVA', base.baseMVA, bare.baseMVA };
allok = true;
for k = 1:size(checks, 1)
    nm = checks{k, 1}; A = checks{k, 2}; B = checks{k, 3};
    d = max(abs(double(A(:)) - double(B(:))));
    ok = (isequal(size(A), size(B))) && (d <= tol);
    allok = allok && ok;
    fprintf('      %-8s size %-14s  maxdiff = %g   %s\n', nm, ...
        sprintf('%dx%d', size(A,1), size(A,2)), d, tf(ok));
end
[fprintf_ok, n_pass, n_fail] = record(allok, ...
    'bus/branch/gen/gencost 与 case33mg 逐字段一致', n_pass, n_fail);

% 字段级：无 DER 时不应凭空多出 genfuel/isolar/iess
extra = setdiff(fieldnames(bare), fieldnames(base));
fprintf('      extra fields when no DER: %s\n', strjoin(extra(:)', ', '));
[fprintf_ok, n_pass, n_fail] = record(isempty(extra), ...
    '无 DER 时不额外引入 genfuel/isolar/iess 字段', n_pass, n_fail);

%% ==================== TEST 2: 光伏 ====================
fprintf('\n--- TEST 2: PV generators ---\n');
[mpc, der] = case33mg_der();
ng0 = size(base.gen, 1);
fprintf('      gen rows: %d -> %d\n', ng0, size(mpc.gen, 1));
[fprintf_ok, n_pass, n_fail] = record(numel(der.pv_buses) == 4, ...
    '4 个光伏节点已加入', n_pass, n_fail);

pv_buses_out = mpc.gen(ng0 + 1 : ng0 + numel(der.pv_buses), GEN_BUS);
[fprintf_ok, n_pass, n_fail] = record(isequal(pv_buses_out(:)', der.pv_buses), ...
    '光伏接在预期母线 [18 22 25 33]', n_pass, n_fail);

pv_pmax = mpc.gen(ng0 + 1 : ng0 + numel(der.pv_buses), PMAX);
pv_pmin = mpc.gen(ng0 + 1 : ng0 + numel(der.pv_buses), PMIN);
fprintf('      PV PMAX each = %.4f MW (total %.4f), PMIN = %s\n', ...
    pv_pmax(1), sum(pv_pmax), mat2str(unique(pv_pmin)'));
[fprintf_ok, n_pass, n_fail] = record(all(abs(pv_pmin) < 1e-12), ...
    '光伏 PMIN=0（可弃光）', n_pass, n_fail);

pen_check = sum(pv_pmax) / sum(base.bus(:, PD));
fprintf('      realized penetration = %.4f (requested %.4f)\n', pen_check, der.pv_penetration);
[fprintf_ok, n_pass, n_fail] = record(abs(pen_check - der.pv_penetration) < 1e-9, ...
    '渗透率 = 光伏装机 / 总负荷', n_pass, n_fail);

pv_q = mpc.gen(ng0 + 1 : ng0 + numel(der.pv_buses), [QMAX QMIN]);
[fprintf_ok, n_pass, n_fail] = record(all(pv_q(:) == 0), ...
    '光伏 QMAX=QMIN=0（单位功率因数）', n_pass, n_fail);

[fprintf_ok, n_pass, n_fail] = record(isfield(mpc, 'isolar') && ...
    numel(mpc.isolar) == 4, 'mpc.isolar 索引已建（MOST 兼容）', n_pass, n_fail);

%% ==================== TEST 3: 储能 ====================
fprintf('\n--- TEST 3: storage ---\n');
[mpc_s, der_s] = case33mg_der('pv_penetration', 0);
g = size(mpc_s.gen, 1);
fprintf('      storage: bus=%d  PMAX=%.3f  PMIN=%.3f  [Pmin<0 => can charge]\n', ...
    mpc_s.gen(g, GEN_BUS), mpc_s.gen(g, PMAX), mpc_s.gen(g, PMIN));
[fprintf_ok, n_pass, n_fail] = record(mpc_s.gen(g, PMIN) < 0 && mpc_s.gen(g, PMAX) > 0, ...
    '储能 PMIN<0 / PMAX>0（可充可放）', n_pass, n_fail);
[fprintf_ok, n_pass, n_fail] = record(isfield(mpc_s, 'iess') && numel(mpc_s.iess) == 1, ...
    'mpc.iess 索引已建（MOST 兼容）', n_pass, n_fail);
[fprintf_ok, n_pass, n_fail] = record(der_s.storage.ecap > 0 && ...
    der_s.storage.ecap / der_s.storage.pcap >= 1, ...
    '储能能量/功率比 >= 1 小时（容量配置合理）', n_pass, n_fail);
fprintf('      storage energy = %.2f MWh, power = %.2f MW, duration = %.1f h\n', ...
    der_s.storage.ecap, der_s.storage.pcap, der_s.storage.ecap / der_s.storage.pcap);

%% ==================== TEST 4: 热限值 ====================
fprintf('\n--- TEST 4: branch thermal limits ---\n');
[mpc_t, der_t] = case33mg_der('thermal', true, 'i_rated', 300);
fprintf('      I_rated = %d A -> RATE_A = %.4f MVA (uniform)\n', der_t.i_rated, der_t.rate_a_mva);
[fprintf_ok, n_pass, n_fail] = record(all(mpc_t.branch(:, RATE_A) == der_t.rate_a_mva), ...
    '所有支路 RATE_A 已置为统一值', n_pass, n_fail);
fprintf('      expected from physics: sqrt(3)*12.66*300/1000 = %.4f MVA\n', ...
    sqrt(3) * 12.66 * 300 / 1000);
[fprintf_ok, n_pass, n_fail] = record( ...
    abs(der_t.rate_a_mva - sqrt(3) * 12.66 * 300 / 1000) < 1e-9, ...
    'RATE_A 换算公式正确', n_pass, n_fail);

% 基态负载率（原始 case33mg 的潮流，与热限值无关）
r0 = runpf('case33mg', mpopt);
S0 = hypot(r0.branch(:, PF), r0.branch(:, PT));
maxload = max(S0) / der_t.rate_a_mva;
fprintf('      base-case max branch loading = %.1f%% (branch %d)\n', ...
    maxload * 100, find(S0 == max(S0), 1));
[fprintf_ok, n_pass, n_fail] = record(maxload < 1.0, ...
    '基态不越限（RATE_A 选得合理）', n_pass, n_fail);
[~, n_pass, n_fail] = record(maxload > 0.5, ...
    '基态负载率 > 50%（限值没有松到无意义）', n_pass, n_fail);

%% ==================== TEST 5: 潮流与光伏注入生效 ====================
fprintf('\n--- TEST 5: power flow & PV injection ---\n');
mf = case33mg_der('pv_penetration', 0);
r_no = runpf(mf, mpopt);
mf = case33mg_der('pv_penetration', 0.5);
r_pv = runpf(mf, mpopt);
fprintf('      no PV : slack=%.1f kW  loss=%.1f kW  Vmin=%.4f\n', ...
    r_no.gen(1, PG) * 1000, sum(r_no.branch(:, PF) + r_no.branch(:, PT)) * 1000, ...
    min(r_no.bus(:, VM)));
fprintf('      50%% PV: slack=%.1f kW  loss=%.1f kW  Vmin=%.4f\n', ...
    r_pv.gen(1, PG) * 1000, sum(r_pv.branch(:, PF) + r_pv.branch(:, PT)) * 1000, ...
    min(r_pv.bus(:, VM)));
[fprintf_ok, n_pass, n_fail] = record(r_pv.gen(1, PG) < r_no.gen(1, PG), ...
    '加 PV 后 Slack 出力下降（注入生效）', n_pass, n_fail);
[fprintf_ok, n_pass, n_fail] = record( ...
    sum(r_pv.branch(:, PF) + r_pv.branch(:, PT)) < sum(r_no.branch(:, PF) + r_no.branch(:, PT)), ...
    '适度 PV 渗透降低网损', n_pass, n_fail);
[fprintf_ok, n_pass, n_fail] = record(r_pv.success && r_no.success, ...
    '两种配置潮流均收敛', n_pass, n_fail);

% pv_available 应该线性影响注入
mf = case33mg_der('pv_penetration', 0.5, 'pv_available', 0.5);
r_half = runpf(mf, mpopt);
i_pv = size(base.gen, 1) + (1:4);            % 光伏机组在 gen 中的行号
pv_inj_full = sum(r_pv.gen(i_pv, PG));
pv_inj_half = sum(r_half.gen(i_pv, PG));
fprintf('      PV injection: available=1.0 -> %.4f MW ; available=0.5 -> %.4f MW\n', ...
    pv_inj_full, pv_inj_half);
[fprintf_ok, n_pass, n_fail] = record(abs(pv_inj_half - pv_inj_full / 2) < 1e-9, ...
    'pv_available 与注入量成正比', n_pass, n_fail);

%% ==================== TEST 6: 选项校验 ====================
fprintf('\n--- TEST 6: option validation ---\n');
threw = false;
try
    case33mg_der('no_such_option', 1);
catch
    threw = true;
end
[fprintf_ok, n_pass, n_fail] = record(threw, '未知选项报错', n_pass, n_fail);

threw = false;
try
    case33mg_der('pv_penetration');
catch
    threw = true;
end
[fprintf_ok, n_pass, n_fail] = record(threw, '选项数量为奇数时报错', n_pass, n_fail);

%% ==================== TEST 7: 电压限值选项 ====================
fprintf('\n--- TEST 7: voltage limit option ---\n');
md = case33mg_der();
fprintf('      default Vmax = %.2f (case33mg original)\n', md.bus(1, VMAX));
[fprintf_ok, n_pass, n_fail] = record(abs(md.bus(1, VMAX) - 1.1) < 1e-12, ...
    '默认 Vmax 保持 case33mg 的 1.1', n_pass, n_fail);

mt = case33mg_der('vmax', 1.05, 'vmin', 0.95);
fprintf('      with option  Vmax = %.2f / Vmin = %.2f\n', mt.bus(1, VMAX), mt.bus(1, VMIN));
[fprintf_ok, n_pass, n_fail] = record( ...
    all(mt.bus(:, VMAX) == 1.05) && all(mt.bus(:, VMIN) == 0.95), ...
    'vmax / vmin 选项生效', n_pass, n_fail);

% VMAX/VMIN 只是 OPF 的约束边界，潮流不校验 -> runpf 结果不应受影响
pen_hi = 1.8;
m1 = case33mg_der('pv_penetration', pen_hi);
m2 = case33mg_der('pv_penetration', pen_hi, 'vmax', 1.05);
r1 = runpf(m1, mpopt);
r2 = runpf(m2, mpopt);
v1 = max(r1.bus(:, VM));
fprintf('      at %.0f%% PV: Vmax = %.4f  (over 1.05? %d ; over 1.10? %d)\n', ...
    pen_hi * 100, v1, v1 > 1.05, v1 > 1.10);
[fprintf_ok, n_pass, n_fail] = record(v1 > 1.05 && v1 < 1.10, ...
    '180% 渗透率下越 1.05 但未越 1.10（限值选择显著影响结论）', n_pass, n_fail);
[fprintf_ok, n_pass, n_fail] = record(abs(v1 - max(r2.bus(:, VM))) < 1e-9, ...
    'VMAX 是约束边界而非潮流参数（runpf 不受影响）', n_pass, n_fail);

%% ==================== 汇总 ====================
fprintf('\n############################################################\n');
fprintf('#  RESULT: %d passed, %d failed                            \n', n_pass, n_fail);
fprintf('############################################################\n\n');

if n_fail > 0
    error('test_p1_regression: %d 项失败', n_fail);
end
end


%% ---------- helper ----------
function s = tf(ok)
if ok, s = '[PASS]'; else, s = '[FAIL]'; end
end

function [dummy, n_pass, n_fail] = record(ok, name, n_pass, n_fail)
fprintf('      %s %s\n', tf(ok), name);
if ok
    n_pass = n_pass + 1;
else
    n_fail = n_fail + 1;
end
dummy = true;
end
