function [v, V] = distflow_eval(L, p_inj, q_inj, v_root)
%DISTFLOW_EVAL  用 LinDistFlow 由节点注入算出各节点电压
%
%   [v, V] = distflow_eval(L, p_inj, q_inj, v_root)
%
%   L      : distflow_lin 返回的模型
%   p_inj  : (nb x 1) 节点有功净注入 MW（发电 - 负荷）
%   q_inj  : (nb x 1) 节点无功净注入 MVAr
%   v_root : 根节点电压平方 p.u.（省略则 1.0）
%
%   v : 各节点电压平方 p.u.     V : 各节点电压 p.u.
%
%   See also distflow_lin, test_p5_distflow.

if nargin < 4 || isempty(v_root)
    v_root = 1.0;
end

% ★ 符号约定（踩过坑）：DistFlow 里 P_k 是"从根流向叶"为正。
%   而节点注入 p_inj = 发电 - 负荷，对负荷节点是负的。
%   辐射网里 支路潮流 = -(其下游所有节点注入之和)，所以必须取负号。
%   漏掉这个负号会让压降符号翻转 —— 表现为"基态电压不降反升"、
%   最小电压恒等于根节点电压，并且支路潮流与 MATPOWER 对不上。
P_br = -L.A_p' * p_inj(:);
Q_br = -L.A_p' * q_inj(:);
dv   = 2 * L.A_p * (L.r .* P_br + L.x .* Q_br);
v    = v_root - dv;
V    = sqrt(max(v, 0));
end
