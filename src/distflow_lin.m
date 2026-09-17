function L = distflow_lin(mpc)
%DISTFLOW_LIN  为辐射状配电网构建线性化 DistFlow (LinDistFlow) 模型
%
%   L = distflow_lin(mpc)
%
%   返回结构体 L，把节点电压表示成节点注入的**线性**函数：
%
%       v_i = v_root - 2 * sum_{k 在 root->i 的路径上} ( r_k * P_k + x_k * Q_k )
%
%   其中 v = V^2（p.u. 的平方），P_k / Q_k 是支路 k 的潮流。
%   对辐射网，支路潮流 = 其下游所有节点注入之和，所以：
%
%       P_br = A_p' * p      ( A_p 是"路径矩阵"：A_p(i,k)=1 当支路 k 在 root->i 路径上 )
%       v    = v_root*1 - 2 * A_p * ( r.*P_br + x.*Q_br )
%
%   于是电压对注入是线性的 —— 这是能把它塞进 LP 的关键。
%
%   ★ 为什么用 v = V^2 而不是直接线性化 V
%     直接写 V_j - V_i = -(rP + xQ) 只在 V ≈ 1 附近才准。本馈线基态末端
%     V = 0.90，偏了 10%，那个近似误差会很大。用 v = V^2 形式（Barber 一脉的
%     标准 LinDistFlow）在同样偏差下精度高一个量级。**误差到底多大不靠猜**
%     —— 见 test_p5_distflow.m，直接用 MATPOWER 全交流潮流对标。
%
%   输出字段：
%     .A_p         (nb x nbr) 路径矩阵（仅含闭合支路）
%     .ci          (nbr x 1)  参与建模的支路在 mpc.branch 中的行号
%     .fb .tb      (nbr x 1)  支路首末母线（外部编号）
%     .r  .x       (nbr x 1)  支路电阻/电抗（p.u.）
%     .nb .nbr     节点数 / 支路数
%     .root        参考母线编号（默认 1）
%     .connected   是否所有母线都连到根（false 说明有孤岛）
%
%   See also test_p5_distflow, case33mg_der.

define_constants;

nb   = size(mpc.bus, 1);
ci   = find(mpc.branch(:, BR_STATUS) == 1);
nbr  = numel(ci);
fb   = mpc.branch(ci, F_BUS);
tb   = mpc.branch(ci, T_BUS);
r    = mpc.branch(ci, BR_R);
x    = mpc.branch(ci, BR_X);

root = 1;                       % case33mg 的平衡母线在 bus 1

% ---- BFS 定根：记录每条母线的上游支路 ----
parent_branch = zeros(nb, 1);
visited = false(nb, 1);
visited(root) = true;
queue = root;
while ~isempty(queue)
    v = queue(1); queue(1) = [];
    for k = 1:nbr
        if fb(k) == v, w = tb(k); else, w = fb(k); end
        if (fb(k) == v || tb(k) == v) && ~visited(w)
            visited(w) = true;
            parent_branch(w) = k;
            queue(end+1) = w;                                       %#ok<AGROW>
        end
    end
end

% ---- 路径矩阵 ----
A_p = zeros(nb, nbr);
for i = 1:nb
    j = i;
    while j ~= root && parent_branch(j) ~= 0
        k = parent_branch(j);
        A_p(i, k) = 1;
        if fb(k) == j, j = tb(k); else, j = fb(k); end
    end
end

L = struct('A_p', A_p, 'ci', ci, 'fb', fb, 'tb', tb, 'r', r, 'x', x, ...
           'nb', nb, 'nbr', nbr, 'root', root, 'connected', all(visited));
end
