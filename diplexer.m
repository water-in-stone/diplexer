clear
clf
%依据Macchiarella的论文《Novel Approach to the Synthesis of Microwave Diplexers》和《Synthesis of Star-junction Multiplexers》实现
%了一个双工器，对于三工器、四工器的情况参照论文类推就可以实现
% created by Baikal Young at 2017.01.24

%双工器的中心频率和带宽
CF_dip = sqrt(2015 * 2305);                                                         %2155.13MHz                                   
BW_dip = 2305 - 2015;


%例子2
RL_rx = 20;
N_rx = 7;
BW_rx = 2125 - 2015;                                                                  %RX频段的带宽
CF_rx = sqrt(2015 * 2125);                                                          %几何平均求中心频率
ftz_rx = [1815, 2199, 2324.3];
ftz_rx = eval_normalize_s(CF_rx, BW_rx, ftz_rx); 
[E_rx, F_rx, P_rx, ep_rx] = eval_E_F_P(CF_rx, BW_rx, N_rx, RL_rx, ftz_rx, CF_dip, BW_dip, 1);  %计算出双工器归一化频率下RX频段的E(s)、F(s)、P(s)、eplicon
S_rx = (E_rx + F_rx) / 2; 

%例子2
RL_tx = 20;
N_tx = 7;
BW_tx = 2305 - 2195;
CF_tx = sqrt(2195 * 2305);
ftz_tx = [2058.8, 2121, 2459.3];
ftz_tx = eval_normalize_s(CF_tx, BW_tx, ftz_tx);                                     %计算归一化的传输零点
[E_tx, F_tx, P_tx, ep_tx] = eval_E_F_P(CF_tx, BW_tx, N_tx, RL_tx, ftz_tx, CF_dip, BW_dip, 2);  %计算出双工器归一化频率下TX频段的E(s)、F(s)、P(s)、eplicon
S_tx = (E_tx + F_tx) / 2;

%type 2型
n0 = -1;
sc0 = 1.5;
precision = 0.000001;                                                                           %迭代中两次计算所得的S_rx和S_tx之间的误差
N = poly([roots(F_tx); roots(F_rx); 1.5]);                                           % 由F_tx和F_rx的根直接得到N(s)，且因为是type 2 型，1.5 为公共腔所代表的反射零点

%第一次迭代计算S_rx和S_tx
[S_rx, S_tx, D, delta1, delta2, p0r, p0t] = eval_S_rx_S_tx(S_rx, S_tx, P_tx, P_rx, N, RL_rx, RL_tx,n0, N_rx, N_tx);
for i = 2 :10
    disp (i);
    if delta1 <= precision && delta2 <= precision
        disp 'converged!';
        break;
    end
    [S_rx, S_tx, D, delta1, delta2, p0r, p0t] = eval_S_rx_S_tx(S_rx, S_tx, P_tx, P_rx, N, RL_rx, RL_tx, n0, N_rx, N_tx);
end

Pr = conv(P_rx, S_tx);
Pt = conv(P_tx, S_rx);
%用自定义的画图函数验证双工器的特征多项式
drawByPoly2 (N, D, Pt, Pr, n0, p0r, p0t, 3);

%由D和N的关系求c0，考虑到精度，c0会有一个极小的虚部，则取其实部
c0 = 2 / real(D(2) - N(2));
p0_rx = c0 * p0r;
p0_tx = c0 * p0t;
A = 1 / c0;

%用最小二乘法计算D_rx和D_tx
root_S_rx = roots(S_rx);
root_S_tx = roots(S_tx);

D_rx = zeros(length(root_S_rx), 1);
for j = 1 : length(root_S_rx)
    D_rx(j) = polyval(D, root_S_rx(j)) / polyval(S_tx, root_S_rx(j)) /  A;
end

D_tx = zeros(length(root_S_tx), 1);
for i = 1: length(root_S_tx)
    D_tx(i) = polyval(D, root_S_tx(i)) / polyval(S_rx, root_S_tx(i)) /  A;
end

D_rx = polyfit(root_S_rx, D_rx, N_rx - 1);                  %用最小二乘法来拟合出D_rx和D_tx
D_tx = polyfit(root_S_tx, D_tx, N_tx - 1);           

% 验证下最小二乘法求得的结果，几乎是吻合的
% f1 = polyval(D_rx2, root_S_rx);
% disp (D_rx - f1);
% f2 = polyval(D_tx2, root_S_tx);
% disp (D_tx - f2);

 
%反推出TX和RX的特征多项式
F_rx = S_rx - [zeros(1, 1), D_rx];
E_rx = S_rx + [zeros(1, 1), D_rx];
F_tx = S_tx - [zeros(1, 1), D_tx];
E_tx = S_tx + [zeros(1, 1), D_tx];

%画图验证下递归后得到的两个频段
drawByPoly(E_rx, F_rx, P_rx, 1 / p0_rx, 1, 4);
drawByPoly(E_tx, F_tx, P_tx, 1 / p0_tx, 1, 5);

%令ftz为归一化双工器频带下的传输零点
ftz_rx = roots(P_rx);
M_rx = getMatrix(E_rx, F_rx, P_rx, 1 / p0_rx, ftz_rx, N_rx, 6);

ftz_tx = roots(P_tx);
M_tx = getMatrix(E_tx, F_tx, P_tx, 1 / p0_tx, ftz_tx, N_tx, 7);

%反归一化出最后的耦合矩阵以及其他参数
[f_dig_rx, K_rx, Qe_rx] = antiNormalize(M_rx, c0, CF_dip, BW_dip);
[f_dig_tx, K_tx, Qe_tx] = antiNormalize(M_tx, c0, CF_dip, BW_dip);
FBW = BW_dip / CF_dip;
Qe = c0 / FBW;
disp('ok');




















