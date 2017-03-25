clear
clf
%����Macchiarella�����ġ�Novel Approach to the Synthesis of Microwave Diplexers���͡�Synthesis of Star-junction Multiplexers��ʵ��
%��һ��˫�������������������Ĺ�������������������ƾͿ���ʵ��
% created by Baikal Young at 2017.01.24

%˫����������Ƶ�ʺʹ���
CF_dip = sqrt(2015 * 2305);                                                         %2155.13MHz                                   
BW_dip = 2305 - 2015;


%����2
RL_rx = 20;
N_rx = 7;
BW_rx = 2125 - 2015;                                                                  %RXƵ�εĴ���
CF_rx = sqrt(2015 * 2125);                                                          %����ƽ��������Ƶ��
ftz_rx = [1815, 2199, 2324.3];
ftz_rx = eval_normalize_s(CF_rx, BW_rx, ftz_rx); 
[E_rx, F_rx, P_rx, ep_rx] = eval_E_F_P(CF_rx, BW_rx, N_rx, RL_rx, ftz_rx, CF_dip, BW_dip, 1);  %�����˫������һ��Ƶ����RXƵ�ε�E(s)��F(s)��P(s)��eplicon
S_rx = (E_rx + F_rx) / 2; 

%����2
RL_tx = 20;
N_tx = 7;
BW_tx = 2305 - 2195;
CF_tx = sqrt(2195 * 2305);
ftz_tx = [2058.8, 2121, 2459.3];
ftz_tx = eval_normalize_s(CF_tx, BW_tx, ftz_tx);                                     %�����һ���Ĵ������
[E_tx, F_tx, P_tx, ep_tx] = eval_E_F_P(CF_tx, BW_tx, N_tx, RL_tx, ftz_tx, CF_dip, BW_dip, 2);  %�����˫������һ��Ƶ����TXƵ�ε�E(s)��F(s)��P(s)��eplicon
S_tx = (E_tx + F_tx) / 2;

%type 2��
n0 = -1;
sc0 = 1.5;
precision = 0.000001;                                                                           %���������μ������õ�S_rx��S_tx֮������
N = poly([roots(F_tx); roots(F_rx); 1.5]);                                           % ��F_tx��F_rx�ĸ�ֱ�ӵõ�N(s)������Ϊ��type 2 �ͣ�1.5 Ϊ����ǻ������ķ������

%��һ�ε�������S_rx��S_tx
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
%���Զ���Ļ�ͼ������֤˫��������������ʽ
drawByPoly2 (N, D, Pt, Pr, n0, p0r, p0t, 3);

%��D��N�Ĺ�ϵ��c0�����ǵ����ȣ�c0����һ����С���鲿����ȡ��ʵ��
c0 = 2 / real(D(2) - N(2));
p0_rx = c0 * p0r;
p0_tx = c0 * p0t;
A = 1 / c0;

%����С���˷�����D_rx��D_tx
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

D_rx = polyfit(root_S_rx, D_rx, N_rx - 1);                  %����С���˷�����ϳ�D_rx��D_tx
D_tx = polyfit(root_S_tx, D_tx, N_tx - 1);           

% ��֤����С���˷���õĽ�����������Ǻϵ�
% f1 = polyval(D_rx2, root_S_rx);
% disp (D_rx - f1);
% f2 = polyval(D_tx2, root_S_tx);
% disp (D_tx - f2);

 
%���Ƴ�TX��RX����������ʽ
F_rx = S_rx - [zeros(1, 1), D_rx];
E_rx = S_rx + [zeros(1, 1), D_rx];
F_tx = S_tx - [zeros(1, 1), D_tx];
E_tx = S_tx + [zeros(1, 1), D_tx];

%��ͼ��֤�µݹ��õ�������Ƶ��
drawByPoly(E_rx, F_rx, P_rx, 1 / p0_rx, 1, 4);
drawByPoly(E_tx, F_tx, P_tx, 1 / p0_tx, 1, 5);

%��ftzΪ��һ��˫����Ƶ���µĴ������
ftz_rx = roots(P_rx);
M_rx = getMatrix(E_rx, F_rx, P_rx, 1 / p0_rx, ftz_rx, N_rx, 6);

ftz_tx = roots(P_tx);
M_tx = getMatrix(E_tx, F_tx, P_tx, 1 / p0_tx, ftz_tx, N_tx, 7);

%����һ����������Ͼ����Լ���������
[f_dig_rx, K_rx, Qe_rx] = antiNormalize(M_rx, c0, CF_dip, BW_dip);
[f_dig_tx, K_tx, Qe_tx] = antiNormalize(M_tx, c0, CF_dip, BW_dip);
FBW = BW_dip / CF_dip;
Qe = c0 / FBW;
disp('ok');




















