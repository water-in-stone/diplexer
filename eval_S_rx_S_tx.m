% ѭ�����S_tx��S_rx����������
function [S_rx, S_tx, D, delta1, delta2, p0r, p0t] = eval_S_rx_S_tx(S_rx, S_tx, P_tx, P_rx, N, RL_rx, RL_tx, n0, N_rx, N_tx)
N1 = abs(polyval(N, 1i))  ^ 2;
N2 = abs(polyval(N, -1i)) ^ 2;
root_S_tx_old = roots(S_tx);
root_S_rx_old = roots(S_rx);
Pt = conv(P_tx, S_rx);
Pr = conv(P_rx, S_tx);
Pt1 = abs(polyval(Pt, 1i)) ^ 2;
Pt2 = abs(polyval(Pt, -1i)) ^ 2;
Pr1 = abs(polyval(Pr, 1i)) ^ 2;
Pr2 = abs(polyval(Pr, -1i)) ^ 2;

%��TX��RXƵ�εĻز����RL�õ� |P0r| ^2 ��|P0t| ^ 2
A = [Pr1, Pt1; Pr2, Pt2];
B = [N1 * (1 / 10 ^ (-RL_tx / 10) - 1); N2 * (1 / 10 ^ (-RL_rx / 10) - 1)];
C = A \ B;
p0r2 = C(1);
p0r = sqrt(p0r2);
p0t2 = C(2);
p0t = sqrt(p0t2);

result_lossless = getLossless(n0, N, p0r2, Pr, p0t2, Pt);       %�����ڲ�����getLossless 
D_roots = roots(result_lossless);
%D(s)Ϊ�ϸ�ĺն�ά�Ķ���ʽ����ȡʵ����ߵĸ����ɵõ�D(s)
D_root = D_roots(real(D_roots) < 0);                            
D = poly(D_root);

% ����֤����õ�D(s)����D(s) * D*(-s)�Ľ��
% temp = conv(D, conj(getNegative(D)));
% delta = result_lossless - temp;
% disp (delta);

roots_S_rx_S_tx = roots((D - N) / 2);
roots_S_rx_S_tx = sortImagDesc(roots_S_rx_S_tx);
root_S_rx_new =  roots_S_rx_S_tx(1 : N_rx, :);                  % ȡS_rx��S_tx�ĳ˻��ĸ���ǰN_rx������S_rx
S_rx = poly(root_S_rx_new);
roots_flipud = flipud(roots_S_rx_S_tx);                             % ȡS_rx��S_tx�ĳ˻��ĸ��ĺ�N_tx������S_tx
root_S_tx_new = roots_flipud(1 :  N_tx , :);   
S_tx = poly(root_S_tx_new);
%����˴ε��������
delta1 = (root_S_rx_new - root_S_rx_old) ;
delta1 = abs(max(delta1));
delta2 = (root_S_tx_new - root_S_tx_old);
delta2 = abs(max(delta2));



%�������Լ���|D(s)| ^ 2 
function result = getLossless(n0, N, P0r, Pr, P0t, Pt)
%��D*(-s)���
N2 =  conv(N, conj(getNegative(N)));
Pr2 = conv(Pr, conj(getNegative(Pr)));
Pt2 = conv(Pt, conj(getNegative(Pt)));
result = abs(n0) * abs(n0) * N2 + P0r * [zeros(1, length(N2) - length(Pr2)), Pr2] + P0t * [zeros(1, length(N2) - length(Pt2)), Pt2] ;





