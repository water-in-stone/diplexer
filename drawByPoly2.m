%由特征多项式得到S11、S21、S31和群时延
function [] = drawByPoly2 (N, D, Pt, Pr, n0, p0r, p0t, num)
w1 = -3;
w2 = 3;
dw = 0.01;
w = w1 : dw : w2;
%将多项式D(s)、N(s)、Pt(s)、Pr(s)转换到w平面
root_N = roots(N) ./ 1i;
root_Pr = roots(Pr) ./ 1i;
root_Pt = roots(Pt) ./ 1i;
root_D = roots(D) ./ 1i;
Nw = poly(root_N);
Dw = poly(root_D);
Ptw = poly(root_Pt);
Prw = poly(root_Pr);

detaPrw=polyder(Prw);
detaPtw=polyder(Ptw);
detaDw=polyder(Dw);
S11p = zeros(1, length(w));
S21p = zeros(1, length(w));
S31p = zeros(1, length(w));
Tgp1 = zeros(1, length(w));
Tgp2 = zeros(1, length(w));
for k=1:1:length(w)
    S11p(k) = n0 * polyval(Nw,w(k)) / polyval(Dw,w(k));
    S21p(k) = p0t *  polyval(Ptw,w(k)) / polyval(Dw,w(k));
    S31p(k) = p0r *  polyval(Prw,w(k)) / polyval(Dw,w(k));
    Tgp1(k) =  polyval(detaPtw, w(k)) / polyval(Ptw, w(k)) - polyval(detaDw, w(k)) / polyval(Dw, w(k));
    Tgp1(k)=-imag(Tgp1(k));
    Tgp2(k) =  polyval(detaPrw, w(k)) / polyval(Prw, w(k)) - polyval(detaDw, w(k)) / polyval(Dw, w(k));
    Tgp2(k)=-imag(Tgp2(k));
end
dBS11p=20*log10(S11p);
dBS21p=20*log10(S21p);
dBS31p=20*log10(S31p);

figure(num);
subplot(2, 1, 1);
plot(w,real(dBS11p),'r', w, real(dBS21p), 'b', w, real(dBS31p), 'g', 'linewidth', 2);
axis([-inf, inf, -180, 20]);
set(gca, 'YTick', -180 : 20 : 10);
grid on
legend('S11','S21','S31',1); 
xlabel('归一化频率(Hz)'); 
ylabel('衰减(dB)') ;
title('双工器的S11、S21、S31参数曲线');

subplot(2, 1, 2);
plot(w, Tgp1, 'g', w, Tgp2, 'b', 'linewidth', 2);
axis([-inf, inf, -10, max(max(Tgp1), max(Tgp2)) + 2]);
grid on
legend('群时延',1);
xlabel('归一化频率(Hz)');
ylabel('群时延(s)');
title('特征多项式-群延时曲线');