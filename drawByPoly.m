%����������ʽ�õ�S11��S21��Ⱥʱ��
function [] = drawByPoly(E, F, P, ep, ep_r, num)
w1 = -3;
w2 = 3;
dw=0.01;
w=w1:dw:w2;

%ת����wƽ��
root_E = roots(E) ./ 1i;
root_F = roots(F) ./ 1i;
root_P = roots(P) ./ 1i;
Ew=poly(root_E);   
Fw=poly(root_F);
Pw=poly(root_P);

detaPw=polyder(Pw);
detaEw=polyder(Ew);
S11p = zeros(1, length(w));
S21p = zeros(1, length(w));
Tgp = zeros(1, length(w));
for k=1:1:length(w)
    S11p(k)=polyval(Fw,w(k)) / polyval(Ew,w(k)) / ep_r;
    S21p(k)=polyval(Pw,w(k)) / polyval(Ew,w(k)) / ep;
    Tgp(k) = polyval(detaPw,w(k))/polyval(Pw,w(k)) - polyval(detaEw,w(k))/polyval(Ew,w(k));
    Tgp(k)=-imag(Tgp(k));
end
dBS11p=20*log10(S11p);
dBS21p=20*log10(S21p);

figure(num)
subplot(2,1,1);
plot(w, real(dBS11p), 'r', w, real(dBS21p), 'b', 'linewidth', 2)
axis([-inf,inf,-160, 20]);
set(gca,'YTick',-160 : 20 : 10);

grid on
legend('S11','S21',1); 
xlabel('��һ��Ƶ��(Hz)'); 
ylabel('˥��(dB)') ;
title('��������ʽ-S��������');

subplot(2,1,2);
plot(w, Tgp, 'g', 'linewidth', 2);
axis([-inf,inf,0,max(Tgp)+0.5]);
set(gca, 'YTick', 0 :10 : max(Tgp) + 2);
grid on
legend('Ⱥʱ��',1);
xlabel('��һ��Ƶ��(Hz)');
ylabel('Ⱥʱ��(s)');
title('��������ʽ-Ⱥ��ʱ����');