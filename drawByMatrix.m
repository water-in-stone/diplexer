%由耦合矩阵绘制S参数曲线和群时延曲线
function []  = drawByMatrix(M, N, num)
w1=-3;
w2=3;
dw=0.01;
w=w1:dw:w2;
R1=1;
RN=1;
S11=zeros(1,length(w));
S21=zeros(1,length(w));
Tg=zeros(1,length(w));
for k=1:1:length(w)
    R=zeros(N+2);
    R(1,1)=R1;
    R(N+2,N+2)=RN;
    I=eye(N+2);
    I(1,1)=0;
    I(N+2,N+2)=0;
    Z=R+1j*M+1j*w(k)*I;
    Y=inv(Z);
    A=Y/(-1j);
    S11(k)=1-2*R1*Y(1,1);
    S21(k)=2*sqrt(R1*RN)*Y(N+2,1);
    for kk=2:1:N+1
        Tg(k)=Tg(k)+(A(N+2,kk)*A(kk,1)/A(N+2,1));
    end
    Tg(k)=imag(Tg(k));
end
dBS11=20*log10(S11);
dBS21=20*log10(S21);

figure(num);
subplot(2, 1, 1);
plot(w,real(dBS11),'r',w,real(dBS21),'b','linewidth',2);
axis([-inf, inf, -160, 20]);
set(gca,'YTick',-160 : 20 : 10);
grid on
legend('S11','S21',1); 
xlabel('归一化频率(Hz)'); 
ylabel('衰减(dB)') ;
title('归一化频率-S参数曲线（由耦合矩阵求S参数）');

subplot(2, 1, 2);
plot(w, Tg, 'g', 'linewidth', 2);
axis([-inf,inf,0,max(Tg)+0.5]);
grid on
legend('群时延',1);
xlabel('归一化频率(Hz)');
ylabel('群时延(s)')
title('由耦合矩阵求S参数后的归一化频率下群时延曲线');