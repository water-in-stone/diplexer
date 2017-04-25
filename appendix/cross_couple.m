%N+2双终端网络综合,拓扑结构为折叠型的耦合矩阵消元，绘制实频曲线,绘制特征多项式曲线 
%%%%%%%%%%%%%%%%%%%%%%%%%%%定义初始条件
clear
ftz=[-1.6058j  1.5943j];  % RL=22   N=4; 
RL=20;                    %回波损耗
N=5;                      %滤波器阶数
%%%%%%%%%%%%%%%%%%%%%%%%%%%切比雪夫多项式综合
syms w;
nz=length(ftz);
w0=sqrt(w^2-1);
ftz=ftz./1i;              %转化为实频率
U(1)=w-1/ftz(1);             
V(1)=w0*sqrt(1-(1/(ftz(1)^2)));
for k=2:1:N
    if k>nz
        U(k)=w*U(k-1)+w0*V(k-1);
        V(k)=w*V(k-1)+w0*U(k-1);
    else
        U(k)=w*U(k-1)-U(k-1)/ftz(k)+w0*sqrt(1-1/ftz(k)^2)*V(k-1);
        V(k)=w*V(k-1)-V(k-1)/ftz(k)+w0*sqrt(1-1/ftz(k)^2)*U(k-1);
    end
end                           %构造切比雪夫多项式,U为未归一化的F,V为带内反射极大值点
F=sym2poly(U(N));             %将多项式表达式转为系数表示
frz=roots(F);                 %求反射零点
P=poly(ftz);                  %构造归一化P
F=poly(frz);                  %构造归一化F
ep=1/sqrt(10^(RL/10)-1)*abs(polyval(P,1)/polyval(F,1)); %计算epcilon
if nz==N
    ep_r=ep/sqrt(ep^2-1);
else
    ep_r=1;
end                                 %计算epcilon_r

%方法一，利用E(s)的根的分布特点
%转换到S平面，通过求解E(s) * E(s)'（即E(s)的共轭）左半平面的根来得到E(s)
% F = poly(1i * frz);
% P = poly(1i * ftz);
% P=[zeros(1, length(F) - length(P)), P];
% poly1 = 1i *P / ep - F / ep_r;
% poly2 = conj(poly1);
% E_roots = roots(conv(poly1, poly2));
% Eroot = E_roots(real(E_roots) < 0);
% E = poly(Eroot);

 %方法二，交替极点法构造E的根，用于满足赫尔维茨条件
P=[zeros(1,length(F)-length(P)),P];
E=P/ep-1i*F/ep_r;
Eroot=roots(E);
r1=Eroot(imag(Eroot)>0);
r2=conj(Eroot(imag(Eroot)<0));
Eroot=[r1;r2];                     

E=poly(1i*Eroot);       %转化为复频率，并构造多项式系数
F=poly(1i*frz);
P=poly(1i*ftz);

if mod(N-nz,2)==0       %正交化条件
    P=1i*P;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%构造导纳矩阵

EF=E+F/ep_r;
for k=N+1:-2:1
    m(k)=real(EF(k));
    n(k)=1i*imag(EF(k));
end
for k=N:-2:1
    m(k)=1i*imag(EF(k));
    n(k)=real(EF(k));
end  
if mod(N,2)==0
    y21n=P/ep;
    y22n=n;
    yd=m;
    if nz==N
        Msl=ep*(ep_r-1)/ep_r;
        y21n=y21n-1i*Msl*yd;
    else
        Msl=0;
    end
else
    y21n=P/ep;
    y22n=m;
    yd=n;
    if nz==N
        Msl=ep*(ep_r-1)/ep_r;
        y21n=y21n-1i*Msl*yd;
    else
        Msl=0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%求留数以及耦合矩阵构造
[r21 lamada]=residue(y21n,yd);
[r22 lamada]=residue(y22n,yd);
Mkk=-imag(lamada);
Mlk=sqrt(real(r22));
Msk=real(r21)./sqrt(real(r22));
M=zeros(N+2);
M=M+diag([0;-imag(lamada);0]);
M(1,:)=[0;Msk;Msl];
M(:,1)=[0;Msk;Msl];
M(N+2,:)=[Msl;Mlk;0];
M(:,N+2)=[Msl;Mlk;0];
M
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%折叠型拓扑结构旋转消元
   for k=1:1:(N+2)/2                          %行消元
       for j=N+2-k:-1:k+2                     %消去第k,j个元素，支点[j-1,j]
           R=eye(N+2);
           theta=-atan(M(k,j)/M(k,j-1));
           R(j,j)=cos(theta);
           R(j-1,j-1)=cos(theta);
           R(j-1,j)=-sin(theta);
           R(j,j-1)=sin(theta);
           M=R*M*R';
       end
       for i=k+2:1:N+2-(k+1)                          %列消元
           R=eye(N+2);                                %消去第i，N+2-k+1个元素，支点[i,i+1]
           theta=atan(M(i,N+2-k+1)/M(i+1,N+2-k+1));
           R(i,i)=cos(theta);
           R(i+1,i+1)=cos(theta);
           R(i,i+1)=-sin(theta);
           R(i+1,i)=sin(theta);
           M=R*M*R';
       end
   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%由耦合矩阵绘制S参数曲线和群时延曲线
M=round(M*10000)/10000;
w1=-4;
w2=4;
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
    A=Y/(-1i);
    S11(k)=1-2*R1*Y(1,1);
    S21(k)=2*sqrt(R1*RN)*Y(N+2,1);
    for kk=2:1:N+1
        Tg(k)=Tg(k)+(A(N+2,kk)*A(kk,1)/A(N+2,1));
    end
    Tg(k)=imag(Tg(k));
end
dBS11=20*log10(S11);
dBS21=20*log10(S21);

%figure(1);
subplot(3,2,1);
plot(w,real(dBS11),'r',w,real(dBS21),'b','linewidth',2);
axis([-inf,inf,-100,10]);
grid on
legend('S11','S21',1); 
xlabel('归一化频率(Hz)'); 
ylabel('衰减(dB)') ;
title('归一化频率-S参数曲线（由耦合矩阵求S参数）');

%figure(2)
subplot(3,2,2);
plot(w,Tg,'linewidth',2);
axis([-inf,inf,0,max(Tg)+0.5]);
grid on
legend('群时延',1);
xlabel('归一化频率(Hz)');
ylabel('群时延(s)')
title('由耦合矩阵求S参数后的归一化频率下群时延曲线');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%实频率曲线
CF=5550;                     %实频率的中心频率，单位MHz
BW=50;                       %带宽，单位MHz
FBW=BW/CF;
m=M*FBW
syms f0 b0 f w0
f=solve('w0=f0/b0*(f/f0-f0/f)','f');    %由带通的频率变换关系求解实频率与归一化频率的关系
f=f(1);                                 %f有两个解，取正解
fstar=subs(f,[f0 b0 w0],[CF BW w(1)]);   
fstop=subs(f,[f0 b0 w0],[CF BW w(length(w))]);  
f=linspace(fstar,fstop,length(w));  %实际物理频率与归一化频率是一一对应的

%figure(3);
subplot(3,2,3);
plot(f,real(dBS11),'r',f,real(dBS21),'b','linewidth',2);
axis([-inf,inf,-100,10]);
grid on
legend('S11','S21',1); 
xlabel('实频率（MHz）'); 
ylabel('衰减(dB)') ;
title('实频-S参数曲线');

%figure(4)
subplot(3,2,4);
plot(f,Tg,'linewidth',2);
axis([-inf,inf,0,max(Tg)+0.5]);
grid on
legend('群时延',1);
xlabel('实频率(MHz)');
ylabel('群时延(s)');
title('实频-群时延曲线');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%由特征多项式绘制曲线
Ew=poly(Eroot);   
Fw=poly(frz);
Pw=poly(ftz);
detaPw=polyder(Pw);
detaEw=polyder(Ew);
for k=1:1:length(w)
    S11p(k)=polyval(Fw,w(k))/polyval(Ew,w(k))/ep_r;
    S21p(k)=polyval(Pw,w(k))/polyval(Ew,w(k))/ep;
    Tgp(k)=polyval(detaPw,w(k))/polyval(Pw,w(k))-polyval(detaEw,w(k))/polyval(Ew,w(k));
    Tgp(k)=-imag(Tgp(k));
end
dBS11p=20*log10(S11p);
dBS21p=20*log10(S21p);

%figure(5)
subplot(3,2,5);
plot(w,real(dBS11p),'r',w,real(dBS21p),'b','linewidth',2)
axis([-inf,inf,-100,10]);
grid on
legend('S11','S21',1); 
xlabel('归一化频率(Hz)'); 
ylabel('衰减(dB)') ;
title('特征多项式-S参数曲线');
% 
% %figure(6)
% subplot(3,2,6);
% plot(w,Tgp,'linewidth',2)
% axis([-inf,inf,0,max(Tg)+0.5]);
% grid on
% legend('群时延',1);
% xlabel('归一化频率(Hz)');
% ylabel('群时延(s)');
% title('特征多项式-群延时曲线');
