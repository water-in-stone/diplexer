%N+2˫�ն������ۺ�,���˽ṹΪ�۵��͵���Ͼ�����Ԫ������ʵƵ����,������������ʽ���� 
%%%%%%%%%%%%%%%%%%%%%%%%%%%�����ʼ����
clear
ftz=[-1.6058j  1.5943j];  % RL=22   N=4; 
RL=20;                    %�ز����
N=5;                      %�˲�������
%%%%%%%%%%%%%%%%%%%%%%%%%%%�б�ѩ�����ʽ�ۺ�
syms w;
nz=length(ftz);
w0=sqrt(w^2-1);
ftz=ftz./1i;              %ת��ΪʵƵ��
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
end                           %�����б�ѩ�����ʽ,UΪδ��һ����F,VΪ���ڷ��伫��ֵ��
F=sym2poly(U(N));             %������ʽ���ʽתΪϵ����ʾ
frz=roots(F);                 %�������
P=poly(ftz);                  %�����һ��P
F=poly(frz);                  %�����һ��F
ep=1/sqrt(10^(RL/10)-1)*abs(polyval(P,1)/polyval(F,1)); %����epcilon
if nz==N
    ep_r=ep/sqrt(ep^2-1);
else
    ep_r=1;
end                                 %����epcilon_r

%����һ������E(s)�ĸ��ķֲ��ص�
%ת����Sƽ�棬ͨ�����E(s) * E(s)'����E(s)�Ĺ�����ƽ��ĸ����õ�E(s)
% F = poly(1i * frz);
% P = poly(1i * ftz);
% P=[zeros(1, length(F) - length(P)), P];
% poly1 = 1i *P / ep - F / ep_r;
% poly2 = conj(poly1);
% E_roots = roots(conv(poly1, poly2));
% Eroot = E_roots(real(E_roots) < 0);
% E = poly(Eroot);

 %�����������漫�㷨����E�ĸ�����������ն�ά������
P=[zeros(1,length(F)-length(P)),P];
E=P/ep-1i*F/ep_r;
Eroot=roots(E);
r1=Eroot(imag(Eroot)>0);
r2=conj(Eroot(imag(Eroot)<0));
Eroot=[r1;r2];                     

E=poly(1i*Eroot);       %ת��Ϊ��Ƶ�ʣ����������ʽϵ��
F=poly(1i*frz);
P=poly(1i*ftz);

if mod(N-nz,2)==0       %����������
    P=1i*P;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���쵼�ɾ���

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�������Լ���Ͼ�����
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�۵������˽ṹ��ת��Ԫ
   for k=1:1:(N+2)/2                          %����Ԫ
       for j=N+2-k:-1:k+2                     %��ȥ��k,j��Ԫ�أ�֧��[j-1,j]
           R=eye(N+2);
           theta=-atan(M(k,j)/M(k,j-1));
           R(j,j)=cos(theta);
           R(j-1,j-1)=cos(theta);
           R(j-1,j)=-sin(theta);
           R(j,j-1)=sin(theta);
           M=R*M*R';
       end
       for i=k+2:1:N+2-(k+1)                          %����Ԫ
           R=eye(N+2);                                %��ȥ��i��N+2-k+1��Ԫ�أ�֧��[i,i+1]
           theta=atan(M(i,N+2-k+1)/M(i+1,N+2-k+1));
           R(i,i)=cos(theta);
           R(i+1,i+1)=cos(theta);
           R(i,i+1)=-sin(theta);
           R(i+1,i)=sin(theta);
           M=R*M*R';
       end
   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����Ͼ������S�������ߺ�Ⱥʱ������
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
xlabel('��һ��Ƶ��(Hz)'); 
ylabel('˥��(dB)') ;
title('��һ��Ƶ��-S�������ߣ�����Ͼ�����S������');

%figure(2)
subplot(3,2,2);
plot(w,Tg,'linewidth',2);
axis([-inf,inf,0,max(Tg)+0.5]);
grid on
legend('Ⱥʱ��',1);
xlabel('��һ��Ƶ��(Hz)');
ylabel('Ⱥʱ��(s)')
title('����Ͼ�����S������Ĺ�һ��Ƶ����Ⱥʱ������');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ʵƵ������
CF=5550;                     %ʵƵ�ʵ�����Ƶ�ʣ���λMHz
BW=50;                       %������λMHz
FBW=BW/CF;
m=M*FBW
syms f0 b0 f w0
f=solve('w0=f0/b0*(f/f0-f0/f)','f');    %�ɴ�ͨ��Ƶ�ʱ任��ϵ���ʵƵ�����һ��Ƶ�ʵĹ�ϵ
f=f(1);                                 %f�������⣬ȡ����
fstar=subs(f,[f0 b0 w0],[CF BW w(1)]);   
fstop=subs(f,[f0 b0 w0],[CF BW w(length(w))]);  
f=linspace(fstar,fstop,length(w));  %ʵ������Ƶ�����һ��Ƶ����һһ��Ӧ��

%figure(3);
subplot(3,2,3);
plot(f,real(dBS11),'r',f,real(dBS21),'b','linewidth',2);
axis([-inf,inf,-100,10]);
grid on
legend('S11','S21',1); 
xlabel('ʵƵ�ʣ�MHz��'); 
ylabel('˥��(dB)') ;
title('ʵƵ-S��������');

%figure(4)
subplot(3,2,4);
plot(f,Tg,'linewidth',2);
axis([-inf,inf,0,max(Tg)+0.5]);
grid on
legend('Ⱥʱ��',1);
xlabel('ʵƵ��(MHz)');
ylabel('Ⱥʱ��(s)');
title('ʵƵ-Ⱥʱ������');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����������ʽ��������
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
xlabel('��һ��Ƶ��(Hz)'); 
ylabel('˥��(dB)') ;
title('��������ʽ-S��������');
% 
% %figure(6)
% subplot(3,2,6);
% plot(w,Tgp,'linewidth',2)
% axis([-inf,inf,0,max(Tg)+0.5]);
% grid on
% legend('Ⱥʱ��',1);
% xlabel('��һ��Ƶ��(Hz)');
% ylabel('Ⱥʱ��(s)');
% title('��������ʽ-Ⱥ��ʱ����');
