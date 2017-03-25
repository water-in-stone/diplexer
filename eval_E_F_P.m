%�б�ѩ�����ʽ�ۺϣ����㵥��Ƶ���µ�E(s)��F(s)��P(s)��������
%ת����˫�������ڵĹ�һ��Ƶ����
function [E, F, P, ep] = eval_E_F_P(CF, BW, N, RL, ftz, CF_dip, BW_dip, num)
syms w;
nz=length(ftz);
w0=sqrt(w^2-1);
ftz=ftz./1i;              %ת��ΪʵƵ��
U(1)=w-1/ftz(1);             
V(1)=w0*sqrt(1-(1/(ftz(1)^2)));
%ʹ�ö���ʽչ���ļ���
for k=2:1:N
    if k>nz
        U(k)=w*U(k-1)+w0*V(k-1);
        V(k)=w*V(k-1)+w0*U(k-1);
    else
        U(k)=w*U(k-1)-U(k-1)/ftz(k)+w0*sqrt(1-1/ftz(k)^2)*V(k-1);
        V(k)=w*V(k-1)-V(k-1)/ftz(k)+w0*sqrt(1-1/ftz(k)^2)*U(k-1);
    end
end                                           %�����б�ѩ�����ʽ,UΪδ��һ����F,VΪ���ڷ��伫��ֵ��
F=sym2poly(U(N));               %������ʽ���ʽתΪϵ����ʾ
frz=roots(F);                           %�������
P=poly(ftz);                            %�����һ��P
F=poly(frz);                            %�����һ��F
ep=1/sqrt(10^(RL/10)-1)*abs(polyval(P,1)/polyval(F,1)); %����epcilon
if nz==N
    ep_r=ep/sqrt(ep^2-1);
else
    ep_r=1;
end                                                                     %����epcilon_r
P=[zeros(1,length(F)-length(P)),P];
E=P/ep-1i*F/ep_r;
Eroot=roots(E);
r1=Eroot(imag(Eroot)>0);
r2=conj(Eroot(imag(Eroot)<0));
Eroot=[r1;r2];                                                  %���漫�㷨����E�ĸ�����������ն�ά������

% �任��������ʽ�ĸ���ʵ��Ƶ���ϣ�Ȼ����ת����˫�������ڵĹ�һ��Ƶ��
syms f0 b0 f w0 
f=solve('w0=f0/b0*(f/f0-f0/f)', 'f');                               %�ɴ�ͨ��Ƶ�ʱ任��ϵ���ʵƵ�����һ��Ƶ�ʵĹ�ϵ
f=f(1);                                                                                 %f�������⣬ȡ����
Eroot_dip = zeros(1, length(Eroot));

%����任��E(w)�ĸ�
for i = 1 : length(Eroot)
    Eroot_dip(i) = subs(f, [f0 b0 w0], [CF BW Eroot(i)]);                                                        %�ɵ�ǰͨ��������Ƶ��CF�ʹ���BWת������Ӧ��ʵƵ��
    Eroot_dip(i) = CF_dip / BW_dip * (Eroot_dip(i) / CF_dip - CF_dip / Eroot_dip(i));    % ��ʵƵ�ϵĵ�ת����˫�����Ĺ�һ��Ƶ����
end

%����任���F(w)�ĸ�
Froot_dip = zeros(1, length(frz));
for i = 1 : length(frz)
    Froot_dip(i) = subs(f, [f0 b0 w0], [CF BW frz(i)]);
    Froot_dip(i) = CF_dip / BW_dip * (Froot_dip(i) / CF_dip - CF_dip / Froot_dip(i));
end

%����任���P(w)�ĸ�
Proot_dip = zeros(1, length(ftz));
for i = 1 : length(ftz)
    Proot_dip(i) = subs(f, [f0 b0 w0], [CF BW ftz(i)]);
    Proot_dip(i) =  CF_dip / BW_dip * (Proot_dip(i) / CF_dip - CF_dip / Proot_dip(i));
end

E = poly(1i * Eroot_dip);
F = poly(1i * Froot_dip);
P = poly(1i * Proot_dip);

% �ٴμ���epcilon������Ҫע�⣬ԭ��ep���ɹ�һ��Ƶ����Ϊs= ��j���Ļز����
% �����ģ���������˱任�󣬶���RXƵ�Σ���w = -1�Ļز���Ĳ�Ϊ22dB
if CF < CF_dip
    ep=1 / sqrt(10 ^ (RL / 10) - 1) * abs(polyval(P, -1j) / polyval(F, -1j)); 
else
    ep=1 / sqrt(10 ^ (RL / 10) - 1) * abs(polyval(P, 1j) / polyval(F, 1j)); 
end

if nz == N
    ep_r=ep/sqrt(ep^2-1);
else
    ep_r=1;
end
%�����Զ���Ļ�ͼ����������֤
drawByPoly(E, F, P, ep, ep_r, num);





