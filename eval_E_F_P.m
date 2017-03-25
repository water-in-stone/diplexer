%切比雪夫多项式综合，计算单个频道下的E(s)、F(s)和P(s)，并将其
%转换到双工器所在的归一化频带上
function [E, F, P, ep] = eval_E_F_P(CF, BW, N, RL, ftz, CF_dip, BW_dip, num)
syms w;
nz=length(ftz);
w0=sqrt(w^2-1);
ftz=ftz./1i;              %转化为实频率
U(1)=w-1/ftz(1);             
V(1)=w0*sqrt(1-(1/(ftz(1)^2)));
%使用多项式展开的技术
for k=2:1:N
    if k>nz
        U(k)=w*U(k-1)+w0*V(k-1);
        V(k)=w*V(k-1)+w0*U(k-1);
    else
        U(k)=w*U(k-1)-U(k-1)/ftz(k)+w0*sqrt(1-1/ftz(k)^2)*V(k-1);
        V(k)=w*V(k-1)-V(k-1)/ftz(k)+w0*sqrt(1-1/ftz(k)^2)*U(k-1);
    end
end                                           %构造切比雪夫多项式,U为未归一化的F,V为带内反射极大值点
F=sym2poly(U(N));               %将多项式表达式转为系数表示
frz=roots(F);                           %求反射零点
P=poly(ftz);                            %构造归一化P
F=poly(frz);                            %构造归一化F
ep=1/sqrt(10^(RL/10)-1)*abs(polyval(P,1)/polyval(F,1)); %计算epcilon
if nz==N
    ep_r=ep/sqrt(ep^2-1);
else
    ep_r=1;
end                                                                     %计算epcilon_r
P=[zeros(1,length(F)-length(P)),P];
E=P/ep-1i*F/ep_r;
Eroot=roots(E);
r1=Eroot(imag(Eroot)>0);
r2=conj(Eroot(imag(Eroot)<0));
Eroot=[r1;r2];                                                  %交替极点法构造E的根，用于满足赫尔维茨条件

% 变换特征多项式的根到实际频率上，然后再转换到双工器所在的归一化频率
syms f0 b0 f w0 
f=solve('w0=f0/b0*(f/f0-f0/f)', 'f');                               %由带通的频率变换关系求解实频率与归一化频率的关系
f=f(1);                                                                                 %f有两个解，取正解
Eroot_dip = zeros(1, length(Eroot));

%计算变换后E(w)的根
for i = 1 : length(Eroot)
    Eroot_dip(i) = subs(f, [f0 b0 w0], [CF BW Eroot(i)]);                                                        %由当前通道的中心频率CF和带宽BW转换到对应的实频上
    Eroot_dip(i) = CF_dip / BW_dip * (Eroot_dip(i) / CF_dip - CF_dip / Eroot_dip(i));    % 将实频上的点转换到双工器的归一化频带中
end

%计算变换后的F(w)的根
Froot_dip = zeros(1, length(frz));
for i = 1 : length(frz)
    Froot_dip(i) = subs(f, [f0 b0 w0], [CF BW frz(i)]);
    Froot_dip(i) = CF_dip / BW_dip * (Froot_dip(i) / CF_dip - CF_dip / Froot_dip(i));
end

%计算变换后的P(w)的根
Proot_dip = zeros(1, length(ftz));
for i = 1 : length(ftz)
    Proot_dip(i) = subs(f, [f0 b0 w0], [CF BW ftz(i)]);
    Proot_dip(i) =  CF_dip / BW_dip * (Proot_dip(i) / CF_dip - CF_dip / Proot_dip(i));
end

E = poly(1i * Eroot_dip);
F = poly(1i * Froot_dip);
P = poly(1i * Proot_dip);

% 再次计算epcilon，这里要注意，原先ep是由归一化频段下为s= ±j处的回波损耗
% 决定的，则这边做了变换后，对于RX频段，则w = -1的回波损耗才为22dB
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
%调用自定义的画图函数进行验证
drawByPoly(E, F, P, ep, ep_r, num);





