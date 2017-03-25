%由E(s)、P(s)、F(s)求得耦合矩阵
function [M] = getMatrix(E, F, P, ep, ftz, N, num)
ftz = real(ftz ./ 1i);
nz = length(ftz);

if mod(N - length(ftz), 2) == 0           %依据传输零点的个数考虑是否进行正交化
    P = 1i * P;
end

%计算epcilon_r
if nz==N
    ep_r=ep/sqrt(ep^2-1);
else
    ep_r=1;
end      

EF = E + F/ep_r;                          
for k=N+1:-2:1
    m(k)=real(EF(k));
    n(k)=1i*imag(EF(k));
end
for k=N:-2:1
    m(k)=1i*imag(EF(k));
    n(k)=real(EF(k));
end  
if mod(N, 2) == 0
    y21n = P/ep;
    y22n = n;
    yd=m;
    if nz == N
        Msl=ep*(ep_r-1)/ep_r;
        y21n=y21n-1i*Msl*yd;
    else
        Msl=0;
    end
else
    y21n = P/ep;
    y22n = m;
    yd = n;
    if nz == N
        Msl=ep*(ep_r-1)/ep_r;
        y21n=y21n-1i*Msl*yd;
    else
        Msl=0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%求留并构造耦合矩阵
[r21, lamada]=residue(y21n,yd);
[r22, lamada]=residue(y22n,yd);
Mkk=-imag(lamada);
Mlk=sqrt(real(r22));
Msk=real(r21)./sqrt(real(r22));
M=zeros(N+2);
M=M+diag([0;Mkk;0]);
M(1,:)=[0;Msk;Msl];
M(:,1)=[0;Msk;Msl];
M(N+2,:)=[Msl;Mlk;0];
M(:,N+2)=[Msl;Mlk;0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%转换为箭头型耦合矩阵
for i = 1: 1: N - 1 
    for j = i  + 2: 1: N + 1                                           %消去第(i, j)个元素，对应在M矩阵中的支点为[i, j - 1]，N+2矩阵中的支点为[i + 1, j]
           R = eye(N+2);
           theta = -atan(M(i, j ) / M(i, i + 1)); 
           R(j, j) = cos(theta);  
           R(i + 1, i + 1) = cos(theta);
           R(i + 1, j) = -sin(theta);                               
           R(j, i + 1) = sin(theta);
           M=R * M * R';
    end
end
M=round(M*10000)/10000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%在箭头型耦合矩阵的基础上，进而演变为级联三角形耦合矩阵，
%这里默认传输零点的个数小于 N / 2，这样易于形成三角结构
for i = 1 : 1 : nz
    for j = N : -1 :  2 * i + 1                                                                                     %确定M矩阵中的支点为[j - 1, j]，N+2矩阵中的支点为[j, j + 1]
        R = eye(N + 2);
        if j == N                                                                                                            %初始化theta的值，并旋转得到三角结构的耦合矩阵
            theta = atan(M(N, N + 1) / (ftz(i) + M(N + 1, N + 1)));                     
        else
           theta = atan(M(j, j + 2) / M(j + 1, j + 2));                                     
        end 
           R(j, j) = cos(theta);  
           R(j + 1, j + 1) = cos(theta);
           R(j, j + 1) = -sin(theta);                               
           R(j + 1, j) = sin(theta);
           M=R * M * R';                
    end
end
M=round(M*10000)/10000;
disp(M);
%验证下求得的耦合矩阵
drawByMatrix(M, N, num);


