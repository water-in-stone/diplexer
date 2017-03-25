%��˫����Ƶ���£������Ӧ����һ��ϵ��
function [f_dig, k , Qe2] = antiNormalize(M, c0, f0, BW)
N = length(M);
FBW = BW / f0;
f_dig = zeros(N, N);
k = zeros(N, N);
for i = 2 : 1 : N - 1
    for j = 2 : 1 : N - 1
        if i == j
             f_dig(i, i) = f0 * (sqrt(1 + (M(i, i) * FBW / 2) ^ 2) - M(i , i) * FBW  / 2);           
        else 
            k(i , j) = FBW * M(i, j);
        end
    end
end
k(1, 2) = FBW * M(1, 2) / sqrt(c0);
%Qe1Ϊ����˵�Qֵ��Qe2Ϊ����˵�Qֵ
Qe2 = 1 / (FBW * M(N - 1, N) ^ 2); 


