%对输入的多项式A(s)，得到A(-s)
function A = getNegative(A)
N = length(A);
isEven = mod(N, 2) == 0;   %判断多项式的最高次幂的奇偶性
for i = 1 : 1 : N
    if isEven                               %若最高次幂为奇数，则奇数项取反
        if mod(i, 2) ~= 0
            A(i) = -A(i);
        end
    else                                        %若最高次幂为奇数，则偶数项取反
        if mod(i ,2) == 0
           A(i) = -A(i);
        end
    end
end
