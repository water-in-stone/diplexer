%������Ķ���ʽA(s)���õ�A(-s)
function A = getNegative(A)
N = length(A);
isEven = mod(N, 2) == 0;   %�ж϶���ʽ����ߴ��ݵ���ż��
for i = 1 : 1 : N
    if isEven                               %����ߴ���Ϊ��������������ȡ��
        if mod(i, 2) ~= 0
            A(i) = -A(i);
        end
    else                                        %����ߴ���Ϊ��������ż����ȡ��
        if mod(i ,2) == 0
           A(i) = -A(i);
        end
    end
end
