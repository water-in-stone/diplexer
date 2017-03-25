function D = sortImagDesc(D)
n = length(D);
for i = 1 : 1 : n
    for j = (i + 1) : 1 : n
        if imag(D(i)) >= imag(D(j))
            temp = D(i);
            D(i) = D(j);
            D(j) = temp;
        end
    end
end