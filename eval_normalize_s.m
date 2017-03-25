%计算复平面上的归一化的频率
function s = eval_normalize_s(CF,  BW,  ftz)
s = zeros(1, length(ftz));
for k = 1 : 1 : length(ftz)
    s(k) = CF / BW * (ftz(k) / CF - CF / ftz(k));
end
s = s * 1i;   % 转化到复平面


