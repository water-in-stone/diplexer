%���㸴ƽ���ϵĹ�һ����Ƶ��
function s = eval_normalize_s(CF,  BW,  ftz)
s = zeros(1, length(ftz));
for k = 1 : 1 : length(ftz)
    s(k) = CF / BW * (ftz(k) / CF - CF / ftz(k));
end
s = s * 1i;   % ת������ƽ��


