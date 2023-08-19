function [peak_idx, p_n] = AV_CFAR(x, train, guard, rate_fa)
%AV_CFAR 此处显示有关此函数的摘要
%   此处显示详细说明
L = length(x);
L_train = round(train/2);
L_guard = round(guard/2);
L_side = L_train + L_guard;

alpha = train * (rate_fa^(-1/train)-1);

peak_idx = [];
p_n = [];
for i = L_side+1: L - L_side

    [~, idx] = max(x(i-L_side: i+L_side));
    if i ~= i - L_side + idx-1
        continue
    end
    p_noise = (sum(x(i-L_side: i+L_side)) - sum(x(i-L_guard: i+L_guard)))/train;
    th = alpha * p_noise;
    if x(i) > th
        peak_idx = [peak_idx; i];
        p_n = [p_noise; p_n];
    end

    
end

end

