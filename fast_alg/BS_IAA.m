function [p] = BS_IAA(x,  N_grid, w, iter_num)
%BS_IAA 此处显示有关此函数的摘要
%   此处显示详细说明

assert(length(w) == N_grid, 'Length is not equal');

% Nt = length(x) / Nr;
% X = reshape(x, [Nr, Nt]);

N_ROI = 0.2 * length(w);
N_sector = 20;
N_samp_sect = N_grid / N_sector;
N_act_sector = 5;
% N = Nr * Nt;
N = length(x);
N_bs = 7;



w = linspace(-pi, pi, N_grid+1);
w(end) = [];
A = exp(-1j * (0: N -1)' * w);


% Inilization

p = abs(A' * x).^2 / N^2;

% train_num = 30;
train_num = 300;
guard_num = train_num/2;



[peak_idx, p_n] = AV_CFAR(p, train_num, guard_num, 1e-5);


ROI_num = length(peak_idx);
ROI_mask_index = zeros(ROI_num, N_ROI);
act_sector_index = zeros(ROI_num, N_act_sector);
for idx = 1: length(peak_idx)
    ROI_mask_index(idx, :) = peak_idx(idx)-N_ROI/2: peak_idx(idx)+N_ROI/2-1;
    tmp = ceil((ROI_mask_index(idx, 1)-0.5)/N_samp_sect): ceil((ROI_mask_index(idx, end)-0.1)/N_samp_sect);
    if length(tmp) == 4
        tmp = [tmp, tmp(end)+1];
    end
    act_sector_index(idx, :) = tmp;
end


%% Bulter matrices Initialization
F = dftmtx(N)/sqrt(N);
B_base = F(:, 1: N_bs);

% keep align with axis
order = (0: N_sector-1)'-N_sector/2;

offset = 2 * pi / N_sector * order;

%% IAA iteration
for idx = 1: ROI_num
    for iter = 1: iter_num
        
%         A_mask = A(:, ROI_mask_index(idx, :));
%         p_mask = p(ROI_mask_index(idx, :));
%         R = A_mask * diag(p_mask) * A_mask' + p_n(idx) * eye(N);
        R = A .* p.' * A' + p_n(idx) * eye(N);
        for k = act_sector_index(idx,:)
            
%             B = diag(exp(-1j*(0:N-1)'*offset(k))) * B_base;
            B = exp(-1j*(0:N-1)'*offset(k)) .* B_base;
            R_iaa_bs = (B' * R * B);
            R_i = inv(R_iaa_bs);
            
            roi_sector_index = (k-1)*N_samp_sect+1: k*N_samp_sect;

            A_sub = B' * A;
            y = B' * x;
            
            for j = roi_sector_index
                s(j) = A_sub(:, j)' * R_i * y  / (A_sub(:, j)' * R_i * A_sub(:, j));
                p(j) =  abs(s(j)).^2;

            end
        end

    end
end


end

