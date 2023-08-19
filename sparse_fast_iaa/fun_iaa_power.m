function [P, s_hat]=fun_iaa_power(y_noisy, A, Iter_no,s_hat_ini)

[M,t_samples] = size(y_noisy);
if nargin == 4 % initial estimate provided
    s_hat     = s_hat_ini;
else
    s_hat     = A'*y_noisy./sum(conj(A).*A, 1).'; % DAS estimate
end
P             = s_hat.*conj(s_hat)/t_samples;
for jj        = 1:Iter_no
%     display(['Iteration ',num2str(jj),'...']);
    R         = A.*P.'*A';
    clear tmp s_hat;
    R_inv     = inv(R);
    R_invA    = R_inv*A;
%     R_invA    = R\A;
    part1     = R_invA'*y_noisy;
    part2     = sum(conj(R_invA).*A, 1).';
    s_hat     = part1./part2;
%     display(['Elapsed time: part 1:', num2str(t), ' part 2: ', num2str(t1)]);
    P_pre=P;
    P        = sum(s_hat.*conj(s_hat),2)/t_samples;
%     if norm(P-P_pre)^2/norm(P)^2 < 1e-5
%         break
%     end
end
