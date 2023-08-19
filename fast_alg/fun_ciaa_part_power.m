function [P_tmp, s_hat]=fun_ciaa_part_power(Y_At_tmp, At, Iter_no,P_tmp)

Mt=size(At,1);
for jj        = 1:Iter_no
%     display(['Iteration ',num2str(jj),'...']);
    R         = At.*P_tmp.'*At' + 1e-2*eye(Mt);
    R_inv     = R \ eye(Mt);
    R_invA    = R_inv*At;
    part1     = sum(conj(R_invA).*Y_At_tmp.').';
    part2     = sum(conj(R_invA).*At, 1).';
    s_hat     = part1./part2;
    P_tmp       = sum(s_hat.*conj(s_hat),2);
end