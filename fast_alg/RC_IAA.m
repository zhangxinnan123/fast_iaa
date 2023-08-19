function [p] = RC_IAA(x, Nr, w, Maxiter)
%RC_IAA 此处显示有关此函数的摘要
%   此处显示详细说明


Nt = length(x) / Nr;

X = reshape(x, [Nr, Nt]);


At = exp(1j * Nr * (0: Nt -1)' * w);
Ar = exp(1j * (0: Nr - 1)' * w);
A = exp(1j * (0: Nt*Nr - 1)' * w);

y = X * conj(At)/Nt;

p = abs(A'*x).^2/(Nt*Nr)^2;


% IAA iteration
for iter = 1: Maxiter

    R = Ar .* p.' * Ar';
    R_i = inv(R);
    R_invA = R_i*Ar;
    
    part1     = sum(conj(R_invA).*y).';
    part2     = sum(conj(R_invA).*Ar, 1).';
    s_hat     = part1./part2;
    p       = sum(s_hat.*conj(s_hat),2);
    
end

end

