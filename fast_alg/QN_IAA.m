function [p_hat, beta] = QN_IAA(x, M, K, Maxiter)
%QN-IAA 此处显示有关此函数的摘要
%   此处显示详细说明
x = x(:);
N = length(x);

x0 = fft(x,M)/N;
p_hat=abs(x0).^2;
for iter = 1: Maxiter
    R = com_R(p_hat, N, M);
    r_n = R(2:end, 1);
    r = R(2: K, 1);
    ww_k = levinson(R(1: K, 1)).';
    a = ww_k(2:end);
    ws_k = [0; conj(flip(a))];
    alpha = R(1) + r' * a;
    
    % construct Q^-1
    a_mb = [1; a] / sqrt(alpha);
    a_nb = [a_mb; zeros([N-K, 1])];

    tt = flip(ww_k).*(1: K)';
    ts = flip(ws_k).*(1: K)';
    
    tmp1=comT(tt, [tt(1);zeros(length(ww_k)-1, 1)], conj(ww_k));
    tmp2=comT(ts, [ts(1);zeros(length(ww_k)-1, 1)], conj(ws_k));
    

    cb=(tmp1-tmp2)/alpha;
    
    c=[flip(conj(cb)); zeros(M-2*K+1,1); cb(1:end-1)];
    cf=fft(c,M);
    b=[cf(1);flip(cf(2:end))];
    
    phi = b + (N-K)*abs(fft(a_mb, M)).^2;
    
    
    aax = comRixA(a_nb, x, N-K+1);
    q_x = [zeros(N-K, 1); comRix(ww_k, ws_k, alpha, x(N-K+1: end))] + aax;

    alpha = fft(q_x, M) ./ phi;
    
    p_hat = abs(alpha).^2;
end
    
    
    
    
end



function R=com_R(p,N,K)
% Compute R using FFT
% R=A*diag(p)*A';
% return the first column of R
x=fft(p,K);
r=x(1:N);
% r=conj(r);
R=toeplitz(r',r); 
% r1(1)=real(r1(1));
end

function y=comT(c,r,x)
% y=A*x  A is Toeplitz  
% A=toeplitz(c,r); c as A's first column and r as A's row
N=length(r);
z=[r; 0; flip(c(2:N))];
Lam=fft(z);
tmp=ifft([x; zeros(N,1)]);
ytmp=fft(Lam.*tmp);
y=ytmp(1:N);
y=y(:);
end

function y=comRix(ww,ws,alpha,x)
% y=inv(R)*x
% inv(R)=(Lw*Lw'-Ls*Ls')/alpha

% Lw=toeplitz(ww,[ww(1);zeros(length(ww)-1,1)]); 
% Ls=toeplitz(ws,[ws(1);zeros(length(ww)-1,1)]);

z1=comT([ww(1);zeros(length(ww)-1,1)],conj(ww),x); 

z2=comT(ww,[ww(1);zeros(length(ww)-1,1)],z1);

z3=comT([ws(1);zeros(length(ww)-1,1)],conj(ws),x);
z4=comT(ws,[ws(1);zeros(length(ww)-1,1)],z3);

y=(z2-z4)./alpha;
end

function y = comRixA(a_n, x, K)
% y = A * A' * x
% A = A = toeplitz(a_nb, zeros([1, N-K+1]));
N = length(a_n);


c1 = zeros(N, 1);
c1(1) = a_n(1);
x2 = comT(c1, conj(a_n), x);
x2 = x2(1: K);

% norm(x2-A'*x)
y = comT(a_n, c1, [x2; zeros(N-K, 1)]);


end

