function [p_hat, beta] = QN_PCG_IAA(x, M, K, Maxiter)
%QN-PCG 此处显示有关此函数的摘要
%   此处显示详细说明
x = x(:);
N = length(x);

x0 = fft(x,M)/N;
p_hat=abs(x0).^2;
a_n_old = zeros([N-1, 1]);
for iter = 1: Maxiter
    R = com_R(p_hat, N, M);
    r_n = R(2:end, 1);
    r_sub = R(1: end-1, 1);
    r = R(2: K, 1);
    ww_k = levinson(R(1:K, 1)).';
    a = ww_k(2:end);
    ws_k = [0; conj(flip(a))];
    alpha = R(1) + r' * a;
    
    
    % construct Q^-1
    a_mb = [1; a] / sqrt(alpha);
    a_nb = [a_mb; zeros([N-K-1, 1])];

    
    %%QN-PCG
    
%     a_n = a_n_old;
    a_n = zeros([N-1, 1]);
    e_n = -r_n - comT(r_sub, conj(r_sub), a_n);
%     z = Qi * e_n;
    z = [zeros(N-1-K,1); comRix(ww_k, ws_k, alpha, e_n(N-K: end))]+comRixA(a_nb, e_n, N-K);
    
    e_old = e_n;
    z_old = z;
    p = z;
    
    w = comT(r_sub, conj(r_sub), p);
    alph = e_n'*z / (p'*w);
    a_n = a_n + alph * p;
    e_n = e_n - alph * w;
    pho = e_n'*e_n;
    k = 1;
    while (pho > 1e-4)

%         z = Qi * e_n;
        z = [zeros(N-K-1, 1); comRix(ww_k, ws_k, alpha, e_n(N-K: end))] + comRixA(a_nb, e_n, N-K);
        beta = e_n' * z / (e_old' * z_old);
        p = z + beta * p;
        w = comT(r_sub, conj(r_sub), p);
        alph = e_n'*z / (p' * w);
        a_n = a_n + alph * p;
        
        e_old = e_n;
        e_n = e_n - alph * w;
        
        z_old = z;
        pho = e_n'*e_n;
        k = k + 1;
    end    

    
    ww = [1; a_n];
    ws = [0; conj(flip(a_n))];
    alpha = R(1) + r_n' * a_n;
    
    tt=flip(ww).*(1:N)';
    ts=flip(ws).*(1:N)';
    

    tmp1=comT(tt,[tt(1);zeros(length(ww)-1,1)],conj(ww));
    tmp2=comT(ts,[ts(1);zeros(length(ww)-1,1)],conj(ws));
    

    cb=(tmp1-tmp2)/alpha;
    
    c=[flip(conj(cb)); zeros(M-2*N+1,1); cb(1:end-1)];
    cf=fft(c,M);
    b=[cf(1);flip(cf(2:end))];

    w=p_hat.*real(b).^2;

    z3=comRix(ww,ws,alpha,x);  % Ri*ye
    
    beta=p_hat.*fft(z3,M);
    %update p_hat
    p_hat=abs(beta)./sqrt(w);

        
    a_n_old = a_n;
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