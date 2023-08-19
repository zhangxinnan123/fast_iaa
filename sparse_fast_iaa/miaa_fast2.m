function [p_hat, s_hat] = miaa_fast2(y, K, d, Maxiter)
%MIAA_FAST2 此处显示有关此函数的摘要
%   此处显示详细说明
N = d(end);
Ng = length(d);
Nm = N-Ng;
yf = zeros(N, 1);
yf(d) = y;
x0 = fft(yf, K)/N;
p_hat = abs(x0).^2;
s_g = d;
s_m = 1:N;
s_m(d) = 0;
s_m(s_m==0) = [];
g_index = zeros(Ng, Ng);
g_in = zeros(Ng, Ng);
for i = 1: Ng
    for j = 1: Ng
        if i >= j
            g_index(i, j) = d(i) - d(j)+N+1;
        else
            g_index(i, j) = d(i) - d(j)+N;
        end
    end
end
for i = 1: Ng
    for j = i: Ng
        g_in(i, j) = d(j)-d(i);
    end
end
        
for iter = 1: Maxiter
    x = fft(p_hat, K);
    r = conj(x(1:N));
    r = [flip(conj(r)); r];
    
    Rng = r(g_index);
    inv_R = inv(Rng);
    c = zeros(N, 1);
    
    for i = 1: Ng
        for j = i: Ng
            c(g_in(i,j)+1) = c(g_in(i,j)+1) + inv_R(i, j);
        end
    end
    cb = [c; zeros(K-2*N+1, 1); conj(flip(c(2: end)))];
    cf = fft(cb, K);
    phi_d = [cf(1); flip(cf(2: end))];
    
    yr = inv_R*y;
    yr2 = zeros(N, 1);
    yr2(d) = yr;

    phi_n = fft(yr2, K);

    s_hat = phi_n ./ phi_d;
    p_hat = abs(s_hat).^2;
    
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

function y=comRix2(ww,ws,alpha,x)
% y=inv(R)*x
% inv(R)=(Lw*Lw'-Ls*Ls')/alpha

% Lw=toeplitz(ww,[ww(1);zeros(length(ww)-1,1)]); 
% Ls=toeplitz(ws,[ws(1);zeros(length(ww)-1,1)]);

z1=comT2([ww(1);zeros(length(ww)-1,1)],conj(ww),x); 

z2=comT2(ww,[ww(1);zeros(length(ww)-1,1)],z1);

z3=comT2([ws(1);zeros(length(ww)-1,1)],conj(ws),x);
z4=comT2(ws,[ws(1);zeros(length(ww)-1,1)],z3);

y=(z2-z4)./alpha;

function R=com_R(p,N,K)
% Compute R using FFT
% R=A*diag(p)*A';
% return the first column of R
x=fft(p,K);
r=x(1:N);
% r=conj(r);
R=toeplitz(r',r); 
% r1(1)=real(r1(1));

function y=comT(c,r,x)
% y=A*x  A is Toeplitz  
% A=toeplitz(c,r); c as A's first column and r as A's row
N=length(r);
z=[r; 0; flip(c(2:N))];
Lam=fft(z);
tmp=ifft([x; zeros(N,1)]);
ytmp=fft(Lam.*tmp);
y=ytmp(1:N);


function y=comT2(c,r,x)
% y=A*X  A is Toeplitz  
% A=toeplitz(c,r); c as A's first column and r as A's row
[~, M] = size(x);
N=length(r);
z=[r; 0; flip(c(2:N))];
Lam=fft(z);
tmp=ifft([x; zeros(N,M)]);
ytmp=fft(Lam.*tmp);
y=ytmp(1:N, :);