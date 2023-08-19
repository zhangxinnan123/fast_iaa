function [p_hat, beta] = miaa_fast(y, K, d, Maxiter)
%MIAA_FAST 此处显示有关此函数的摘要
%   此处显示详细说明
N = d(end);
Ng = length(d);
Nm = N-Ng;
yf = zeros(N, 1);
yf(d) = y;
x0 = fft(yf, K)/N;
p_hat = abs(x0).^2;
s_m = 1:N;
s_m(d) = 0;
s_m(s_m==0) = [];


for iter = 1: Maxiter
    x = fft(p_hat, K);
    r_t = conj(x(1:N));

    r = r_t(2: end);
    wt = levinson(r_t).';
    ws = [0; conj(flip(wt(2: end)))];
    alpha = r_t(1)+r'*wt(2: end);
    
    tt = flip(wt).*(1: N)';
    ts = flip(ws).*(1: N)';
    
    tmp1 = comT(tt, [tt(1); zeros(length(wt)-1, 1)], conj(wt));
    tmp2 = comT(ts, [ts(1); zeros(length(wt)-1, 1)], conj(ws));


    cb = (tmp1-tmp2)/alpha;
    
    c = [flip(conj(cb)); zeros(K-2*N+1,1); cb(1:end-1)];
    cf = fft(c, K);
    phi_d = [cf(1); flip(cf(2: end))];

    z3 = comRix(wt, ws, alpha, yf);  % Ri*ye
    
    phi_n = fft(z3, K);
    
    % step 6, inv(R)
    Lw = toeplitz(wt,[wt(1);zeros(length(wt)-1,1)]); 
    Ls = toeplitz(ws,[ws(1);zeros(length(wt)-1,1)]);
    Lwm = Lw(s_m, :);
    Lsm = Ls(s_m, :);
    srs = (Lwm*Lwm'-Lsm*Lsm')/alpha;
    tmp = (srs+srs')/2;
    
%     tmp2 = (Lw*Lw'-Ls*Ls')/alpha;
%     tmp3 = tmp2(s_m,s_m);
    
    L = inv(chol(tmp));
    SL = zeros(N, Nm);
    SL(s_m, :) = L;

    X = comRix2(wt, ws, alpha, SL);
%     X2 = tmp2(:, s_m);
%     X = X2*L;
%     norm(X-X3)
    phi_n_gamma = fft(X*(L'*z3(s_m)), K);
    db = zeros(N, 1);
    for i = 1: Nm
        db = db + comT(X(:,i),[X(1, i); zeros(N-1,1)], flip(conj(X(:,i))));
    end
    d_gamma = [flip(db); zeros(K-2*N+1, 1); conj(db(1: end-1))];
    df = fft(d_gamma, K);
    phi_d_gamma = [df(1); flip(df(2: end))];
    
    beta = (phi_n-phi_n_gamma)./(phi_d-phi_d_gamma);
    p_hat = abs(beta).^2;
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
[N, M] = size(x);
z=[r; 0; flip(c(2:N))];
Lam=fft(z);
tmp=ifft([x; zeros(N,M)]);
ytmp=fft(Lam.*tmp);
y=ytmp(1:N, :);



