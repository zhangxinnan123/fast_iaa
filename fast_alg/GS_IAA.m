% GS-IAA algorithm
function [p_hat, beta]=GS_IAA(x, M, Maxiter)

x=x(:);
N=length(x);
if nargin==1
    M=N*10;
end

x0=fft(x, M)/N;
p_hat=abs(x0).^2;
% Maxiter=10;

for iter = 1: Maxiter
    
    % covraince update;
    R = com_R(p_hat, N, M);  

    %(Update weights)
    r = R(2: end,1);
    ww = levinson(R(:, 1)).';
    ws = [0; conj(flip(ww(2: end)))];
    alpha = R(1)+r'*ww(2: end);
    
    
    tt = flip(ww).*(1: N)';
    ts = flip(ws).*(1: N)';
    

    tmp1 = comT(tt, [tt(1); zeros(length(ww)-1, 1)], conj(ww));
    tmp2 = comT(ts, [ts(1); zeros(length(ww)-1, 1)], conj(ws));
    

    cb=(tmp1-tmp2)/alpha;
    
    c=[flip(conj(cb)); zeros(M-2*N+1,1); cb(1:end-1)];
    cf=fft(c, M);
    b=[cf(1); flip(cf(2: end))];

    w=p_hat.*real(b).^2;

    z3=comRix(ww, ws, alpha, x);  % Ri*ye

    beta=p_hat.*fft(z3, M);
    %update p_hat
    p_hat=abs(beta)./sqrt(w);

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
y=y(:);


function y=comT2(c,r,x)
% y=A*x  A is Toeplitz  
% A=toeplitz(c,r); c as A's first column and r as A's row
N=length(r);
z=[c; 0; flip(r(2:N))];
Lam=fft(z);
tmp=fft([x; zeros(N,1)]);
ytmp=ifft(Lam.*tmp);
y=ytmp(1:N);
y=y(:);


