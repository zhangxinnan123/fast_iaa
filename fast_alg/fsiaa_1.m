function [p, beta] = FSIAA_1(x, M, Ns)
%FSIAA_1 此处显示有关此函数的摘要
%   此处显示详细说明
x = x(:);
N = length(x);
L = N / Ns;
X = reshape(x, [Ns, L]);
x0=fft(x,M)/N;
p_hat=abs(x0).^2;
p_hat = repmat(p_hat, [1, L]);
Maxiter = 10;
iter = 1;
R = zeros([Ns, Ns, L]);
while(iter <= Maxiter)
    l = 1;
    while (l <= L)
        R(:, :, l) = com_R(p_hat(:, l), Ns, M);
        
        r = R(2:end, 1, l);
        ww = levinson(R(:, 1, l)).';
        ws=[0; conj(flip(ww(2:end)))];
        alpha=R(1,1,1)+r'*ww(2:end);
        
        tt=flip(ww).*(1:Ns)';
        ts=flip(ws).*(1:Ns)';
        
        tmp1=comT(tt,[tt(1);zeros(length(ww)-1,1)],conj(ww));
        tmp2=comT(ts,[ts(1);zeros(length(ww)-1,1)],conj(ws));


        cb=(tmp1-tmp2)/alpha;

        c=[flip(conj(cb)); zeros(M-2*Ns+1,1); cb(1:end-1)];
        cf=fft(c,M);
        b=[cf(1);flip(cf(2:end))];
    %     b = ifft(c, M) * M;
        w=p_hat(:, l).*real(b).^2;

        z3=comRix(ww,ws,alpha,X(:, l));  % Ri*ye

        beta=p_hat(:, l).*(fft(z3,M).*exp(-1j*2*pi*(l-1)*Ns*(0:M-1)'/M));
        %update p_hat
        p_hat(:, l)=abs(beta)./sqrt(w);
        l = l + 1;
    end
    iter = iter + 1;
end
p = mean(p_hat, 2);

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