% time cost between IAA and GS-IAA
clc;
clear all;
close all;

angle = [-0.25, 0.5];
M = 2;
N = 1;
SNR = 10;
L = 100;


S = exp(1j * random('unif', -pi, pi, M, N));
Noise = (randn(L, N) + 1j * randn(L, N))/sqrt(2)*sqrt(1 / 10^(SNR / 10));

K = 10 * L;

A_s = exp(1j * (0: L-1)' * pi *sind(angle));

Y = A_s * S + Noise;


w_grid = linspace(0, 2*pi*(1-1/K), K);

A = exp(1j * (0: L-1)' * w_grid);

% IAA
p1 = IAA(Y, A, 10);
plot(10*log10(fftshift(p1)));

% GS-IAA
[p2,  beta] = GS_IAA(Y, K, 10);
plot(10*log10(fftshift(p2)));

[p3,  beta1] = QN_PCG_IAA(Y, K, 20, 10);
plot(10*log10(fftshift(p3)));

[p4] = QN_IAA(Y, K, 30, 10);
plot(10*log10(fftshift(p4)));

p5 = RC_IAA(Y, 10, w_grid, 10);
plot(10*log10(fftshift(p5)));

p6 = BS_IAA(Y, K, w_grid, 10);
plot(10*log10(p6));

p7 = C_IAA(Y, 10, 10, 10);
plot(10*log10(p7));

% FSIAA-1
[p8, beta8] = FSIAA_1(Y, K, L/4);
plot(10*log10(fftshift(p8)));
% FSIAA-2
[p9, beta9] = FSIAA_1(Y, K, L/4);
plot(10*log10(fftshift(p9)));