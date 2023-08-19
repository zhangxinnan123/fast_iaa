%%
% GS-IAA test
%%

clc;
clear all;
close all;

rng(6);
% N = 200;
% L = 150;
% d = randperm(N-1);
% d = sort(d(1:L-1));
% d = [d,N];
L = 16;
% d = [0, 1, 2, 4, 5, 8, 10,11, 13, 15,17, 19]+1;
% d = [1, 2, 4, 7, 14, 21, 28, 32, 36, 37];
% d = [0 1 4 10 16 22 28 30 33 35] + 1;
d = [0 9 12 20 26 35 38 46 65 74 77 78 85 87 90 98]+1;
angle = [-0.75, 0.75];
% angle = [10, 20];
w = angle / 180 * pi;
d_lambda = 0.5;
SNR = 20;
N = 1;
M = length(w);
% K = 10*(d(end)-d(1));
K = 500;
theta = linspace(0, 2*pi*(1-1/K), K);

%% Noncoherent
S = zeros(M ,N);
p1 = 1;
p2 = 1;

for i = 1: M
    S(i, :) = exp(1j * random('unif', -pi, pi, 1, N));
end

% S(2,:) = S(1,:);
% S(2,:) = sqrt(0.01)*S(2,:);
%% Received siganl
Noise = (randn(L, N) + 1j * randn(L, N));
Noise = sqrt(1/ 10^(SNR / 10)) * Noise;

A = exp(1j * d.' * pi * sind(angle));
Y = A * S + Noise;
t_iaa = 0;
t_fmiaa1 = 0;
t_fmiaa2 = 0;
t_das = 0;
%% IAA
N_monto = 20;
for i = 1: N_monto
    A_i = exp(1j * d' * theta);
    t1 = clock;
    [p_iaa, s_iaa] = fun_iaa_power(Y, A_i, 10);
    t2 = clock;
    [p_f, s_f] = miaa_fast2(Y, K, d, 10);
    t3 = clock;
    p_das = abs(A_i' * Y).^2 / L^2;
    t4 = clock;
    [p_f2, s_f2] = miaa_fast(Y, K, d, 10);
    t5 = clock;
    t_iaa = t_iaa + etime(t2, t1);
    t_fmiaa1 = t_fmiaa1 + etime(t3, t2);
    t_das = t_das + etime(t4, t3);
    t_fmiaa2 = t_fmiaa2 + etime(t5, t4);
end
p_iaa = circshift(p_iaa, K/2);
p_f = circshift(p_f, K/2);
p_das = circshift(p_das, K/2);
p_f2 = circshift(p_f2, K/2);

t_iaa/N_monto
t_fmiaa1/N_monto
t_das/N_monto
t_fmiaa2/N_monto

w_grid = linspace(-pi, pi, K+1);
w_grid(end) = [];
theta_iaa = asind(w_grid/pi);

%% plot
figure;
hold on

plot(theta_iaa, 10*log10(p_iaa), 'r', 'LineWidth',2);
plot(theta_iaa, 10*log10(p_f), '--b', 'LineWidth',2);
% plot(theta_iaa, 10*log10(p_f2), '--','Color','#77AC30', 'LineWidth',1.5);
plot(theta_iaa, 10*log10(abs(p_das)), 'k', 'LineWidth',2);
for i = 1: length(w)
    plot([angle(i), angle(i)], [-30, 10], 'b--');
end
xlabel(['\theta']);
ylabel(['Power(dB)'])
axis([-90, 90, -30, 10])
set(gca,'FontSize',18,'FontWeight','bold');
legend('IAA','FMIAA-1', 'DAS')
%% plot2
figure;
hold on

plot(theta_iaa, 10*log10(p_iaa), 'r', 'LineWidth',2);
% plot(theta_iaa, 10*log10(p_f), '--b', 'LineWidth',1.5);
plot(theta_iaa, 10*log10(p_f2), '--','Color','#77AC30', 'LineWidth',2);
plot(theta_iaa, 10*log10(abs(p_das)), 'k', 'LineWidth',2);
for i = 1: length(w)
    plot([angle(i), angle(i)], [-30, 10], 'b--');
end
xlabel(['\theta']);
ylabel(['Power(dB)'])
axis([-90, 90, -30, 10])
set(gca,'FontSize',18,'FontWeight','bold');
legend('IAA','FMIAA-2', 'DAS')
    

    

