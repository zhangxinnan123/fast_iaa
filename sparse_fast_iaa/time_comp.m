%%
% GS-IAA test
%%

clc;
clear all;
close all;

% rng(6);
N = 100;
p = 0.05*N;
L_set = p:p:N-p;
len = length(L_set);

angle = [-0.7, 0.7];
w = angle / 180 * pi;

SNR = 20;
snap = 1;
M = length(w);

%% Noncoherent
S = zeros(M ,snap);
p1 = 1;
p2 = 1;

for i = 1: M
    S(i, :) = exp(1j * random('unif', -pi, pi, 1, snap));
end
N_monto = 50;
t_iaa = zeros(len, N_monto);
t_fmiaa1 = zeros(len, N_monto);
t_fmiaa2 = zeros(len, N_monto);

%%

for i = 1: len
    L = L_set(i);
    d = randperm(N-2)+1;
    d = sort(d(2:L-1));
    d = [1,d,N];

    K = 10*(d(end)-d(1));

    theta = linspace(0, 2*pi*(1-1/K), K);


    % Received siganl
    Noise = (randn(L, snap) + 1j * randn(L, snap));
    Noise = sqrt(1/ 10^(SNR / 10)) * Noise;

    A = exp(1j * d.' * pi * sind(angle));
    Y = A * S + Noise;

%% IAA


    for j = 1: N_monto
        A_i = exp(1j * d' * theta);
        t1 = clock;
        [p_iaa, s_iaa] = fun_iaa_power(Y, A_i, 10);
        t2 = clock;
        [p_f, s_f] = miaa_fast2(Y, K, d, 10);
        t3 = clock;
%         p_das = abs(A_i' * Y).^2 / L^2;
%         t4 = clock;
        [p_f2, s_f2] = miaa_fast(Y, K, d, 10);
        t4 = clock;
        
        t_iaa(i, j) = etime(t2, t1);
        t_fmiaa1(i, j) = etime(t3, t2);
        t_fmiaa2(i, j) = etime(t4, t3);
    end

end

% save('N_100_10');
t_iaa = mean(t_iaa, 2);
t_fmiaa1 = mean(t_fmiaa1, 2);
t_fmiaa2 = mean(t_fmiaa2, 2);

%% plot
x_axis = 5:5:95;
figure;
hold on

plot(x_axis, t_iaa./t_iaa, 'b', 'LineWidth',1.5);
plot(x_axis, t_iaa./t_fmiaa1, 'r', 'LineWidth',1.5);
plot(x_axis, t_iaa./t_fmiaa2, 'Color','#77AC30', 'LineWidth',1.5);

xlabel(['阵元数占最大孔径之比']);
ylabel(['运行速度比率'])
set(gca,'FontSize',16);
legend('IAA','FMIAA-1', 'FMIAA-2')

figure(2);
hold on
plot(x_axis, t_iaa, 'b', 'LineWidth',1.5);
plot(x_axis, t_fmiaa1, 'r', 'LineWidth',1.5);
plot(x_axis, t_fmiaa2, 'Color','#77AC30', 'LineWidth',1.5);
xlabel(['阵元数占最大孔径之比']);
ylabel(['运行时间(s)'])
set(gca,'FontSize',16);
legend('IAA','FMIAA-1', 'FMIAA-2')

    
    

    

