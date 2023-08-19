clc;
close all;
load('N_100_30_2')
t_iaa = mean(t_iaa, 2);
t_fmiaa1 = mean(t_fmiaa1, 2);
t_fmiaa2 = mean(t_fmiaa2, 2);

%% plot
x_axis = 5:5:95;
figure;
hold on

plot(x_axis, t_iaa./t_iaa, 'r-o', 'LineWidth',2);
plot(x_axis, t_iaa./t_fmiaa1, 'b-h', 'LineWidth',2);
plot(x_axis, t_iaa./t_fmiaa2, '-+','Color','#77AC30', 'LineWidth',2);

xlabel(['阵元数比率 %']);
ylabel(['运行速度比率'])
set(gca,'FontSize',18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
legend('IAA','FMIAA-1', 'FMIAA-2')

figure(2);
hold on
plot(x_axis, t_iaa, 'r-o', 'LineWidth',2);
plot(x_axis, t_fmiaa1, 'b-h', 'LineWidth',2);
plot(x_axis, t_fmiaa2,'-+', 'Color','#77AC30','LineWidth',2);
xlabel(['阵元数比率 %']);
ylabel(['运行时间(s)'])
set(gca,'FontSize',18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
legend('IAA','FMIAA-1', 'FMIAA-2')