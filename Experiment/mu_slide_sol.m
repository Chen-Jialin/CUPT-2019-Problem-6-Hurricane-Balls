close all;clear;clc
% 加载小球位移x随时间t变化的数据
load mu_slide_data;

% 做最小二乘法拟合x = v * t + 1/2 * a * t^2
A = [t,t.^2];
va = (A' * A) \ A' * x;
v = va(1);
a = 2 * va(2);
plot(t,x,'.','MarkerSize',12);
hold on;
grid on;
plot(t,v * t + 0.5 * a * t.^2,'LineWidth',1.5);
xlabel('时间 t / s','Interpreter','tex','FontSize',14);
ylabel('小球在斜面上滑动的距离 x / m','Interpreter','tex','FontSize',14)
legend('数据点 (t,x)','拟合曲线 x= v_0t + at^2 / 2, v_0 = 0.03050, a = 1.215','Interpreter','tex','FontSize',10);


% 计算滑动摩擦系数a = g * sin(theta) - mu_1 * g * cos(theta)
m = 32.95 * 10^(-3);
g = 9.7964;
R = 10.001 * 10^(-3);
theta = 23 / 180 * pi;
mu_1 = (g * sin(theta) - a) / (g * cos(theta));