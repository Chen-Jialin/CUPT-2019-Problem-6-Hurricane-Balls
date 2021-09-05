close all;clear;clc;
d2 = [0.09907
    0.107
    0.12
    0.133
    0.14
    0.148];
d2_ = [0.05734
    0.05915
    0.06131
    0.04981
    0.05779
    0.05113];
d = 100 * d2 / 2;
% d = 100 * max(d2_,d2 - d2_);
phi = [184.0776945
    121.7397062
    99.87451272
    83.50450939
    79.2102153
    63.40829209];
p1 = -1.018;
p2 = 75.27;
q1 = -0.5424;
v = (p1 .* d + p2) ./ (d + q1);
plot(v,phi / 2 / pi,'.','MarkerSize',20);
hold on;
syms phi1;
g = 9.7964;
i = 1;
r = 0.0095;
m = 28.15 * 10^(-3);
mu = 6.5635e-5;
theta = acos(2 / 5 - g / r / phi1^2);
for v = 9.6:0.01:16.5
%     Fp = -4.746e-6 * (v - r * phi1 * sin(theta))^2 + 0.0003986 * (v - r * phi1 * sin(theta));
    Fr = -4.746e-6 * (r * phi1 * sin(theta))^2 + 0.0003986 * (r * phi1 * sin(theta));
%     equ = 8 * Fp * r * sin(theta) - 2 * Fr * 2 * pi * r * sin(theta) - 2 * mu * m * g * 2 * pi;
    equ = 2 * (4 * (-4.746e-6 * v^2 + 0.0003986 * v) + 0 * (4.746e-6  * v * r * phi1 * sin(theta) - 0.0003986 * r * phi1 * sin(theta)) + 8 / 3 * (-4.746e-6 * r^2 * phi1^2)) * r * sin(theta) - 2 * Fr * 2 * pi * r * sin(theta) - 2 * mu * m * g * 2 * pi;
    phi1val(i) = vpasolve(equ,phi1);
    i = i + 1;
end
plot((9.6:0.01:16.5),phi1val / 2 / pi,'LineWidth',4)
grid on;
ylim([0,100]);
xlabel('风速v/(m\cdot s^{-1})','interpreter','tex','fontsize',18);
ylabel('飓风球转速/(转\cdot s^{-1})','interpreter','tex','fontsize',18);
legend('数据点','拟合线','fontsize',14);