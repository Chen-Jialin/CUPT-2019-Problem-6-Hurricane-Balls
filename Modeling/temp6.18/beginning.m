% calculation and show the begining motion (with slipping) of the hurricane balls
close all;clear;clc;
T_max = 10;
[T,Y] = ode23('f1',[0,T_max],[0,0,0,0.33796357228538765,0,0,0,20 * pi,0,-20 * pi]);
% for i = 1:size(Y,1)
%     if Y(i,4) > (pi / 2)
%         Y(i,4) = pi / 2;
%     end
% end
g = 9.8;
mu = 0.25;
m = 32.93 * 10^(-3);
R = 10 * 10^(-3);
x = Y(:,1);
y = Y(:,2);
phi = Y(:,3);
theta = Y(:,4);
psi = Y(:,5);
x1 = Y(:,6);
y1 = Y(:,7);
phi1 = Y(:,8);
theta1 = Y(:,9);
psi1 = Y(:,10);
vx = x1 + R * phi1 .* sin(theta) .* sin(phi) - R * psi1 .* sin(theta) .* sin(phi) - 2 * R * theta1 .* cos(theta / 2).^2 .* cos(phi);
vy = y1 - R * phi1 .* sin(theta) .* cos(phi) + R * psi1 .* sin(theta) .* cos(phi) - 2 * R * theta1 .* cos(theta / 2).^2 .* sin(phi);
xcomponent = vx ./ sqrt(vx.^2 + vy.^2);
ycomponent = vy ./ sqrt(vx.^2 + vy.^2);
N = (2*(7*g*m - 7*R*m*theta1.^2.*cos(theta) - 5*R*m*phi1.^2.*cos(theta).*sin(theta).^2 + 2*R*m.*phi1.*psi1.*sin(theta).^2))./(5*sin(theta).^2 + 10*mu*xcomponent.*cos(theta/2).^2.*cos(phi).*sin(theta) + 10*mu.*ycomponent.*cos(theta/2).^2.*sin(phi).*sin(theta) + 7);

figure(1);
subplot(1,2,1);plot(Y(:,1),Y(:,2));title('质心移动轨迹');xlabel('$x$','Interpreter','latex');ylabel('y','Interpreter','latex');axis equal;
subplot(1,2,2);plot(T,Y(:,6));hold on;
plot(T,Y(:,7));hold off;title('质心速度x,y分量');xlabel('T','Interpreter','latex');ylabel('v','Interpreter','latex');legend('$\dot{x}$','$\dot{y}$','Interpreter','latex');

figure(2);
%     subplot(1,2,1);plot(x - R * sin(theta) .* cos(phi),y - R * sin(theta) .* sin(phi));title('触地点移动轨迹');xlabel('$x$','Interpreter','latex');ylabel('$y$','Interpreter','latex');axis equal;
%     subplot(1,2,2);
plot(T,x1 + R * phi1 .* sin(theta) .* sin(phi) - R * psi1 .* sin(theta) .* sin(phi) - 2 * R * theta1 .* cos(theta / 2).^2 .* cos(phi),'LineWidth',2);hold on;
plot(T,y1 - R * phi1 .* sin(theta) .* cos(phi) + R * psi1 .* sin(theta) .* cos(phi) - 2 * R * theta1 .* cos(theta / 2).^2 .* sin(phi),'LineWidth',2);hold off;
legend('v_x','v_y','Interpreter','tex')
title('触地点速度x,y分量','Interpreter','tex');
xlabel('时间t / s','Interpreter','tex');
ylabel('速度/(m\cdot s^{-1})','Interpreter','tex')
grid on;

figure(3);
plot(T,Y(:,8),'LineWidth',2);hold on;
plot(T,Y(:,10),'LineWidth',2);hold off,xlabel('t','Interpreter','tex'),ylabel('\omega','Interpreter','tex'),legend('$\dot{\phi}$','$\dot{\psi}$','Interpreter','latex')
grid on;
xlabel('时间t / s','Interpreter','tex');
ylabel('角速度/s^{-1}','Interpreter','tex');
legend('公转角速度','自转角速度','Interpreter','tex');

figure(4);
plot(T,Y(:,4),'LineWidth',2);
xlabel('时间t / s','Interpreter','tex'),ylabel('两球心连线与竖直方向夹角\theta','Interpreter','tex');
grid on;

% figure(5);
% fI = mu * N .* xcomponent .* sin(phi) - mu * N .* ycomponent .* cos(phi);
% plot(T,fI,'LineWidth',2);
% xlabel('时间t / s','Interpreter','tex');
% ylabel('侧视图中垂直纸面摩擦力f_{I}','Interpreter','tex');
% grid on;
% 
% figure(6);
% fII = -mu * N .* xcomponent .* cos(phi) - mu * N .* ycomponent .* sin(phi);
% plot(T,fII),xlabel('时间t / s','Interpreter','tex'),ylabel('侧视图中沿纸面摩擦力f_{II}','Interpreter','tex');
% grid on;
% 
% figure(7);
% plot(T,N),xlabel('时间t/s','Interpreter','tex'),ylabel('N','Interpreter','tex');
% grid on;