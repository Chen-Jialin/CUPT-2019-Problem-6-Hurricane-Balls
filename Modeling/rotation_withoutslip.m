% calculation and show the begining motion (with slipping) of the hurricane balls
close all;clear;clc;
[T,Y] = ode45('f1',[0,1],[0,0,0,pi / 2,0,0,0,2 * pi * 20,0,0]);

figure(1);
subplot(1,2,1),plot(Y(:,1),Y(:,2)),xlabel('$x$','Interpreter','latex'),ylabel('y','Interpreter','latex'),axis equal;
subplot(1,2,2),plot(T,Y(:,6)),hold on;
plot(T,Y(:,7)),hold off,xlabel('T','Interpreter','latex'),ylabel('v','Interpreter','latex'),legend('$\dot{x}$','$\dot{y}$','Interpreter','latex');

figure(2);
plot(T,Y(:,8)),hold on;
plot(T,Y(:,10)),hold off,xlabel('$t$','Interpreter','latex'),ylabel('$\omega$','Interpreter','latex'),legend('$\dot{phi}$','$\dot{psi}$','Interpreter','latex')

figure(3);
subplot(1,2,1),plot(T,Y(:,4));
subplot(1,2,2),plot(T,Y(:,9));

figure(4);
plot(T,N),xlabel('t','Interpreter','latex'),ylabel('N','Interpreter','latex')