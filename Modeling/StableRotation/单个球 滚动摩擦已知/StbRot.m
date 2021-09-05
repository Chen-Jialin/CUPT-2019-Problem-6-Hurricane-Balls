% calculation and show stable rotation (without rotation) of the hurricane balls
close all;clear;clc;
T_max = 10;
Omega = 30 * pi;
theta = acos(2 / 5 - 9.7964 / Omega^2 / 0.01);
[T,Y] = ode45('f1',[0,T_max],[0,0,0,theta,0,0,0,Omega,0,Omega]);
figure(1)
plot(T,Y(:,4))
figure(2)
plot(T,Y(:,8))
hold on;
plot(T,Y(:,10))