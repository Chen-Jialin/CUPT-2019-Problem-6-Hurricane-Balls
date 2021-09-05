% calculation and show stable rotation (without rotation) of the hurricane balls
close all;clear;clc;
T_max = 10;
Omega = 20 * pi;
theta = acos(2 / 5 - 9.7964 / Omega^2 / 0.01);
[T,Y] = ode23('f1',[0,T_max],[0,0,0,theta,0,0,0,Omega,0,Omega]);
figure(1)
plot(T,Y(:,4))
figure(2)
plot(T,Y(:,8))
figure(3)
plot(T,Y(:,9))