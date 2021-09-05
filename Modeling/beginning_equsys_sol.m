% Formalization of the equation system of the beginning process
clear;clc;
syms x y phi theta psi x1 y1 phi1 theta1 psi1 x2 y2 phi2 theta2 psi2 m g mu1 R N xcomponent ycomponent
f = mu1 * N;
fx = -f * xcomponent;
fy = -f * ycomponent;
I1 = 14 / 5 * m * R^2;
I2 = 14 / 5 * m * R^2;
I3 = 4 / 5 * m * R^2;
omega1 = phi1 * sin(theta) * sin(psi) + theta1 * cos(psi);
omega2 = phi1 * sin(theta) * cos(psi) - theta1 * sin(psi);
omega3 = phi1 * cos(theta) + psi1;
omega11 = phi2 * sin(theta) * sin(psi) + phi1 * theta1 * cos(theta) * sin(psi) + phi1 * psi1 * sin(theta) * cos(psi) + theta2 * cos(psi) - theta1 * psi1 * sin(psi);
omega21 = phi2 * sin(theta) * cos(psi) + phi1 * theta1 * cos(theta) * cos(psi) - phi1 * psi1 * sin(theta) * sin(psi) - theta2 * sin(psi) - theta1 * psi1 * cos(psi);
omega31 = phi2 * cos(theta) - phi1 * theta1 * sin(theta) + psi2;
M10 = -R * (fx * cos(phi) + fy * sin(phi)) * (1 + cos(theta)) + R * N * sin(theta);
M20 = 2 * R * (fx * sin(phi) - fy * cos(phi)) * cos(theta / 2)^2;
M1 = M10 * cos(psi) + M20 * sin(psi);
M2 = -M10 * sin(psi) + M20 * cos(psi);
M3 = -R * (fx * sin(phi) - fy * cos(phi)) * sin(theta);
equ1 = 2 * m * x2 - fx;
equ2 = 2 * m * y2 - fy;
equ3 = 2 * m * R * (-theta2 * sin(theta) - theta1^2 * cos(theta)) - N + 2 * m * g;
equ4 = I1 * omega11 - (I2 - I3) * omega2 * omega3 - M1;
equ5 = I2 * omega21 - (I3 - I1) * omega3 * omega1 - M2;
equ6 = I3 * omega31 - (I1 - I2) * omega1 * omega2 - M3;
[x2,y2,phi2,theta2,psi2,N] = solve(equ1,equ2,equ3,equ4,equ5,equ6,x2,y2,phi2,theta2,psi2,N);