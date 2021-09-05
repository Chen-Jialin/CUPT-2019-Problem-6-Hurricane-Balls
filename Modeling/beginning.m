% calculation and show the begining motion (with slipping) of the hurricane balls
close all;clear;clc;
T_max = 3;
[T,Y] = ode45('f1',[0,T_max],[0,0,0,pi / 2,0,0,0,25,0,0]);
for i = 1:size(Y,1)
    if Y(i,4) > (pi / 2)
        Y(i,4) = pi / 2;
    end
end
option =0;
g = 9.7964;
mu = 0.29;
m = 56.60 * 10^(-3);
R = 12 * 10^(-3);
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
N = -(2*(7*g*m - 7*R*m*theta1.^2.*cos(theta) - 5*R*m*phi1.^2.*cos(theta).*sin(theta).^2 + 2*R*m*phi1.*psi1.*sin(theta).^2))./(5*mu*xcomponent.*cos(phi).*sin(theta) - 5*sin(theta).^2 + 5*mu*ycomponent.*sin(phi).*sin(theta) + 5*mu*xcomponent.*cos(phi).*cos(theta).*sin(theta) + 5*mu*ycomponent.*cos(theta).*sin(phi).*sin(theta) - 7);
if option == 1
    [X,Y,Z] = sphere(16);
    [X,Y,Z] = multiplyEuMat(EuMat(0,pi / 2,0),X,Y,Z);
    [X1,Y1,Z1] = multiplyEuMat(EuMat(phi(1),pi / 2 - theta(1),psi(1)),X,Y,Z);
    [X2,Y2,Z2] = multiplyEuMat(EuMat(phi(1),pi / 2 - theta(1),psi(1)),X,Y,Z);
    X1 = R * X1;Y1 = R * (Y1 - 1);Z1 = R * (Z1 + 1);
    X2 = R * X2;Y2 = R * (Y2 + 1);Z2 = R * (Z2 + 1);
    S1 = surf(X1,Y1,Z1);hold on;
    view(90,0);
    material shiny;
    S2 = surf(X2,Y2,Z2);hold off;
    material shiny;
    xlabel('x');ylabel('y');zlabel('z');
    axis([-3 * R;3 * R;-3 * R;3 * R;-3 * R;3 * R]);
    axis square;
    step = 0.002;
    I = [];
    i = 0;
    for t = 0:step:T_max
        i = i + 1;
        [temp,I(i)] = min(abs(T - t));
    end
    for i = I(2:end)
        X0 = x(i);Y0 = y(i);Z0 = R * (1 + cos(theta(i)));
        X1 = R * X;Y1 = R * (Y + 1);Z1 = R * Z;
        X2 = R * X;Y2 = R * (Y - 1);Z2 = R * Z;
        [X1,Y1,Z1] = multiplyEuMat(EuMat(phi(i),pi / 2 - theta(i),psi(i)),X1,Y1,Z1);
        [X2,Y2,Z2] = multiplyEuMat(EuMat(phi(i),pi / 2 - theta(i),psi(i)),X2,Y2,Z2);
        X1 = X1 + X0;Y1 = Y1 + Y0;Z1 = Z1 + Z0;
        X2 = X2 + X0;Y2 = Y2 + Y0;Z2 = Z2 + Z0;
        set(S1,'XData',X1,'YData',Y1,'ZData',Z1);
        set(S2,'XData',X2,'YData',Y2,'ZData',Z2);
        
        drawnow;
        pause((T(i) - T(i - 1)) * 10);
        disp(T(i))
    end
    
else
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
    legend('公转角速度','自转角速度');
    
    figure(4);
    plot(T,Y(:,4),'LineWidth',2);
    xlabel('时间t / s','Interpreter','tex'),ylabel('两球心连线与竖直方向夹角\theta','Interpreter','tex');
    grid on;
    
    figure(5);
    fI = mu * N .* xcomponent .* sin(phi) - mu * N .* ycomponent .* cos(phi);
    plot(T,fI,'LineWidth',2);
    xlabel('时间t / s','Interpreter','tex');
    ylabel('侧视图中垂直纸面摩擦力f_{I}','Interpreter','tex');
    grid on;
    
    figure(6);
    fII = -mu * N .* xcomponent .* cos(phi) - mu * N .* ycomponent .* sin(phi);
    plot(T,fII),xlabel('时间t / s','Interpreter','tex'),ylabel('侧视图中沿纸面摩擦力f_{II}','Interpreter','tex');
    grid on;
    
    figure(7);
    plot(T,N),xlabel('时间t/s','Interpreter','tex'),ylabel('N','Interpreter','tex');
    grid on;
end