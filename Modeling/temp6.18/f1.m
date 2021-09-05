function dy = f1(t,y)
    % Formalized normal differential equation system of the beginning (without slipping)
    g = 9.7964;
    mu = 0.25;
    m = 32.93 * 10^(-3);
    R = 10 * 10^(-3);
    vp = 60;
    S = pi * R^2;
    p00 = -0.003193;
    p01 = -0.001244;
    p02 = -8.636e-05;
    p03 = 2.055e-06;
    p10 = 4.146;
    p12 = 0.008178;
    dy = zeros(10,1);
    phi = y(3);
%     if  (y(4) > pi / 2) && (y(9) > 0)
%         y(4) = pi / 2;
%         dy(4) = -0.5 * y(9);
%     else
%         dy(4) = y(9);
%     end
    theta = y(4);
    psi = y(5);
    x1 = y(6);
    y1 = y(7);
    phi1 = y(8);
    theta1 = y(9);
    psi1 = y(10);
    vx = x1 + R * phi1 * sin(theta) * sin(phi) - R * psi1 * sin(theta) * sin(phi) - 2 * R * theta1 * cos(theta / 2)^2 * cos(phi);
    vy = y1 - R * phi1 * sin(theta) * cos(phi) + R * psi1 * sin(theta) * cos(phi) - 2 * R * theta1 * cos(theta / 2)^2 * sin(phi);
%     if (abs(vx) > 10^(-3)) || (abs(vy) > 10^(-3))
%         xcomponent = vx / sqrt(vx^2 + vy^2);
%         ycomponent = vy / sqrt(vx^2 + vy^2);
%     else
        xcomponent = 0;
        ycomponent = 0;
%     end
    fp = p00 + p10 * S + p01 * (vp - R * phi1 * sin(theta)) + p02 * (vp - R * phi1 * sin(theta))^2 + p12 * S * (vp - R * phi1 * sin(theta))^2 + p03 * (vp - R * phi1 * sin(theta))^3;
    fr = p00 + p10 * S + p01 * (R * phi1 * sin(theta)) + p02 * (R * phi1 * sin(theta))^2 + p12 * S * (R * phi1 * sin(theta))^2 + p03 * (R * phi1 * sin(theta))^3;
    dy(1) = y(6);
    dy(2) = y(7);
    dy(3) = y(8);
    dy(4) = y(9);
    dy(5) = y(10);
    dy(6) = (mu*xcomponent*(5*fp*cos(phi)*cos(theta)*sin(theta) - 7*g*m + 7*R*m*theta1^2*cos(theta) + 5*R*m*phi1^2*cos(theta)*sin(theta)^2 - 2*R*m*phi1*psi1*sin(theta)^2))/(m*(5*sin(theta)^2 + 10*mu*xcomponent*cos(theta/2)^2*cos(phi)*sin(theta) + 10*mu*ycomponent*cos(theta/2)^2*sin(phi)*sin(theta) + 7));
    dy(7) = (mu*ycomponent*(5*fp*cos(phi)*cos(theta)*sin(theta) - 7*g*m + 7*R*m*theta1^2*cos(theta) + 5*R*m*phi1^2*cos(theta)*sin(theta)^2 - 2*R*m*phi1*psi1*sin(theta)^2))/(m*(5*sin(theta)^2 + 10*mu*xcomponent*cos(theta/2)^2*cos(phi)*sin(theta) + 10*mu*ycomponent*cos(theta/2)^2*sin(phi)*sin(theta) + 7));
    dy(8) = -(35*fr*sin(phi) - 35*R*fp*cos(phi) + 25*fr*sin(phi)*sin(theta)^2 - 25*R*fp*cos(phi)*sin(theta)^2 - 14*R^2*m*psi1*theta1 + 84*R^2*m*phi1*theta1*cos(theta) - 10*R^2*m*psi1*theta1*sin(theta)^2 + 50*fr*mu*ycomponent*cos(theta/2)^2*sin(phi)^2*sin(theta) + 60*R^2*m*phi1*theta1*cos(theta)*sin(theta)^2 - 50*R*fp*mu*xcomponent*cos(theta/2)^2*cos(phi)^2*sin(theta) - 70*R*g*m*mu*ycomponent*cos(theta/2)^2*cos(phi) + 50*fr*mu*xcomponent*cos(theta/2)^2*cos(phi)*sin(phi)*sin(theta) + 70*R*g*m*mu*xcomponent*cos(theta/2)^2*sin(phi) - 50*R*fp*mu*ycomponent*cos(theta/2)^2*cos(phi)*sin(phi)*sin(theta) + 70*R^2*m*mu*theta1^2*ycomponent*cos(theta/2)^2*cos(phi)*cos(theta) + 50*R*fp*mu*ycomponent*cos(theta/2)^2*cos(phi)^2*cos(theta)*sin(theta) - 70*R^2*m*mu*theta1^2*xcomponent*cos(theta/2)^2*cos(theta)*sin(phi) + 20*R^2*m*mu*phi1*psi1*xcomponent*cos(theta/2)^2*sin(phi)*sin(theta)^2 - 50*R*fp*mu*xcomponent*cos(theta/2)^2*cos(phi)*cos(theta)*sin(phi)*sin(theta) - 20*R^2*m*mu*psi1*theta1*xcomponent*cos(theta/2)^2*cos(phi)*sin(theta) - 20*R^2*m*mu*psi1*theta1*ycomponent*cos(theta/2)^2*sin(phi)*sin(theta) + 50*R^2*m*mu*phi1^2*ycomponent*cos(theta/2)^2*cos(phi)*cos(theta)*sin(theta)^2 - 50*R^2*m*mu*phi1^2*xcomponent*cos(theta/2)^2*cos(theta)*sin(phi)*sin(theta)^2 - 20*R^2*m*mu*phi1*psi1*ycomponent*cos(theta/2)^2*cos(phi)*sin(theta)^2 + 120*R^2*m*mu*phi1*theta1*xcomponent*cos(theta/2)^2*cos(phi)*cos(theta)*sin(theta) + 120*R^2*m*mu*phi1*theta1*ycomponent*cos(theta/2)^2*cos(theta)*sin(phi)*sin(theta))/(7*R^2*m*sin(theta)*(5*sin(theta)^2 + 10*mu*xcomponent*cos(theta/2)^2*cos(phi)*sin(theta) + 10*mu*ycomponent*cos(theta/2)^2*sin(phi)*sin(theta) + 7));
    dy(9) = (5*g*m*sin(theta) + 5*fp*cos(phi)*cos(theta) + 5*R*m*phi1^2*cos(theta)*sin(theta) - 5*R*m*theta1^2*cos(theta)*sin(theta) - 2*R*m*phi1*psi1*sin(theta) + 10*g*m*mu*xcomponent*cos(theta/2)^2*cos(phi) + 10*g*m*mu*ycomponent*cos(theta/2)^2*sin(phi) - 10*R*m*mu*theta1^2*xcomponent*cos(theta/2)^2*cos(phi)*cos(theta) - 10*R*m*mu*theta1^2*ycomponent*cos(theta/2)^2*cos(theta)*sin(phi))/(R*m*(5*sin(theta)^2 + 10*mu*xcomponent*cos(theta/2)^2*cos(phi)*sin(theta) + 10*mu*ycomponent*cos(theta/2)^2*sin(phi)*sin(theta) + 7));
    dy(10) = (70*fr*cos(theta)*sin(phi) - 70*R*fp*cos(phi)*cos(theta) + 50*fr*cos(theta)*sin(phi)*sin(theta)^2 - 50*R*fp*cos(phi)*cos(theta)*sin(theta)^2 - 28*R^2*m*psi1*theta1*cos(theta) + 168*R^2*m*phi1*theta1*cos(theta)^2 + 98*R^2*m*phi1*theta1*sin(theta)^2 + 70*R^2*m*phi1*theta1*sin(theta)^4 - 20*R^2*m*psi1*theta1*cos(theta)*sin(theta)^2 + 120*R^2*m*phi1*theta1*cos(theta)^2*sin(theta)^2 + 100*fr*mu*ycomponent*cos(theta/2)^2*cos(theta)*sin(phi)^2*sin(theta) + 175*R*fp*mu*ycomponent*cos(phi)^2*cos(theta)*sin(theta)^3 - 245*R*g*m*mu*ycomponent*cos(phi)*sin(theta)^2 + 245*R*g*m*mu*xcomponent*sin(phi)*sin(theta)^2 + 175*R^2*m*mu*phi1^2*ycomponent*cos(phi)*cos(theta)*sin(theta)^4 + 245*R^2*m*mu*theta1^2*ycomponent*cos(phi)*cos(theta)*sin(theta)^2 - 175*R^2*m*mu*phi1^2*xcomponent*cos(theta)*sin(phi)*sin(theta)^4 - 245*R^2*m*mu*theta1^2*xcomponent*cos(theta)*sin(phi)*sin(theta)^2 - 70*R^2*m*mu*phi1*psi1*ycomponent*cos(phi)*sin(theta)^4 + 70*R^2*m*mu*phi1*psi1*xcomponent*sin(phi)*sin(theta)^4 - 100*R*fp*mu*xcomponent*cos(theta/2)^2*cos(phi)^2*cos(theta)*sin(theta) - 140*R*g*m*mu*ycomponent*cos(theta/2)^2*cos(phi)*cos(theta) + 100*fr*mu*xcomponent*cos(theta/2)^2*cos(phi)*cos(theta)*sin(phi)*sin(theta) + 140*R*g*m*mu*xcomponent*cos(theta/2)^2*cos(theta)*sin(phi) - 175*R*fp*mu*xcomponent*cos(phi)*cos(theta)*sin(phi)*sin(theta)^3 + 140*R^2*m*mu*theta1^2*ycomponent*cos(theta/2)^2*cos(phi)*cos(theta)^2 + 100*R*fp*mu*ycomponent*cos(theta/2)^2*cos(phi)^2*cos(theta)^2*sin(theta) - 140*R^2*m*mu*theta1^2*xcomponent*cos(theta/2)^2*cos(theta)^2*sin(phi) + 140*R^2*m*mu*phi1*theta1*ycomponent*cos(theta/2)^2*sin(phi)*sin(theta)^3 - 100*R*fp*mu*ycomponent*cos(theta/2)^2*cos(phi)*cos(theta)*sin(phi)*sin(theta) + 100*R^2*m*mu*phi1^2*ycomponent*cos(theta/2)^2*cos(phi)*cos(theta)^2*sin(theta)^2 - 100*R^2*m*mu*phi1^2*xcomponent*cos(theta/2)^2*cos(theta)^2*sin(phi)*sin(theta)^2 - 100*R*fp*mu*xcomponent*cos(theta/2)^2*cos(phi)*cos(theta)^2*sin(phi)*sin(theta) + 140*R^2*m*mu*phi1*theta1*xcomponent*cos(theta/2)^2*cos(phi)*sin(theta)^3 - 40*R^2*m*mu*psi1*theta1*xcomponent*cos(theta/2)^2*cos(phi)*cos(theta)*sin(theta) - 40*R^2*m*mu*psi1*theta1*ycomponent*cos(theta/2)^2*cos(theta)*sin(phi)*sin(theta) - 40*R^2*m*mu*phi1*psi1*ycomponent*cos(theta/2)^2*cos(phi)*cos(theta)*sin(theta)^2 + 240*R^2*m*mu*phi1*theta1*xcomponent*cos(theta/2)^2*cos(phi)*cos(theta)^2*sin(theta) + 40*R^2*m*mu*phi1*psi1*xcomponent*cos(theta/2)^2*cos(theta)*sin(phi)*sin(theta)^2 + 240*R^2*m*mu*phi1*theta1*ycomponent*cos(theta/2)^2*cos(theta)^2*sin(phi)*sin(theta))/(14*R^2*m*sin(theta)*(5*sin(theta)^2 + 10*mu*xcomponent*cos(theta/2)^2*cos(phi)*sin(theta) + 10*mu*ycomponent*cos(theta/2)^2*sin(phi)*sin(theta) + 7));
end