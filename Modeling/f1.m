function dy = f1(t,y)
    % Formalized normal differential equation system of the beginning (without slipping)
    g = 9.7964;
    mu1 = 0.29;
    m = 56.60 * 10^(-3);
    R = 12 * 10^(-3);
    dy = zeros(10,1);
    phi = y(3);
    if  (y(4) > pi / 2) && (y(9) > 0)
        y(4) = pi / 2;
        dy(4) = -0.2 * y(9);
    else
        dy(4) = y(9);
    end
    theta = y(4);
    psi = y(5);
    x1 = y(6);
    y1 = y(7);
    phi1 = y(8);
    theta1 = y(9);
    psi1 = y(10);
    vx = x1 + R * phi1 * sin(theta) * sin(phi) - R * psi1 * sin(theta) * sin(phi) - 2 * R * theta1 * cos(theta / 2)^2 * cos(phi);
    vy = y1 - R * phi1 * sin(theta) * cos(phi) + R * psi1 * sin(theta) * cos(phi) - 2 * R * theta1 * cos(theta / 2)^2 * sin(phi);
    if (vx > 2 * 10^(-3)) || (vy > 2 * 10^(-3))
        xcomponent = vx / sqrt(vx^2 + vy^2);
        ycomponent = vy / sqrt(vx^2 + vy^2);
    else
        xcomponent = 0;
        ycomponent = 0;
    end
    dy(1) = y(6);
    dy(2) = y(7);
    dy(3) = y(8);
%     dy(4) = y(9);
    dy(5) = y(10);
    dy(6) = -(mu1*xcomponent*(7*g - 7*R*theta1^2*cos(theta) - 5*R*phi1^2*cos(theta)*sin(theta)^2 + 2*R*phi1*psi1*sin(theta)^2))/(5*sin(theta)^2 + 5*mu1*xcomponent*cos(phi)*sin(theta) + 5*mu1*ycomponent*sin(phi)*sin(theta) + 5*mu1*xcomponent*cos(phi)*cos(theta)*sin(theta) + 5*mu1*ycomponent*cos(theta)*sin(phi)*sin(theta) + 7);
    dy(7) = -(mu1*ycomponent*(7*g - 7*R*theta1^2*cos(theta) - 5*R*phi1^2*cos(theta)*sin(theta)^2 + 2*R*phi1*psi1*sin(theta)^2))/(5*sin(theta)^2 + 5*mu1*xcomponent*cos(phi)*sin(theta) + 5*mu1*ycomponent*sin(phi)*sin(theta) + 5*mu1*xcomponent*cos(phi)*cos(theta)*sin(theta) + 5*mu1*ycomponent*cos(theta)*sin(phi)*sin(theta) + 7);
    dy(8) = (2*(7*R*psi1*theta1 - 42*R*phi1*theta1*cos(theta) + 5*R*psi1*theta1*sin(theta)^2 + 35*g*mu1*ycomponent*cos(theta/2)^2*cos(phi) - 35*g*mu1*xcomponent*cos(theta/2)^2*sin(phi) - 30*R*phi1*theta1*cos(theta)*sin(theta)^2 + 5*R*mu1*psi1*theta1*xcomponent*cos(phi)*sin(theta) + 5*R*mu1*psi1*theta1*ycomponent*sin(phi)*sin(theta) - 35*R*mu1*theta1^2*ycomponent*cos(theta/2)^2*cos(phi)*cos(theta) + 35*R*mu1*theta1^2*xcomponent*cos(theta/2)^2*cos(theta)*sin(phi) + 25*R*mu1*phi1^2*xcomponent*cos(theta/2)^2*cos(theta)*sin(phi)*sin(theta)^2 - 30*R*mu1*phi1*theta1*xcomponent*cos(phi)*cos(theta)*sin(theta) + 5*R*mu1*psi1*theta1*xcomponent*cos(phi)*cos(theta)*sin(theta) - 30*R*mu1*phi1*theta1*ycomponent*cos(theta)*sin(phi)*sin(theta) + 5*R*mu1*psi1*theta1*ycomponent*cos(theta)*sin(phi)*sin(theta) + 10*R*mu1*phi1*psi1*ycomponent*cos(theta/2)^2*cos(phi)*sin(theta)^2 - 10*R*mu1*phi1*psi1*xcomponent*cos(theta/2)^2*sin(phi)*sin(theta)^2 - 30*R*mu1*phi1*theta1*xcomponent*cos(phi)*cos(theta)^2*sin(theta) - 30*R*mu1*phi1*theta1*ycomponent*cos(theta)^2*sin(phi)*sin(theta) - 25*R*mu1*phi1^2*ycomponent*cos(theta/2)^2*cos(phi)*cos(theta)*sin(theta)^2))/(7*R*sin(theta)*(5*sin(theta)^2 + 5*mu1*xcomponent*cos(phi)*sin(theta) + 5*mu1*ycomponent*sin(phi)*sin(theta) + 5*mu1*xcomponent*cos(phi)*cos(theta)*sin(theta) + 5*mu1*ycomponent*cos(theta)*sin(phi)*sin(theta) + 7));
    dy(9) = (5*g*sin(theta) + 5*R*phi1^2*cos(theta)*sin(theta) - 5*R*theta1^2*cos(theta)*sin(theta) - 2*R*phi1*psi1*sin(theta) + 5*g*mu1*xcomponent*cos(phi) + 5*g*mu1*ycomponent*sin(phi) + 5*g*mu1*xcomponent*cos(phi)*cos(theta) + 5*g*mu1*ycomponent*cos(theta)*sin(phi) - 5*R*mu1*theta1^2*xcomponent*cos(phi)*cos(theta) - 5*R*mu1*theta1^2*ycomponent*cos(theta)*sin(phi) - 5*R*mu1*theta1^2*xcomponent*cos(phi)*cos(theta)^2 - 5*R*mu1*theta1^2*ycomponent*cos(theta)^2*sin(phi))/(R*(5*sin(theta)^2 + 5*mu1*xcomponent*cos(phi)*sin(theta) + 5*mu1*ycomponent*sin(phi)*sin(theta) + 5*mu1*xcomponent*cos(phi)*cos(theta)*sin(theta) + 5*mu1*ycomponent*cos(theta)*sin(phi)*sin(theta) + 7));
    dy(10) = (168*R*phi1*theta1*cos(theta)^2 - 28*R*psi1*theta1*cos(theta) + 98*R*phi1*theta1*sin(theta)^2 + 70*R*phi1*theta1*sin(theta)^4 + 120*R*phi1*theta1*cos(theta)^2*sin(theta)^2 - 20*R*psi1*theta1*cos(theta)*sin(theta)^2 - 245*g*mu1*ycomponent*cos(phi)*sin(theta)^2 + 245*g*mu1*xcomponent*sin(phi)*sin(theta)^2 - 140*g*mu1*ycomponent*cos(theta/2)^2*cos(phi)*cos(theta) + 140*g*mu1*xcomponent*cos(theta/2)^2*cos(theta)*sin(phi) + 140*R*mu1*theta1^2*ycomponent*cos(theta/2)^2*cos(phi)*cos(theta)^2 - 140*R*mu1*theta1^2*xcomponent*cos(theta/2)^2*cos(theta)^2*sin(phi) + 175*R*mu1*phi1^2*ycomponent*cos(phi)*cos(theta)*sin(theta)^4 + 245*R*mu1*theta1^2*ycomponent*cos(phi)*cos(theta)*sin(theta)^2 - 175*R*mu1*phi1^2*xcomponent*cos(theta)*sin(phi)*sin(theta)^4 - 245*R*mu1*theta1^2*xcomponent*cos(theta)*sin(phi)*sin(theta)^2 - 70*R*mu1*phi1*psi1*ycomponent*cos(phi)*sin(theta)^4 + 70*R*mu1*phi1*theta1*xcomponent*cos(phi)*sin(theta)^3 + 70*R*mu1*phi1*psi1*xcomponent*sin(phi)*sin(theta)^4 + 70*R*mu1*phi1*theta1*ycomponent*sin(phi)*sin(theta)^3 - 20*R*mu1*psi1*theta1*xcomponent*cos(phi)*cos(theta)*sin(theta) - 20*R*mu1*psi1*theta1*ycomponent*cos(theta)*sin(phi)*sin(theta) + 100*R*mu1*phi1^2*ycomponent*cos(theta/2)^2*cos(phi)*cos(theta)^2*sin(theta)^2 - 100*R*mu1*phi1^2*xcomponent*cos(theta/2)^2*cos(theta)^2*sin(phi)*sin(theta)^2 + 120*R*mu1*phi1*theta1*xcomponent*cos(phi)*cos(theta)^2*sin(theta) + 70*R*mu1*phi1*theta1*xcomponent*cos(phi)*cos(theta)*sin(theta)^3 + 120*R*mu1*phi1*theta1*xcomponent*cos(phi)*cos(theta)^3*sin(theta) - 20*R*mu1*psi1*theta1*xcomponent*cos(phi)*cos(theta)^2*sin(theta) + 120*R*mu1*phi1*theta1*ycomponent*cos(theta)^2*sin(phi)*sin(theta) + 70*R*mu1*phi1*theta1*ycomponent*cos(theta)*sin(phi)*sin(theta)^3 + 120*R*mu1*phi1*theta1*ycomponent*cos(theta)^3*sin(phi)*sin(theta) - 20*R*mu1*psi1*theta1*ycomponent*cos(theta)^2*sin(phi)*sin(theta) - 40*R*mu1*phi1*psi1*ycomponent*cos(theta/2)^2*cos(phi)*cos(theta)*sin(theta)^2 + 40*R*mu1*phi1*psi1*xcomponent*cos(theta/2)^2*cos(theta)*sin(phi)*sin(theta)^2)/(14*R*sin(theta)*(5*sin(theta)^2 + 5*mu1*xcomponent*cos(phi)*sin(theta) + 5*mu1*ycomponent*sin(phi)*sin(theta) + 5*mu1*xcomponent*cos(phi)*cos(theta)*sin(theta) + 5*mu1*ycomponent*cos(theta)*sin(phi)*sin(theta) + 7));
end