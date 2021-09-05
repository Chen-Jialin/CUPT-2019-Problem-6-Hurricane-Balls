function dy = f1(t,y)
    % formalized normal differential equation system of stable rotation (without slipping)
    m = 33.93 * 10^(-3);
    g = 9.7964;
    R = 0.01;
    S = pi * R^2;
%     mu = 0.0011;
    vp = 20;
    p00 = -0.003193;
    p01 = 0.001244;
    p02 = -8.636e-05;
    p03 = 2.055e-06;
    p10 = 4.146;
    p12 = 0.008178;
    dy = zeros(10,1);
%     x = y(1);
%     y = y(2);
    phi = y(3);
    theta = y(4);
    psi = y(5);
%     x1 = y(6);
%     y1 = y(7);
    phi1 = y(8);
    theta1 = y(9);
    psi1 = y(10);
    fp = (p00 + p10 * S + p01 * (vp - R * phi1 * sin(theta))...
        + p02 * (vp - R * phi1 * sin(theta))^2 ...
        + p12 * S * (vp - R * phi1 * sin(theta))^2 ...
        + p03 * (vp - R * phi1 * sin(theta))^3);
    fr = p00 + p10 * S + p01 * (R * phi1 * sin(theta))...
        + p02 * (R * phi1 * sin(theta))^2 ...
        + p12 * S * (R * phi1 * sin(theta))^2 ...
        + p03 * (R * phi1 * sin(theta))^3;
    dy(1) = y(6);
    dy(2) = y(7);
    dy(3) = y(8);
    dy(4) = y(9);
    dy(5) = y(10);
    dy(6) = -(175*fp*sin(phi)^3*sin(theta)^4 + 50*fp*sin(phi)*sin(theta)^2 - 100*fr*sin(phi)*sin(theta)^2 + 70*R*fp*sin(phi) - 140*R*fr*sin(phi) + 175*fp*cos(phi)^2*sin(phi)*sin(theta)^4 + 245*R*fp*sin(phi)^3*sin(theta)^2 + 50*fp*cos(theta)*sin(phi)*sin(theta)^2 - 100*fr*cos(theta)*sin(phi)*sin(theta)^2 + 70*R*fp*cos(theta)*sin(phi) - 140*R*fr*cos(theta)*sin(phi) + 100*fp*cos(theta/2)^2*sin(phi)^3*sin(theta)^2 + 140*R*fp*cos(theta/2)^2*sin(phi)^3 + 200*R*fp*cos(theta/2)^4*sin(phi)^3 + 400*R*fp*cos(theta/2)^6*sin(phi)^5 - 400*R*fr*cos(theta/2)^4*sin(phi)^3 + 196*R^2*m*phi1^2*cos(phi)*sin(theta) + 140*R*m*phi1^2*cos(phi)*sin(theta)^3 + 196*R^2*m*theta1^2*cos(phi)*sin(theta) + 140*R*m*theta1^2*cos(phi)*sin(theta)^3 + 800*R*fp*cos(theta/2)^6*cos(phi)^2*sin(phi)^3 + 700*R*fp*cos(theta/2)^4*sin(phi)^5*sin(theta)^2 + 56*R^2*m*psi1*theta1*sin(phi) + 350*R*m*phi1^2*cos(phi)^3*sin(theta)^5 + 350*R*m*theta1^2*cos(phi)^3*sin(theta)^5 + 245*R*fp*cos(phi)^2*sin(phi)*sin(theta)^2 + 100*fp*cos(theta/2)^2*cos(phi)^2*sin(phi)*sin(theta)^2 + 100*fp*cos(theta/2)^2*cos(theta)*sin(phi)^3*sin(theta)^2 + 140*R*fp*cos(theta/2)^2*cos(phi)^2*sin(phi) + 200*R*fp*cos(theta/2)^4*cos(phi)^2*sin(phi) + 400*R*fp*cos(theta/2)^6*cos(phi)^4*sin(phi) - 400*R*fr*cos(theta/2)^4*cos(phi)^2*sin(phi) + 490*R^2*m*phi1^2*cos(phi)^3*sin(theta)^3 + 140*R*fp*cos(theta/2)^2*cos(theta)*sin(phi)^3 + 200*R*fp*cos(theta/2)^4*cos(theta)*sin(phi)^3 + 400*R*fp*cos(theta/2)^6*cos(theta)*sin(phi)^5 - 400*R*fr*cos(theta/2)^4*cos(theta)*sin(phi)^3 + 490*R^2*m*theta1^2*cos(phi)^3*sin(theta)^3 + 280*R^2*m*theta1^2*cos(theta/2)^2*cos(phi)^3*sin(theta) + 200*R*m*theta1^2*cos(theta/2)^2*cos(phi)^3*sin(theta)^3 + 700*R*fp*cos(theta/2)^4*cos(phi)^4*sin(phi)*sin(theta)^2 - 196*R^2*m*phi1*psi1*cos(phi)*sin(theta) - 140*R*m*phi1*psi1*cos(phi)*sin(theta)^3 - 140*R^2*m*phi1*theta1*cos(theta)*sin(phi) - 140*R^2*m*psi1*theta1*cos(theta)*sin(phi) - 140*R*m*phi1*theta1*sin(phi)*sin(theta)^4 + 40*R*m*psi1*theta1*sin(phi)*sin(theta)^2 - 400*R*g*m*cos(theta/2)^4*cos(phi)^3*sin(theta) + 1120*R^2*m*phi1*theta1*cos(theta/2)^6*sin(phi)^3 + 160*R^2*m*psi1*theta1*cos(theta/2)^4*sin(phi)^3 + 350*R*m*phi1^2*cos(phi)*sin(phi)^2*sin(theta)^5 + 350*R*m*theta1^2*cos(phi)*sin(phi)^2*sin(theta)^5 + 1400*R*fp*cos(theta/2)^4*cos(phi)^2*sin(phi)^3*sin(theta)^2 - 350*R*m*phi1*psi1*cos(phi)^3*sin(theta)^5 - 336*R^2*m*phi1*theta1*cos(theta)^2*sin(phi) - 196*R^2*m*phi1*theta1*sin(phi)*sin(theta)^2 + 100*fp*cos(theta/2)^2*cos(phi)^2*cos(theta)*sin(phi)*sin(theta)^2 - 700*R*g*m*cos(theta/2)^2*cos(phi)^3*sin(theta)^3 + 140*R*fp*cos(theta/2)^2*cos(phi)^2*cos(theta)*sin(phi) + 200*R*fp*cos(theta/2)^4*cos(phi)^2*cos(theta)*sin(phi) + 400*R*fp*cos(theta/2)^6*cos(phi)^4*cos(theta)*sin(phi) - 400*R*fr*cos(theta/2)^4*cos(phi)^2*cos(theta)*sin(phi) + 490*R^2*m*phi1^2*cos(phi)*sin(phi)^2*sin(theta)^3 + 490*R^2*m*theta1^2*cos(phi)*sin(phi)^2*sin(theta)^3 - 280*R*g*m*cos(theta/2)^2*cos(phi)*sin(theta) - 490*R^2*m*phi1*psi1*cos(phi)^3*sin(theta)^3 + 392*R^2*m*phi1*theta1*cos(theta/2)^2*sin(phi) + 800*R*fp*cos(theta/2)^6*cos(phi)^2*cos(theta)*sin(phi)^3 + 280*R^2*m*phi1^2*cos(theta/2)^2*cos(phi)^3*sin(theta) + 200*R*m*phi1^2*cos(theta/2)^2*cos(phi)^3*sin(theta)^3 - 100*R*m*phi1*theta1*cos(theta)*sin(phi)*sin(theta)^2 - 100*R*m*psi1*theta1*cos(theta)*sin(phi)*sin(theta)^2 - 400*R*g*m*cos(theta/2)^4*cos(phi)^3*cos(theta)*sin(theta) + 1120*R^2*m*phi1*theta1*cos(theta/2)^4*cos(phi)^2*sin(phi) - 1120*R^2*m*phi1*theta1*cos(theta/2)^6*cos(phi)^2*sin(phi) + 160*R^2*m*psi1*theta1*cos(theta/2)^4*cos(phi)^2*sin(phi) - 280*R^2*m*phi1*psi1*cos(theta/2)^2*cos(phi)^3*sin(theta) - 200*R*m*phi1*psi1*cos(theta/2)^2*cos(phi)^3*sin(theta)^3 + 160*R^2*m*phi1*psi1*cos(theta/2)^4*cos(phi)^3*sin(theta) - 400*R^2*m*phi1*theta1*cos(theta/2)^4*cos(theta)*sin(phi)^3 - 400*R^2*m*psi1*theta1*cos(theta/2)^4*cos(theta)*sin(phi)^3 - 400*R*g*m*cos(theta/2)^4*cos(phi)*sin(phi)^2*sin(theta) + 280*R*m*theta1^2*cos(theta/2)^2*cos(phi)*cos(theta)*sin(theta) - 700*R^2*m*phi1^2*cos(theta/2)^2*cos(phi)^3*cos(theta)*sin(theta)^3 - 400*R^2*m*phi1^2*cos(theta/2)^4*cos(phi)^3*cos(theta)^2*sin(theta) - 350*R*m*phi1*psi1*cos(phi)*sin(phi)^2*sin(theta)^5 - 240*R*m*phi1*theta1*cos(theta)^2*sin(phi)*sin(theta)^2 + 280*R^2*m*phi1*psi1*cos(theta/2)^2*cos(phi)^3*sin(theta)^3 - 960*R^2*m*phi1*theta1*cos(theta/2)^4*cos(theta)^2*sin(phi)^3 - 280*R^2*m*phi1^2*cos(theta/2)^2*cos(phi)*cos(theta)*sin(theta) - 700*R*g*m*cos(theta/2)^2*cos(phi)*sin(phi)^2*sin(theta)^3 + 400*R*m*theta1^2*cos(theta/2)^4*cos(phi)^3*cos(theta)*sin(theta) - 560*R^2*m*phi1*theta1*cos(theta/2)^4*sin(phi)^3*sin(theta)^2 + 112*R^2*m*phi1*psi1*cos(theta/2)^2*cos(phi)*sin(theta) - 490*R^2*m*phi1*psi1*cos(phi)*sin(phi)^2*sin(theta)^3 + 280*R*m*phi1*theta1*cos(theta/2)^2*sin(phi)*sin(theta)^2 + 280*R^2*m*phi1^2*cos(theta/2)^2*cos(phi)^3*cos(theta)*sin(theta) + 200*R*m*phi1^2*cos(theta/2)^2*cos(phi)^3*cos(theta)*sin(theta)^3 - 400*R^2*m*phi1^2*cos(theta/2)^4*cos(phi)^3*cos(theta)*sin(theta) + 280*R^2*m*theta1^2*cos(theta/2)^2*cos(phi)^3*cos(theta)*sin(theta) + 900*R*m*theta1^2*cos(theta/2)^2*cos(phi)^3*cos(theta)*sin(theta)^3 + 400*R*m*theta1^2*cos(theta/2)^4*cos(phi)^3*cos(theta)^2*sin(theta) + 280*R^2*m*phi1^2*cos(theta/2)^2*cos(phi)*sin(phi)^2*sin(theta) + 200*R*m*phi1^2*cos(theta/2)^2*cos(phi)*sin(phi)^2*sin(theta)^3 + 280*R^2*m*theta1^2*cos(theta/2)^2*cos(phi)*sin(phi)^2*sin(theta) + 200*R*m*theta1^2*cos(theta/2)^2*cos(phi)*sin(phi)^2*sin(theta)^3 - 700*R^2*m*phi1^2*cos(theta/2)^2*cos(phi)*cos(theta)*sin(phi)^2*sin(theta)^3 - 400*R^2*m*phi1^2*cos(theta/2)^4*cos(phi)*cos(theta)^2*sin(phi)^2*sin(theta) - 960*R^2*m*phi1*theta1*cos(theta/2)^4*cos(phi)^2*cos(theta)^2*sin(phi) + 280*R^2*m*phi1*psi1*cos(theta/2)^2*cos(phi)*sin(phi)^2*sin(theta)^3 + 1960*R^2*m*phi1*theta1*cos(theta/2)^2*cos(phi)^2*sin(phi)*sin(theta)^2 - 560*R^2*m*phi1*theta1*cos(theta/2)^4*cos(phi)^2*sin(phi)*sin(theta)^2 + 400*R*m*theta1^2*cos(theta/2)^4*cos(phi)*cos(theta)*sin(phi)^2*sin(theta) + 280*R^2*m*phi1^2*cos(theta/2)^2*cos(phi)*cos(theta)*sin(phi)^2*sin(theta) + 200*R*m*phi1^2*cos(theta/2)^2*cos(phi)*cos(theta)*sin(phi)^2*sin(theta)^3 - 400*R^2*m*phi1^2*cos(theta/2)^4*cos(phi)*cos(theta)*sin(phi)^2*sin(theta) + 280*R^2*m*theta1^2*cos(theta/2)^2*cos(phi)*cos(theta)*sin(phi)^2*sin(theta) + 900*R*m*theta1^2*cos(theta/2)^2*cos(phi)*cos(theta)*sin(phi)^2*sin(theta)^3 + 400*R*m*theta1^2*cos(theta/2)^4*cos(phi)*cos(theta)^2*sin(phi)^2*sin(theta) + 720*R^2*m*phi1*theta1*cos(theta/2)^4*cos(phi)^2*cos(theta)*sin(phi) - 400*R^2*m*psi1*theta1*cos(theta/2)^4*cos(phi)^2*cos(theta)*sin(phi) - 280*R^2*m*phi1*psi1*cos(theta/2)^2*cos(phi)^3*cos(theta)*sin(theta) - 200*R*m*phi1*psi1*cos(theta/2)^2*cos(phi)^3*cos(theta)*sin(theta)^3 + 160*R^2*m*phi1*psi1*cos(theta/2)^4*cos(phi)^3*cos(theta)*sin(theta) - 400*R*g*m*cos(theta/2)^4*cos(phi)*cos(theta)*sin(phi)^2*sin(theta) - 280*R^2*m*phi1*psi1*cos(theta/2)^2*cos(phi)*sin(phi)^2*sin(theta) - 200*R*m*phi1*psi1*cos(theta/2)^2*cos(phi)*sin(phi)^2*sin(theta)^3 + 160*R^2*m*phi1*psi1*cos(theta/2)^4*cos(phi)*sin(phi)^2*sin(theta) + 1400*R*m*phi1*theta1*cos(theta/2)^2*cos(phi)^2*sin(phi)*sin(theta)^4 + 800*R*m*phi1*theta1*cos(theta/2)^4*cos(phi)^2*sin(phi)*sin(theta)^2 - 280*R^2*m*phi1*psi1*cos(theta/2)^2*cos(phi)*cos(theta)*sin(phi)^2*sin(theta) - 200*R*m*phi1*psi1*cos(theta/2)^2*cos(phi)*cos(theta)*sin(phi)^2*sin(theta)^3 + 160*R^2*m*phi1*psi1*cos(theta/2)^4*cos(phi)*cos(theta)*sin(phi)^2*sin(theta) + 800*R*m*phi1*theta1*cos(theta/2)^4*cos(phi)^2*cos(theta)*sin(phi)*sin(theta)^2)/(2*m*(7*R + 5*sin(theta)^2 + 20*R*cos(theta/2)^4*cos(phi)^2 + 20*R*cos(theta/2)^4*sin(phi)^2)*(20*cos(theta/2)^2*cos(phi)^2 + 20*cos(theta/2)^2*sin(phi)^2 + 35*cos(phi)^2*sin(theta)^2 + 35*sin(phi)^2*sin(theta)^2 + 20*cos(theta/2)^2*cos(phi)^2*cos(theta) + 20*cos(theta/2)^2*cos(theta)*sin(phi)^2 + 14));
    dy(7) = -(100*fr*cos(phi)*sin(theta)^2 - 50*fp*cos(phi)*sin(theta)^2 - 175*fp*cos(phi)^3*sin(theta)^4 - 70*R*fp*cos(phi) + 140*R*fr*cos(phi) - 175*fp*cos(phi)*sin(phi)^2*sin(theta)^4 - 245*R*fp*cos(phi)^3*sin(theta)^2 - 50*fp*cos(phi)*cos(theta)*sin(theta)^2 + 100*fr*cos(phi)*cos(theta)*sin(theta)^2 - 70*R*fp*cos(phi)*cos(theta) + 140*R*fr*cos(phi)*cos(theta) - 100*fp*cos(theta/2)^2*cos(phi)^3*sin(theta)^2 - 140*R*fp*cos(theta/2)^2*cos(phi)^3 - 200*R*fp*cos(theta/2)^4*cos(phi)^3 - 400*R*fp*cos(theta/2)^6*cos(phi)^5 + 400*R*fr*cos(theta/2)^4*cos(phi)^3 + 196*R^2*m*phi1^2*sin(phi)*sin(theta) + 140*R*m*phi1^2*sin(phi)*sin(theta)^3 + 196*R^2*m*theta1^2*sin(phi)*sin(theta) + 140*R*m*theta1^2*sin(phi)*sin(theta)^3 - 800*R*fp*cos(theta/2)^6*cos(phi)^3*sin(phi)^2 - 700*R*fp*cos(theta/2)^4*cos(phi)^5*sin(theta)^2 - 56*R^2*m*psi1*theta1*cos(phi) - 245*R*fp*cos(phi)*sin(phi)^2*sin(theta)^2 + 350*R*m*phi1^2*sin(phi)^3*sin(theta)^5 + 350*R*m*theta1^2*sin(phi)^3*sin(theta)^5 - 100*fp*cos(theta/2)^2*cos(phi)^3*cos(theta)*sin(theta)^2 - 140*R*fp*cos(theta/2)^2*cos(phi)^3*cos(theta) - 200*R*fp*cos(theta/2)^4*cos(phi)^3*cos(theta) - 400*R*fp*cos(theta/2)^6*cos(phi)^5*cos(theta) + 400*R*fr*cos(theta/2)^4*cos(phi)^3*cos(theta) - 100*fp*cos(theta/2)^2*cos(phi)*sin(phi)^2*sin(theta)^2 - 140*R*fp*cos(theta/2)^2*cos(phi)*sin(phi)^2 - 200*R*fp*cos(theta/2)^4*cos(phi)*sin(phi)^2 - 400*R*fp*cos(theta/2)^6*cos(phi)*sin(phi)^4 + 400*R*fr*cos(theta/2)^4*cos(phi)*sin(phi)^2 + 490*R^2*m*phi1^2*sin(phi)^3*sin(theta)^3 + 490*R^2*m*theta1^2*sin(phi)^3*sin(theta)^3 - 700*R*fp*cos(theta/2)^4*cos(phi)*sin(phi)^4*sin(theta)^2 + 140*R^2*m*phi1*theta1*cos(phi)*cos(theta) + 140*R^2*m*psi1*theta1*cos(phi)*cos(theta) + 280*R^2*m*phi1^2*cos(theta/2)^2*sin(phi)^3*sin(theta) + 200*R*m*phi1^2*cos(theta/2)^2*sin(phi)^3*sin(theta)^3 + 280*R^2*m*theta1^2*cos(theta/2)^2*sin(phi)^3*sin(theta) + 200*R*m*theta1^2*cos(theta/2)^2*sin(phi)^3*sin(theta)^3 + 140*R*m*phi1*theta1*cos(phi)*sin(theta)^4 - 40*R*m*psi1*theta1*cos(phi)*sin(theta)^2 - 196*R^2*m*phi1*psi1*sin(phi)*sin(theta) - 140*R*m*phi1*psi1*sin(phi)*sin(theta)^3 + 1120*R^2*m*phi1*theta1*cos(theta/2)^6*cos(phi)^3 - 160*R^2*m*psi1*theta1*cos(theta/2)^4*cos(phi)^3 - 400*R*g*m*cos(theta/2)^4*sin(phi)^3*sin(theta) + 350*R*m*phi1^2*cos(phi)^2*sin(phi)*sin(theta)^5 + 350*R*m*theta1^2*cos(phi)^2*sin(phi)*sin(theta)^5 - 1400*R*fp*cos(theta/2)^4*cos(phi)^3*sin(phi)^2*sin(theta)^2 + 336*R^2*m*phi1*theta1*cos(phi)*cos(theta)^2 + 196*R^2*m*phi1*theta1*cos(phi)*sin(theta)^2 - 350*R*m*phi1*psi1*sin(phi)^3*sin(theta)^5 - 100*fp*cos(theta/2)^2*cos(phi)*cos(theta)*sin(phi)^2*sin(theta)^2 - 140*R*fp*cos(theta/2)^2*cos(phi)*cos(theta)*sin(phi)^2 - 200*R*fp*cos(theta/2)^4*cos(phi)*cos(theta)*sin(phi)^2 - 400*R*fp*cos(theta/2)^6*cos(phi)*cos(theta)*sin(phi)^4 + 400*R*fr*cos(theta/2)^4*cos(phi)*cos(theta)*sin(phi)^2 - 700*R*g*m*cos(theta/2)^2*sin(phi)^3*sin(theta)^3 + 490*R^2*m*phi1^2*cos(phi)^2*sin(phi)*sin(theta)^3 + 490*R^2*m*theta1^2*cos(phi)^2*sin(phi)*sin(theta)^3 + 392*R^2*m*phi1*theta1*cos(theta/2)^2*cos(phi) - 280*R*g*m*cos(theta/2)^2*sin(phi)*sin(theta) - 490*R^2*m*phi1*psi1*sin(phi)^3*sin(theta)^3 - 800*R*fp*cos(theta/2)^6*cos(phi)^3*cos(theta)*sin(phi)^2 + 400*R^2*m*phi1*theta1*cos(theta/2)^4*cos(phi)^3*cos(theta) + 400*R^2*m*psi1*theta1*cos(theta/2)^4*cos(phi)^3*cos(theta) + 1120*R^2*m*phi1*theta1*cos(theta/2)^4*cos(phi)*sin(phi)^2 - 1120*R^2*m*phi1*theta1*cos(theta/2)^6*cos(phi)*sin(phi)^2 - 160*R^2*m*psi1*theta1*cos(theta/2)^4*cos(phi)*sin(phi)^2 - 400*R*g*m*cos(theta/2)^4*cos(phi)^2*sin(phi)*sin(theta) - 400*R*g*m*cos(theta/2)^4*cos(theta)*sin(phi)^3*sin(theta) - 280*R^2*m*phi1*psi1*cos(theta/2)^2*sin(phi)^3*sin(theta) - 200*R*m*phi1*psi1*cos(theta/2)^2*sin(phi)^3*sin(theta)^3 + 160*R^2*m*phi1*psi1*cos(theta/2)^4*sin(phi)^3*sin(theta) + 280*R*m*theta1^2*cos(theta/2)^2*cos(theta)*sin(phi)*sin(theta) - 700*R^2*m*phi1^2*cos(theta/2)^2*cos(theta)*sin(phi)^3*sin(theta)^3 - 400*R^2*m*phi1^2*cos(theta/2)^4*cos(theta)^2*sin(phi)^3*sin(theta) + 240*R*m*phi1*theta1*cos(phi)*cos(theta)^2*sin(theta)^2 - 350*R*m*phi1*psi1*cos(phi)^2*sin(phi)*sin(theta)^5 + 960*R^2*m*phi1*theta1*cos(theta/2)^4*cos(phi)^3*cos(theta)^2 + 560*R^2*m*phi1*theta1*cos(theta/2)^4*cos(phi)^3*sin(theta)^2 - 700*R*g*m*cos(theta/2)^2*cos(phi)^2*sin(phi)*sin(theta)^3 + 280*R^2*m*phi1*psi1*cos(theta/2)^2*sin(phi)^3*sin(theta)^3 - 280*R^2*m*phi1^2*cos(theta/2)^2*cos(theta)*sin(phi)*sin(theta) + 400*R*m*theta1^2*cos(theta/2)^4*cos(theta)*sin(phi)^3*sin(theta) + 280*R*m*phi1*theta1*cos(theta/2)^2*cos(phi)*sin(theta)^2 - 490*R^2*m*phi1*psi1*cos(phi)^2*sin(phi)*sin(theta)^3 + 112*R^2*m*phi1*psi1*cos(theta/2)^2*sin(phi)*sin(theta) + 280*R^2*m*phi1^2*cos(theta/2)^2*cos(phi)^2*sin(phi)*sin(theta) + 200*R*m*phi1^2*cos(theta/2)^2*cos(phi)^2*sin(phi)*sin(theta)^3 + 280*R^2*m*theta1^2*cos(theta/2)^2*cos(phi)^2*sin(phi)*sin(theta) + 200*R*m*theta1^2*cos(theta/2)^2*cos(phi)^2*sin(phi)*sin(theta)^3 + 280*R^2*m*phi1^2*cos(theta/2)^2*cos(theta)*sin(phi)^3*sin(theta) + 200*R*m*phi1^2*cos(theta/2)^2*cos(theta)*sin(phi)^3*sin(theta)^3 - 400*R^2*m*phi1^2*cos(theta/2)^4*cos(theta)*sin(phi)^3*sin(theta) + 280*R^2*m*theta1^2*cos(theta/2)^2*cos(theta)*sin(phi)^3*sin(theta) + 900*R*m*theta1^2*cos(theta/2)^2*cos(theta)*sin(phi)^3*sin(theta)^3 + 400*R*m*theta1^2*cos(theta/2)^4*cos(theta)^2*sin(phi)^3*sin(theta) + 100*R*m*phi1*theta1*cos(phi)*cos(theta)*sin(theta)^2 + 100*R*m*psi1*theta1*cos(phi)*cos(theta)*sin(theta)^2 - 700*R^2*m*phi1^2*cos(theta/2)^2*cos(phi)^2*cos(theta)*sin(phi)*sin(theta)^3 - 400*R^2*m*phi1^2*cos(theta/2)^4*cos(phi)^2*cos(theta)^2*sin(phi)*sin(theta) + 960*R^2*m*phi1*theta1*cos(theta/2)^4*cos(phi)*cos(theta)^2*sin(phi)^2 + 280*R^2*m*phi1*psi1*cos(theta/2)^2*cos(phi)^2*sin(phi)*sin(theta)^3 + 1960*R^2*m*phi1*theta1*cos(theta/2)^2*cos(phi)*sin(phi)^2*sin(theta)^2 + 560*R^2*m*phi1*theta1*cos(theta/2)^4*cos(phi)*sin(phi)^2*sin(theta)^2 + 400*R*m*theta1^2*cos(theta/2)^4*cos(phi)^2*cos(theta)*sin(phi)*sin(theta) + 280*R^2*m*phi1^2*cos(theta/2)^2*cos(phi)^2*cos(theta)*sin(phi)*sin(theta) + 200*R*m*phi1^2*cos(theta/2)^2*cos(phi)^2*cos(theta)*sin(phi)*sin(theta)^3 - 400*R^2*m*phi1^2*cos(theta/2)^4*cos(phi)^2*cos(theta)*sin(phi)*sin(theta) + 280*R^2*m*theta1^2*cos(theta/2)^2*cos(phi)^2*cos(theta)*sin(phi)*sin(theta) + 900*R*m*theta1^2*cos(theta/2)^2*cos(phi)^2*cos(theta)*sin(phi)*sin(theta)^3 + 400*R*m*theta1^2*cos(theta/2)^4*cos(phi)^2*cos(theta)^2*sin(phi)*sin(theta) + 1520*R^2*m*phi1*theta1*cos(theta/2)^4*cos(phi)*cos(theta)*sin(phi)^2 + 400*R^2*m*psi1*theta1*cos(theta/2)^4*cos(phi)*cos(theta)*sin(phi)^2 - 400*R*g*m*cos(theta/2)^4*cos(phi)^2*cos(theta)*sin(phi)*sin(theta) - 280*R^2*m*phi1*psi1*cos(theta/2)^2*cos(phi)^2*sin(phi)*sin(theta) - 200*R*m*phi1*psi1*cos(theta/2)^2*cos(phi)^2*sin(phi)*sin(theta)^3 + 160*R^2*m*phi1*psi1*cos(theta/2)^4*cos(phi)^2*sin(phi)*sin(theta) + 1400*R*m*phi1*theta1*cos(theta/2)^2*cos(phi)*sin(phi)^2*sin(theta)^4 + 800*R*m*phi1*theta1*cos(theta/2)^4*cos(phi)*sin(phi)^2*sin(theta)^2 - 280*R^2*m*phi1*psi1*cos(theta/2)^2*cos(theta)*sin(phi)^3*sin(theta) - 200*R*m*phi1*psi1*cos(theta/2)^2*cos(theta)*sin(phi)^3*sin(theta)^3 + 160*R^2*m*phi1*psi1*cos(theta/2)^4*cos(theta)*sin(phi)^3*sin(theta) - 280*R^2*m*phi1*psi1*cos(theta/2)^2*cos(phi)^2*cos(theta)*sin(phi)*sin(theta) - 200*R*m*phi1*psi1*cos(theta/2)^2*cos(phi)^2*cos(theta)*sin(phi)*sin(theta)^3 + 160*R^2*m*phi1*psi1*cos(theta/2)^4*cos(phi)^2*cos(theta)*sin(phi)*sin(theta) + 800*R*m*phi1*theta1*cos(theta/2)^4*cos(phi)*cos(theta)*sin(phi)^2*sin(theta)^2)/(2*m*(7*R + 5*sin(theta)^2 + 20*R*cos(theta/2)^4*cos(phi)^2 + 20*R*cos(theta/2)^4*sin(phi)^2)*(20*cos(theta/2)^2*cos(phi)^2 + 20*cos(theta/2)^2*sin(phi)^2 + 35*cos(phi)^2*sin(theta)^2 + 35*sin(phi)^2*sin(theta)^2 + 20*cos(theta/2)^2*cos(phi)^2*cos(theta) + 20*cos(theta/2)^2*cos(theta)*sin(phi)^2 + 14));
    dy(8) = (5*g*sin(theta) - 5*theta1^2*cos(theta)*sin(theta) + 5*R*phi1^2*cos(theta)*sin(theta) - 2*R*phi1*psi1*sin(theta) + 10*R*phi1^2*cos(theta/2)^2*cos(phi)^2*sin(theta) + 10*R*theta1^2*cos(theta/2)^2*cos(phi)^2*sin(theta) + 10*R*phi1^2*cos(theta/2)^2*sin(phi)^2*sin(theta) + 10*R*theta1^2*cos(theta/2)^2*sin(phi)^2*sin(theta) - 10*R*phi1*psi1*cos(theta/2)^2*cos(phi)^2*sin(theta) - 10*R*phi1*psi1*cos(theta/2)^2*sin(phi)^2*sin(theta) + 40*R*phi1*theta1*cos(theta/2)^4*cos(phi)*sin(phi))/(7*R + 5*sin(theta)^2 + 20*R*cos(theta/2)^4*cos(phi)^2 + 20*R*cos(theta/2)^4*sin(phi)^2);
    dy(10) = -(5*fp*cos(theta) - 10*fr*cos(theta) + 5*fp*cos(phi)^2*sin(theta)^2 + 25*fr*cos(phi)^2*sin(theta)^2 + 5*fp*sin(phi)^2*sin(theta)^2 + 25*fr*sin(phi)^2*sin(theta)^2 + 10*fp*cos(theta/2)^2*cos(phi)^2*cos(theta) + 10*fp*cos(theta/2)^2*cos(theta)*sin(phi)^2 - 24*R*m*phi1*theta1*cos(theta)^2 - 14*R*m*phi1*theta1*sin(theta)^2 + 4*R*m*psi1*theta1*cos(theta) - 10*R*m*psi1*theta1*cos(phi)^2*sin(theta)^2 - 10*R*m*psi1*theta1*sin(phi)^2*sin(theta)^2 - 20*R*m*phi1*theta1*cos(theta/2)^2*cos(phi)^2*cos(theta)^2 + 20*R*m*psi1*theta1*cos(theta/2)^2*cos(phi)^2*cos(theta)^2 + 50*R*m*phi1*theta1*cos(theta/2)^2*cos(phi)^2*sin(theta)^2 - 20*R*m*phi1*theta1*cos(theta/2)^2*cos(theta)^2*sin(phi)^2 + 20*R*m*psi1*theta1*cos(theta/2)^2*cos(theta)^2*sin(phi)^2 - 90*R*m*phi1*theta1*cos(theta/2)^2*sin(phi)^2*sin(theta)^2 + 25*R*m*phi1*theta1*cos(phi)^2*cos(theta)*sin(theta)^2 + 35*R*m*psi1*theta1*cos(phi)^2*cos(theta)*sin(theta)^2 + 25*R*m*phi1*theta1*cos(theta)*sin(phi)^2*sin(theta)^2 + 35*R*m*psi1*theta1*cos(theta)*sin(phi)^2*sin(theta)^2 + 40*R*m*phi1*theta1*cos(theta/2)^4*cos(phi)^2*cos(theta) - 40*R*m*phi1*theta1*cos(theta/2)^4*cos(theta)*sin(phi)^2)/(R*m*sin(theta)*(20*cos(theta/2)^2*cos(phi)^2 + 20*cos(theta/2)^2*sin(phi)^2 + 35*cos(phi)^2*sin(theta)^2 + 35*sin(phi)^2*sin(theta)^2 + 20*cos(theta/2)^2*cos(phi)^2*cos(theta) + 20*cos(theta/2)^2*cos(theta)*sin(phi)^2 + 14));
end