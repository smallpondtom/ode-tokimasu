function dydt = test05(t,y)
% y1(0) = 1,  y2(0) = -8
    dydt = [0; 0];
    dydt(1) = y(2);
    dydt(2) = -16*cos(pi*t/2)*y(2) - (64*pi^2 + 64*cos(pi*t/2)^2 - 4*pi*sin(pi*t/2))*y(1);
end