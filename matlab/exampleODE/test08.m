function dydt = test08(t,y)
% IC: y1(0) = 0.4, y2(0) = 0, y3(0) = 0, y4(0) = 2 
    %assert((0 <= t) && (t <= 2*pi), "0 <= t <= 2pi");
    dydt = [0; 0; 0; 0];
    dydt(1) = y(2);
    dydt(2) = -y(1)/(y(1)^2 + y(3)^2)^1.5;
    dydt(3) = y(4);
    dydt(4) = -y(3)/(y(1)^2 + y(3)^2)^1.5;
end