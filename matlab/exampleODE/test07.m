function dydt = test07(t,y)
% IC: y1(0) = 3, y2(0) = 0, y3(0) = 0 
    %assert((0 <= t) && (t <= 2*pi), "0 <= t <= 2pi");
    dydt = [0; 0; 0];
    dydt(1) = -y(2)-y(1)*y(2) / sqrt(y(1)^2 + y(2)^2);
    dydt(2) = y(1) - y(2)*y(3) / sqrt(y(1)^2 + y(2)^2);
    dydt(3) = y(1) / sqrt(y(1)^2 + y(2)^2);
end