function dydt = test13(t,y)
% IC: y1(0) = 1, y2(0) = 0.5
    dydt = [0;0];
    
    dydt(1) = y(2);
    dydt(2) = -4*t*y(2) - (2 + 4*t^2)*y(1);
end