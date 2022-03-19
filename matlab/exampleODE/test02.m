function dydt = test02(t,y)
% IC: y1(0) = 0.5
    dydt = [0];
    dydt(1) = -y(1)*cos(t);
end