function dydt = test03(t,y)
% IC: y1(0) = 2, y2(0) = 1
    dydt = [(cos(t) - sin(t).*y(1))./y(2); sin(t)];
end