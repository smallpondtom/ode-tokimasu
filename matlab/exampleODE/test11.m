function dydt = test11(t,y)
% IC: y1(-1) = 0.5, y2(-1) = -1, y3(-1) = 1/3, y4(-1) = -1/8, y5(-1) = -35,
    assert((-1 <= t), "-1 <= t");

    dydt = zeros(5,1);
    dydt(1) = y(2);
    dydt(2) = 1;
    dydt(3) = y(1)+y(2);
    dydt(4) = y(3);
    dydt(5) = 3*t^2 - 20*t + 24;
end