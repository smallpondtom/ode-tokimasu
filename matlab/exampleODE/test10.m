function dydt = test10(t,y)
% IC: y1(0) = 1, y2(0) = 1, y3(0) = 1, y4(0) = 1
%     y5(0) = 2, y6(0) = 0, y7(0) = 1
    assert((0 <= t), "0 <= t");

    dydt = zeros(7,1);
    dydt(1) = -y(1);
    dydt(2) = -y(2)^1.5;
    dydt(3) = y(3)*cos(t);
    dydt(4) = y(4)*(1 - y(4)/20) / 4;
    dydt(5) = -y(5) + y(6);
    dydt(6) = y(5) - 2*y(6) + y(7);
    dydt(7) = y(6) - y(7);
end