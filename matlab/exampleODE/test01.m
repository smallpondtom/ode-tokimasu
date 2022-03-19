function dydt = test01(t,y)
% IC: y1(0) = 0, y2(0) = 0
    m = 0.25;
    w = 8;
    k = 2;
    H = 10;
    dydt = [0, 0]';

    assert((0 <= t)&&(t <= H/2),"0 <= t <= H/2");

    dydt(1) = y(2);
    dydt(2) = (w - k*y(2))/m;
end