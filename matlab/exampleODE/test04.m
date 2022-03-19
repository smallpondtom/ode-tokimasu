function dydt = test04(t,y)
% IC: y1(0) = 0, y2(0) = 0
    w = 200;
    k1 = 1/150;
    k2 = 4/150;
    %H = 1000;
    g = 32;
    m = w / g;
    t1 = 5;
    dydt = [0, 0]';

    assert((0 <= t)&&(t <= 20),"0 <= t <= 20");

    dydt(1) = y(2);
    if t <= t1
        dydt(2) = (w - k1*y(2)^2)/m;
    else
        dydt(2) = (w - k2*y(2)^2)/m;
    end
end