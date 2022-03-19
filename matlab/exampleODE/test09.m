function dydt = test09(t,y)
% IC: y1(0) = -4
    assert((1 <= t), "1 <= t");
    dydt = 2*y / t + 5;
end