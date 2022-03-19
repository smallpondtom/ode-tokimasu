function dydt = test06(t,y)
% IC: y1(1) = 1
    assert((1 <= t) && (t <= 12), "1 <= t <= 12");
    dydt = ((2 * log(y) + 8) / t - 5)*y;
end