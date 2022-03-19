function dydt = test03(t,y)
    dydt = [(cos(t) - sin(t).*y(1))./y(2); sin(t)];
end