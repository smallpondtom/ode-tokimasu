function Mass = MassFcn1(t)
Mass = zeros(4);
Mass(1,1) = 1;
Mass(1,2) = 1e-7*t;
Mass(2,2) = 2;
Mass(3,3) = 3;
Mass(4,4) = 4;