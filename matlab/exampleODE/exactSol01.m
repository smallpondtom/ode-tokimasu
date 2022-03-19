function y = exactSol01(t)
    y = zeros(length(t),2);
    y(:,1) = 4*(t + 1/8*exp(-8*t)-1/8);
    y(:,2) = 4*(1 - exp(-8*t));
end