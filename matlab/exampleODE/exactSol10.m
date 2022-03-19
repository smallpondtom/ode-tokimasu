function y = exactSol10(t)
    N = length(t);
    y = zeros(N,7);
    y(:,1) = exp(-t);
    y(:,2) = 4./(t + 2).^2;
    y(:,3) = exp(sin(t));
    y(:,4) = 20./(1 + 19*exp(-t/4));
    y(:,5) = 1 + (exp(-t) + exp(-3*t))/2;
    y(:,6) = 1 - exp(-3*t);
    y(:,7) = 1 - (exp(-t) - exp(-3*t))/2;
end