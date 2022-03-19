function y = exactSol11(t)
    N = length(t);
    y = zeros(N,5);
    y(:,1) = t.^2/2;
    y(:,2) = t;
    y(:,3) = t.^2 .* (1 + t/3)/2;
    y(:,4) = t.^3 .* (1 + t/4)/6;
    y(:,5) = t .* (t-6) .* (t - 4);
end