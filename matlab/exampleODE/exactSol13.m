function y = exactSol13(t)
    y(:,1) = 0.5*exp(-1 * t.^2) .* (t + 2);
    y(:,2) = (-t.^2 - 2*t + 0.5) .* exp(-1 * t.^2);
end