function y = exactSol07(t)
    N = length(t);
    y = zeros(N,3);
    y(:,1) = cos(t) .* (2 + cos(t));
    y(:,2) = sin(t) .* (2 + cos(t));
    y(:,3) = sin(t);
end