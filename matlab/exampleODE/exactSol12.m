function y = exactSol12(t)
    N = length(t);

    epsilon = 1e-6;
    y(:,1) = t ./ (epsilon + t.^2).^0.5;
    y(:,2) = epsilon ./ (epsilon + t.^2).^1.5;
end