function y = exactSol04(t)
    w = 200;
    k1 = 1/150;
    k2 = 4/150;
    %H = 1000;
    g = 32;
    t1 = 5;
    y = zeros(length(t),2);

    idx1 = find(t <= t1);
    idx2 = find(t > t1);
    
    % 0 <= t <= t1
    T1 = t(idx1); 
    y(idx1,1) = 1/g * w/k1 .* log(cosh(g * (k1/w)^0.5 .* T1));
    y(idx1,2) = (w/k1)^0.5 .* tanh(g * (k1/w)^0.5 .* T1); 

    % t1 < t
    T2 = t(idx2);
    C = w/g * (1/k1 - 1/k2) * log(cosh(g * (k1/w)^0.5 * t1));
    y(idx2,1) = 1/g * w/k2 * log(cosh(g/sqrt(w) * (sqrt(k2)*(T2-t1) + sqrt(k1)*t1))) + C;
    y(idx2,2) = (w/k2)^0.5 * tanh(g/sqrt(w) * (sqrt(k2)*(T2-t1) + sqrt(k1)*t1));
end