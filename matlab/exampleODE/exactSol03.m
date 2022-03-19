function y = exactSol03(t)
    y = zeros(length(t),2);
    y(:,1) = (sin(t)+2)./(-cos(t)+2);
    y(:,2) = -cos(t)+2;
end