function J = jacobian07(t,y)
% This if for the ESDIRK23 algorithm
    J = zeros(3,3);
    
    K = (y(1)^2 + y(2))^1.5;
    L = sqrt(y(1)^2 + y(2)^2);
    
    J(1,1) = y(1)^2*y(2)/K - y(2)/L;
    J(1,2) = y(1)*y(2)^2/K - y(1)/L - 1;
    J(1,3) = 0;
    J(2,1) = y(1)*y(2)*y(3)/K + 1;
    J(2,2) = y(2)^2*y(3)/K - y(3)/L;
    J(2,3) = -y(2)/L;
    J(3,1) = 1/L - y(1)^2/K;
    J(3,2) = -y(1)*y(2)/K;
    J(3,3) = 0;
end