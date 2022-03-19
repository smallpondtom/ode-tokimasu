function dydt = test12(t,y)
% IC: y1(-0.1) = -0.1/(epsilon + 0.01)^0.5
%     y2(-0.1) = epsilon / (epsilon + 0.01)^0.5
% 
% where epsilon = 1e-6
    assert((-0.1 <= t) && (t <= 0.1), "-0.1 <= t <= 0.1");
    dydt = [0;0];
    
    epsilon = 1e-6;
    dydt(1) = y(2);
    dydt(2) = (-3*epsilon/(epsilon + t^2)^2)*y(1);
end