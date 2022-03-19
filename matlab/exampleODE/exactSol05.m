function y = exactSol05(t)
    y = zeros(length(t),2);
    y(:,1) = exp(-16/pi * sin(pi*t/2)) .* cos(8*pi*t);
    y(:,2) = -8*exp(-16/pi * sin(pi*t/2))  .* (pi*sin(8*pi*t) + cos(pi*t/2).*cos(8*pi*t));
end