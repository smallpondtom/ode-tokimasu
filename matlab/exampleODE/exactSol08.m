function y = exactSol08(t)
    N = length(t);
    y = zeros(N,4);
    
    i = 1;
    for ti = t
        u = solve_transcendental(ti, 0.6, 1e-3);
        y(i,1) = cos(u) - 0.6;
        y(i,2) = -sin(u) / (1 - 0.6*cos(u));
        y(i,3) = 0.8*sin(u);
        y(i,4) = 0.8*cos(u)/(1 - 0.6*cos(u));
        i = i + 1;
    end
end


function u = solve_transcendental(t, c, tol)
 %{
       NAME:    SOLVE_TRANSCENDENTAL
       AUTHOR:  TOMOKI KOIKE
       INPUTS:  (1) t:    TIME  
                (2) c:    COEFFICIENT
                (3) tol:  TOLERANCE
       OUTPUTS: (1) u:    OUTPUT
       DESCRIPTION: Solving the transcendental equation u - c*sin(u) = t.
 %}

    %Checking for user inputed tolerance
    if nargin == 2
        %using default value
        tol = 10^-8;
    elseif nargin > 3
        error('Too many inputs.')
    elseif nargin < 2
        error('Too few inputs.')
    end
   
    
    % Check tolerance 
    if tol > 0.01
        error('Set a tolerance smaller than 0.001')
    end 
    
    % Bessel function method 
    del1 = 1; del2 = 1; del3 = 1;
    N = 1;
    ustore = [];
    while (N < 10)
        f = 0;
        for m = 1:N
            f = f + (1 / m) * besselj(m,m*c) * sin(m*t);
        end
        ustore = [ustore, t + 2*f];
        if N > 4
            del1 = abs(ustore(N) - ustore(N-1));
            del2 = abs(ustore(N-1) - ustore(N-2));
            del3 = abs(ustore(N-2) - ustore(N-3));

            if (del1 < tol && del2 < tol && del3 < tol)
                break;
            end
        end
        N = N + 1;
    end
    u = ustore(end);
end