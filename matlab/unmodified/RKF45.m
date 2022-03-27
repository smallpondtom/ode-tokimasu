function [t,x,E] = RKF45(fun,tspan,x0,AbsTol,RelTol,varargin)
% Solves x' = f(t,x) using the Runge -Kutta - Fehlberg method
%
% INPUTS
% fun : function handle
% tspan : span of time in which to approximate solution
% or time vector in which to approximate solution
% x0 : initial conditions
% AbsTol : Absolute tolerance
% RelTol : Relative tolerance
%
% OUTPUTS
% t : time vector
% x : approximation vector

    %% INITIALIZATION
    % RKF45 Butcher Tableau
    A = [0, 0, 0, 0, 0, 0;
         1/4, 0, 0, 0, 0, 0;
         3/32, 9/32, 0, 0, 0, 0;
         1932/2197, -7200/2197, 7296/2197, 0, 0, 0;
         439/216, -8, 3680/513, -845/4104, 0, 0;
         -8/27, 2, -3544/2565, 1859/4104, -11/40, 0];
    b = [25/216, 0, 1408/2565, 2197/4104, -1/5, 0];
    bhat = [16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55];
    c = [0,1/4,3/8,12/13,1,1/2];
    d = b - bhat;

    % Various parameters
    epsilon = 0.8;
    p = 4;  % order
    phat = p+1;
    kp = 0.4/phat;
    kI = 0.3/phat;
    h = 1e-3;  % Initial step size if using step size controller
    fixedstepsize = length(tspan) > 2;
    x = zeros(length(tspan),length(x0));
    x(1,:) = x0;

    %% INTEGRATE
    if(fixedstepsize)
        h = (tspan(end) - tspan(1))/(length(tspan) -1);
        t = tspan;
        for i = 1:length(t)-1
            % RKF45 steps
            fn = feval(fun,t(i),x(i,:),varargin{:})';
            x(i+1,:) = RKF45Step(fun,fn,t(i),x(i,:),h,A,bhat,c,d,varargin{:});
        end
    else
        t = zeros(1,1);
        t(1) = tspan(1);
        firststep = true ;

        %% RK4 steps
        i = 1;
        while (t(i) < tspan (end))
            % Make sure endpoint is included
            if(h > tspan ( end) - t(i))
                h = tspan (end) - t(i);
            end

            % Full step
            fn = feval(fun,t(i),x(i,:),varargin{:})';
            [xfull,e] = RKF45Step(fun,fn,t(i),x(i,:),h,A,bhat,c,d,varargin{:});
    
            % Error estimate
            E = max(1e-10, max(abs(e)./(AbsTol+abs(xfull)*RelTol)));
    
            if(E<=1)
                % Next step
                t(i+1) = t(i)+h;
                x(i+1,:) = xfull ;
                if(firststep)
                    % New asymptotic step size
                    h = h*(epsilon/E)^(1/phat);
                    firststep = false ;
                else
                    % New PI step size
                    h = h*(epsilon/E)^kI*(Eold/E)^kp;
                end
                % Save the error for use in the PI step size controller
                Eold = E;
                i = i+1;
            else
                % New asymptotic step size
                h = h*(epsilon/E)^(1/phat);
            end
        end
    end
end

function [xnp1,e] = RKF45Step(fun,fn,tn,xn,h,A,b,c,d,varargin)
    X2 = xn + h*A(2,1)*fn;
    f2 = feval(fun,tn + c(2)*h,X2,varargin{:})';
    
    X3 = xn + h*(A(3,1)*fn + A(3,2)*f2);
    f3 = feval(fun,tn + c(3)*h,X3,varargin{:})';
    
    X4 = xn + h*(A(4,1)*fn + A(4,2)*f2 + A(4,3)*f3);
    f4 = feval(fun,tn + c(4)*h,X4,varargin{:})';
    
    X5 = xn + h*(A(5,1)*fn + A(5,2)*f2 + A(5,3)*f3 + A(5,4)*f4);
    f5 = feval(fun,tn + c(5)*h,X5,varargin{:})';
    
    X6 = xn + h*(A(6,1)*fn + A(6,2)*f2 + A(6,3)*f3 + A(6,4)*f4 +A(6,5)*f5);
    f6 = feval(fun,tn + c(6)*h,X6,varargin{:})';
    
    xnp1 = xn + h*(b(1)*fn + b(3)*f3 + b(4)*f4 + b(5)*f5 + b(6)*f6);
    %xnp1 = xn + h*(bhat(1)*fn + bhat(3)*f3 + bhat(4)*f4 + bhat(5)*f5 + bhat(6)*f6);
    e = h*(d(1)*fn + d(3)*f3 + d(4)*f4 + d(5)*f5 + d(6)*f6);
end