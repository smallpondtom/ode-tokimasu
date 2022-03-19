function [t,x] = RK4(fun,tspan,x0,AbsTol,RelTol,varargin)
% Solves x' = f(t,x) using the Classic Runge - Kutta method
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
    % Butcher Tableau
    A = [ 0,   0, 0, 0;
        1/2,   0, 0, 0;
          0, 1/2, 0, 0;
          0,   0, 1, 0];
    b = [1/6,1/3,1/3,1/6];
    c = [0,1/2,1/2,1];
    
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
    
    %% INTEGRATION
    if(fixedstepsize)
        h = (tspan(end) - tspan(1))/(length(tspan) -1);
        t = tspan;
    
        for i = 1:length(t)-1
            % RK4 steps
            fn = feval(fun,t(i),x(i,:),varargin{:})';
            x(i+1,:) = RK4Step(fun,fn,t(i),x(i,:),h,A,b,c,varargin{:});
        end
    else
        t = zeros (1 ,1);
        t(1) = tspan (1);
        firststep = true ;
    
        % RK4 steps
        i = 1;
        while(t(i)<tspan(end))
            % Make sure endpoint is included
            if(h > tspan(end) - t(i))
                h = tspan(end) - t(i);
            end
            % Full step
            fn = feval(fun,t(i),x(i,:),varargin{:})';
            xfull = RK4Step(fun,fn,t(i),x(i,:),h,A,b,c,varargin{:});
            
            % Double step
            xhalf = RK4Step(fun,fn,t(i),x(i,:),h/2,A,b,c,varargin{:});
            fn = feval(fun,t(i)+h/2,xhalf,varargin{:})';
            xdouble = RK4Step(fun,fn,t(i)+h/2,xhalf,h/2,A,b,c,varargin{:});
    
            % Error estimate
            e = xfull - xdouble ;
            E = max (1e-10,max(abs(e)./(AbsTol + abs(xfull)*RelTol)));
    
            if(E<=1)
                % Next step
                t(i+1) = t(i)+h;
                x(i+1,:) = xdouble;
    
                if(firststep)
                    % New asymptotic step size
                    h = h*(epsilon/E)^(1/phat);
                    firststep = false;
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

function xnp1 = RK4Step(fun,fn,tn,xn,h,A,b,c,varargin)
    X2 = xn + A(2,1)*h*fn;
    f2 = feval(fun,tn+c(2)*h,X2,varargin{:})';
    X3 = xn + A(3,2)*h*f2;
    f3 = feval(fun,tn + c (3)*h,X3,varargin{:})';
    X4 = xn + A(4,3)*h*f3;
    f4 = feval(fun,tn + h,X4,varargin{:})';
    xnp1 = xn + h*(b(1)*fn + b(2)*f2 + b(3)*f3 + b(4)*f4);
end