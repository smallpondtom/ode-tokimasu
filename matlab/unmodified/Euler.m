function [t,x] = Euler(fun,tspan,x0,AbsTol,RelTol,varargin)
% Solves x' = f(t,x) using Euler 's method
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
    % Various parameters
    epsilon = 0.8;
    p = 1;
    phat = p+1;
    kp = 0.4/phat;
    kI = 0.3/phat;
    h = 1e-3;  % Initial step size if using step size controller
    fixedstepsize = length(tspan) > 2;
    x = zeros(length(tspan),length(x0));
    x(1,:) = x0;
    
    %% LOOP
    % Fixed step size
    if( fixedstepsize )
        h = (tspan(end)-tspan(1))/(length(tspan)-1);
        t = tspan;
        for i = 1:length(t) -1
            % Euler steps
            fn = feval(fun,t(i),x(i ,:),varargin{:})';
            x(i+1,:) = EulerStep(fun,fn,t(i),x(i,:),h,varargin{:});
        end
    else
        % Step size control
        t = zeros (1,1);
        t(1) = tspan(1);
        firststep = true ;
        i = 1;
        while (t(i) < tspan(end))
            % Make sure endpoint is included
            if(h > tspan(end) - t(i))
                h = tspan(end) - t(i);
            end
            % Full step
            fn = feval(fun,t(i),x(i,:) ,varargin{:})';
            xfull = EulerStep(fun,fn,t(i),x(i,:),h,varargin{:});
            % Double step
            xhalf = EulerStep(fun,fn,t(i),x(i,:),h/2,varargin{:});
            fn = feval(fun,t(i)+h/2,xhalf,varargin{:})';
            xdouble = EulerStep(fun,fn,t(i)+h/2,xhalf,h/2,varargin{:});
            % Error estimate
            e = xfull - xdouble ;
            E = max(1e-10,max(abs(e)./(AbsTol+abs(xfull)*RelTol)));
            % Accept or dismiss step
            if(E<=1)
                % Next step
                t(i+1) = t(i)+h;
                x(i+1,:) = xdouble;
                if(firststep)
                    % If it is the first step use asymptotic step size controller
                    h = h*(epsilon/E)^(1/phat);
                    firststep = false;
                else
                    % New PI step size
                    h = h*(epsilon/E)^kI *(Eold/E)^kp;
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

function xnp1 = EulerStep(fun,fn,tn,xn,h,varargin)
    xnp1 = xn + h*fn;
end