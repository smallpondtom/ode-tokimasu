function [t,x] = EulerMod(fun,gfun,gJac,tspan,x0,AbsTol,RelTol,varargin)
% Solves g(x)' = f(t,x) using Euler 's method
%
% fun : function handle
% gfun : function handle
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
    epsilon = 0.8;
    p = 1;
    phat = p+1;
    kp = 0.4/ phat;
    kI = 0.3/ phat;
    h = 1e-3; % Initial step size if using step size controller
    hmin = 1e-6;
    hmax = 1e1;
    fixedstepsize = length(tspan) > 2;
    x = zeros(length(tspan),length(x0));
    x(1,:) = x0;
    gnp1 = feval(gfun,x0,varargin{:})';

    if( fixedstepsize )
        h = (tspan(end)-tspan(1))/(length(tspan)-1);
        t = tspan;
        for i = 1:length(t)-1
            % Euler steps
            fn = feval(fun,t(i),x(i,:),varargin{:})';
            [x(i+1,:),gnp1] = EulerModStep(fun,gfun,gJac,fn,gnp1,t(i),x(i ,:),h,AbsTol,RelTol,varargin{:});
        end
    else
        t = zeros(1,1);
        t(1) = tspan(1);
        firststep = true;
    
        % Euler steps
        i = 1;
        while (t(i) < tspan(end))
            % Make sure endpoint is included
            if(h > tspan(end) - t(i))
                h = tspan(end) - t(i);
            end

            % Full step
            fn = feval (fun,t(i),x(i ,:),varargin{:})';
            [xfull,~,Convergence,Divergence] = EulerModStep(fun,gfun,gJac,fn,gnp1,t(i),x(i,:),h,AbsTol,RelTol,varargin{:});
        
            % Double step
            [xhalf,ghalf,Con,Div] = EulerModStep(fun,gfun,gJac,fn,gnp1,t(i),x(i ,:),h/2,AbsTol,RelTol,varargin{:});
            Convergence = Convergence && Con;
            Divergence = Divergence || Div;
            fn = feval(fun,t(i)+h/2,xhalf,varargin{:})';
            [xdouble,gdouble,Con,Div,alpharatio] = EulerModStep(fun,gfun,gJac,fn,ghalf,t(i)+h/2,xhalf,h/2,AbsTol,RelTol,varargin{:});
            Convergence = Convergence && Con;
            Divergence = Divergence || Div;
        
            % Error estimate
            e = xfull - xdouble;
            E = max (1e-10 , max(abs(e)./(AbsTol+abs(xfull)*RelTol)));
        
            % fprintf ('h = %.4f, E = %.4f\n',h,E);
            % disp ( xfull )
            % disp ( xdouble )
            if(Convergence)
                if(E <= 1)
                    % disp ('Step ');
                    % Next step
                    t(i+1) = t(i)+h;
                    x(i+1,:) = xdouble;
                    gnp1 = gdouble;
        
                    if(firststep)
                        % New asymptotic step size
                        h = max(hmin,min(hmax,h*(epsilon/E)^(1/phat)));
                        firststep = false ;
                    else
                        % New PI step size
                        h = max(hmin,min(hmax,h*(epsilon/E)^kI*(Eold/E)^kp));
                    end
                    % Save the error for use in the PI step size controller
                    Eold = E;
                    i = i+1;
                else
                    % disp ('Fail ');
                    % New asymptotic step size
                    h = max(hmin,min(hmax,h*(epsilon/E)^(1/phat)));
                end
        
                if(alpharatio<1)
                    h = h* alpharatio;
                end
            elseif(Divergence)
                % disp ('Div ');
                halpha = h*alpharatio;
                h = max(0.5*h,halpha);
            else
                % disp ('Slow ');
                if(alpharatio<1)
                    halpha = h*alpharatio;
                    h = max (0.5*h,halpha);
                else
                    h = h/2;
                end
            end
        end
    end
end

function [xnp1,gnp1,Convergence,Divergence,alpharatio] = EulerModStep(fun,gfun,gJac,fn,gn,tn,xn,h,AbsTol,RelTol,varargin)
    % Full Euler step
    gnp1 = gn + h*fn;
    [xnp1,Convergence,Divergence,alpharatio] = NewtonSolve(gfun,gJac,gnp1,xn,AbsTol,RelTol,varargin{:});
end