function [t,x] = DOPRI54(fun,tspan,x0,AbsTol,RelTol,varargin)
% Solves x' = f(t,x) using the Dormand - Prince method
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
    % DOPRI54 Butcher Tableau
    A = [0,0,0,0,0,0,0;
         1/5,0,0,0,0,0,0;
         3/40,9/40,0,0,0,0,0;
         44/45,-56/15,32/9,0,0,0,0;
         19372/6561,-25360/2187,64448/6561,-212/729,0,0,0;
         9017/3168,-355/33,46732/5247,49/176,-5103/18656,0,0;
         35/384,0,500/1113,125/192,-2187/6784,11/84,0];
    b = [5179/57600,0,7571/16695,393/640,-92097/339200,187/2100,1/40];
    bhat = [35/384,0,500/1113,125/192,-2187/6784,11/84,0];
    c = [0,1/5,3/10,4/5,8/9,1,1];
    d = b-bhat ;

    % Various parameters
    epsilon = 0.8;
    p = 4;  % order
    kp = 0.4/(p+1);
    kI = 0.3/(p+1);
    h = 1e-3;  % Initial step size if using step size controller
    fixedstepsize = length(tspan) > 2;
    x = zeros(length(tspan),length(x0));
    x(1,:) = x0;
    fn = feval(fun,tspan(1),x(1,:),varargin{:})';
    
    %% INTEGRATE
    if(fixedstepsize)
        h = (tspan(end) - tspan(1))/(length(tspan) -1);
        t = tspan;

        for i = 1:length(t)-1
            % DOPRI54 steps
            [x(i+1,:),~,fn] = DOPRI54Step(fun,fn,t(i),x(i,:),h,A,c,d,varargin{:});
        end
    else
        t = zeros(1,1);
        t(1) = tspan(1);
        firststep = true;

        % DOPRI54 steps
        i = 1;
        while (t(i) < tspan(end))
            % Make sure endpoint is included
            if(h > tspan(end) - t(i))
                h = tspan(end) - t(i);
            end

            % Full step
            [xfull,e,fnp1] = DOPRI54Step(fun,fn,t(i),x(i,:),h,A,c,d,varargin{:});
            
            % Error estimate
            E = max(1e-10,max(abs(e)./(AbsTol+abs(xfull)*RelTol)));
            if(E<=1)
                % Next step
                t(i+1) = t(i)+h;
                x(i+1,:) = xfull;
                fn = fnp1;
                if(firststep)
                    % New asymptotic step size
                    h = h*(epsilon/E)^(1/(p+1));
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
            h = h*(epsilon/E)^(1/(p+1));
            end
        end
    end
end

function [xnp1 ,e, fnp1 ] = DOPRI54Step (fun ,fn ,tn ,xn ,h,A,c,d,varargin )
    X2 = xn + A(2,1)*h*fn;
    f2 = feval (fun,tn + c(2)*h,X2,varargin{:})';
    
    X3 = xn + h*(A(3,1)*fn + A(3,2)*f2);
    f3 = feval(fun,tn + c(3)*h,X3,varargin{:})';
    
    X4 = xn + h*(A(4,1)*fn + A(4,2)*f2 + A(4,3)*f3);
    f4 = feval(fun,tn + c(4)*h,X4,varargin{:})';
    
    X5 = xn + h*(A(5,1)*fn + A(5,2)*f2 + A(5,3)*f3 + A(5,4)*f4);
    f5 = feval(fun,tn + c(5)*h,X5,varargin{:})';
    
    X6 = xn + h*(A(6,1)*fn + A(6,2)*f2 + A(6,3)*f3 + A(6,4)*f4 + A(6,5)*f5);
    f6 = feval(fun,tn + c(6)*h,X6, varargin{:})';
    
    X7 = xn + h*(A(7,1)*fn + A(7,3)*f3 + A(7,4)*f4 + A(7,5)*f5 + A(7,6)*f6);
    f7 = feval(fun,tn + c(7)*h,X7,varargin{:})';
    
    xnp1 = X7;
    e = h*(d(1)*fn + d(3)*f3 + d(4)*f4 + d(5)*f5 + d(6)*f6 + d(7)*f7);
    fnp1 = f7;
end