function [t,x] = ESDIRK23(fun,Jac,tspan,x0,AbsTol,RelTol,varargin)
% Solves x' = f(t,x) using the ESDIRK23 method
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
    % ESDIRK23 Parameters
    gamma = 1-1/sqrt(2);
    a31 = (1-gamma)/2;
    c = [0;2*gamma;1];
    b = [a31;a31;gamma];
    bhat = [(6*gamma-1)/(12*gamma); 1/(12*gamma*(1-2*gamma)); (1-3*gamma)/(3*(1-2*gamma))];
    d = b-bhat;
    
    % Various parameters
    epsilon = 0.8;
    p = 2;  % order
    phat = p+1;
    kp = 1/phat;
    kI = 1/phat;
    h = 1e-3;  % Initial step size if using step size controller
    hmin = 1e-6;
    hmax = 1e1;
    fixedstepsize = length(tspan) > 2;
    x = zeros(length(tspan),length(x0));
    x(1,:) = x0;
    fn = feval(fun,tspan(1),x(1,:),varargin{:})';
    if(fixedstepsize)
        h = (tspan(end) - tspan(1))/(length(tspan) -1);
        t = tspan;
        for i = 1:length(t)-1
            % ESDIRK23 steps
            [x(i+1,:),~,fn] = ESDIRK23Step(fun,Jac,fn,t(i),x(i,:),h,AbsTol, ...
                RelTol,b,c,d,gamma,varargin{:});
        end
    else
        t = zeros(1,1);
        t(1) = tspan(1);
        firststep = true;

        % Main ESDIRK Integrator
        i = 1;
        while (t(i) < tspan(end))
            % Make sure endpoint is included
            if(h > tspan(end) - t(i))
                h = tspan(end) - t(i);
            end
            % Full step
            [xfull,e,fnp1,Convergence,Divergence,alpharatio] = ESDIRK23Step(fun,Jac,fn,t(i),x(i,:), ...
                h,AbsTol,RelTol,b,c,d,gamma,varargin{:});

            % Error estimate
            E = max (1e-10,max(abs(e)./(AbsTol + abs(xfull)*RelTol)));
            % Step size control
            if(Convergence)
                if(E<=1)
                    % Next step
                    t(i+1) = t(i)+h;
                    x(i+1,:) = xfull ;
                    fn = fnp1;
                        if(firststep)
                            % New asymptotic step size
                            h = h*(epsilon/E)^(1/phat);
                            firststep = false;
                        else
                            % New PI step size
                            hrat = h/h_old;
                            Erat1 = (epsilon/E)^kI;
                            Erat2 = (Eold/E)^kp;
                            h = max(hmin,min(hmax,h*hrat*Erat1*Erat2));
                        end
                    % Save the error for use in PI step size controller
                    Eold = E;
                    h_old = h;
                    i = i+1;
                else
                    % New asymptotic step size
                    h = max(hmin,min(hmax,h*(epsilon/E)^(1/phat)));
                end

                if(alpharatio < 1)
                    h = h*alpharatio;
                end
            elseif (Divergence)
            halpha = h*alpharatio;
            h = max(0.5*h,halpha);
            else
                if(alpharatio < 1)
                    halpha = h*alpharatio;
                    h = max(0.5*h,halpha);
                else
                    h = h/2;
                end
            end
        end
    end
end

function [xnp1,e,fnp1,Convergence,Divergence,alpharatio] = ESDIRK23Step(fun,Jac,fn,tn,xn,h,AbsTol,RelTol,b,c,d,gamma,varargin)
    % Various parameters
    SlowConvergence = false;
    Divergence = false;
    a21 = gamma;
    itermax = 1e1;
    epsilon = 0.8;
    tau = 0.1*epsilon;
    alpha = 0;
    alpharef = 0.4;
    I = eye(length(xn));

    % Jacobian Update
    J = feval(Jac,tn,xn,varargin{:});
    dRdx = I - h*gamma*J;
    [L,U,pivot] = lu(dRdx,'vector');
    
    % Stage 2 of the ESDIRK23 Method
    psi2 = xn + h*a21*fn;
    
    % Initial guess for the state
    T2 = tn + c(2)*h;
    X2 = xn + c(2)*h*fn;
    
    f2 = feval(fun,T2,X2,varargin{:})';
    R2 = (X2 - h*gamma*f2 - psi2 )';
    rNewton = norm (R2'./(AbsTol + abs(X2).*RelTol),inf);
    rNewtonOld = rNewton;
    iter = 0;
    Convergence = false;
    while (~Convergence && ~SlowConvergence && ~Divergence)
        dX = U\(L\(-R2(pivot,1)));
        X2 = X2 + dX';
        f2 = feval(fun,T2,X2,varargin{:})';
        R2 = (X2 - h*gamma*f2 - psi2)';
        rNewton = norm(R2'./(AbsTol + abs(X2).*RelTol),inf);
        alpha = max(alpha,rNewton/rNewtonOld);
        Convergence = rNewton < tau;
        SlowConvergence = iter > itermax;
        Divergence = alpha > 1;
        rNewtonOld = rNewton;
        iter = iter+1;
    end
    
    % Stage 3
    psi3 = xn + h*(b(1)*fn+b(2)*f2);
    
    % Initial guess for the state
    T3 = tn + c(3)*h;
    X3 = xn + c(3)*h*fn;
    
    f3 = feval(fun,T3,X3,varargin{:})';
    R3 = (X3 - h*gamma*f3-psi3)';
    rNewton = norm(R3'./(AbsTol + abs(X3).*RelTol),inf);
    rNewtonOld = rNewton;
    iter = 0;
    Convergence = false;
    while (~Convergence && ~SlowConvergence && ~Divergence)
        dX = U\(L\(-R3(pivot,1)));
        X3 = X3 + dX';
        f3 = feval(fun,T3,X3,varargin{:})';
        R3 = (X3 - h*gamma*f3 - psi3 )';
        rNewton = norm(R3'./(AbsTol + abs(X3).*RelTol),inf);
        alpha = max(alpha, rNewton / rNewtonOld);
        Convergence = rNewton < tau;
        SlowConvergence = iter > itermax ;
        Divergence = alpha > 1;
        rNewtonOld = rNewton;
        iter = iter+1;
    end
    
    xnp1 = X3;
    fnp1 = f3;
    e = h*(d(1)*fn + d(2)*f2 + d(3)*f3);
    alpharatio = alpharef / alpha ;
end