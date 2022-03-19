function [xnp1,Convergence,Divergence,alpharatio] = NewtonSolve(gfun,gJac,gnp1,x0,AbsTol,RelTol,varargin)
    % Various parameters
    SlowConvergence = false;
    Divergence = false;
    itermax = 1e1;
    epsilon = 0.8;
    tau = 0.1*epsilon;
    alpha = 0;
    alpharef = 0.4;
    
    %% Obtain Solution
    % Start guess
    xnp1 = x0;
    g = feval(gfun,xnp1,varargin{:})';
    dg = feval(gJac,xnp1,varargin{:});
    [L,U,pivot] = lu(dg,'vector');

    % Residual
    R = g-gnp1;
    rNewton = norm(R./(AbsTol+abs(gnp1)*RelTol),inf);
    rNewtonOld = rNewton;
    iter = 0;
    Convergence = rNewton < tau;

    % Newton steps
    while (~Convergence && ~SlowConvergence && ~Divergence)
        R = R';
        dxn = U\(L\(R(pivot,1)));
        xnp1 = xnp1 - dxn';
        g = feval(gfun,xnp1,varargin{:})';
        R = g - gnp1;
        rNewton = norm(R./(AbsTol+abs(gnp1)*RelTol),'inf');
        alpha = max(alpha,rNewton/rNewtonOld);
        Convergence = rNewton < tau;
        SlowConvergence = iter > itermax ;
        Divergence = alpha > 1;
        rNewtonOld = rNewton;
        iter = iter+1;
    end

    alpharatio = alpharef / alpha ;
end