%% Test ODE solver
clear; close all; clc;

addpath("exampleODE\");
addpath("unmodified\");
addpath("modified\")

% Setting
%tspan = linspace(0,30,10000);  % fixed stepsize
tspan = [0,20];  % adaptive stepsize
ic = [1,-8];
atol = 1e-10;
rtol = 1e-8;

% Simulate
[t, ys] = DOPRI54(@test05,tspan,ic,atol,rtol);  % ode solver
ye = exactSol05(t);  % exact solution

% Reshape result array thin matrix/column vector
if length(ys(1,:)) > 20  % 20 is an arbitrary number to detect fat matrix
    ys = ys';
end
if length(ye(1,:)) > 20
    ye = ye';
end

err = abs(ye-ys);  % compute error

% Plot
N = length(ys(1,:));

subplot(N,1,1)
    plot(t,ys(:,1),'.',DisplayName="sim")
    grid on; grid minor; box on; hold on;
    plot(t,ye(:,1),DisplayName="exact")
    plot(t,err(:,1),DisplayName="error")
    hold off; legend(Location="best"); ylabel("y");

if N > 1
    for i = 2:N
        subplot(N,1,i)
            plot(t,ys(:,i),'.',DisplayName="sim")
            grid on; grid minor; box on; hold on;
            plot(t,ye(:,i),DisplayName="exact")
            plot(t,err(:,i),DisplayName="error")
            hold off; legend(Location="best"); ylabel("y");
    end
end

xlabel("t")

