%% Test ODE solver
clear; close all; clc;

% Get user input to decide which test to run
which_test = input("Which test would you like choose a number from 1 to 13 -> ");
setup_test;

% Get user input to decide which method to use
flag = false;
while ~flag
    which_method = input("Which method would you like? (1) Euler, (2) RK4, (3) RKF45, or (4) DOPRI54 -> ");
    if which_method == 1
        flag = length(ic) == 1;
        if ~flag
            fprintf("Euler can only be used for 1st-order examples. \nPlease select " + ...
                "a different method or start over and select a different test.\n");
        end
    else
        flag = true;
    end
end
setup_solver;

tic;
% Simulate
[t, ys] = solver(fs,tspan,ic,1e-10,1e-8);  % ode solver
ye = fe(t);  % exact solution 

% Reshape result array thin matrix/column vector
if length(ys(1,:)) > 20  % 20 is an arbitrary number to detect fat matrix
    ys = ys';
end
if length(ye(1,:)) > 20
    ye = ye';
end
err = abs(ye-ys);  % compute error
toc;

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
