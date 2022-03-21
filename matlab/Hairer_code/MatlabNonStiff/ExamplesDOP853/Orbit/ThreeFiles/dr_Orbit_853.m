% dr_Orbit_54
%ORBITODE  Restricted three-body problem

clc
clear all
close all

addpath D:\RDPSolver

mu = 1 / 82.45;
mustar = 1 - mu;
y0 = [1.2; 0; 0; -1.04935750983031990726];
tspan = [0 7];

options = rdpset('RelTol',1e-5,'AbsTol',1e-4,...
                 'Events',@OrbiteEvents);
tic
[t,y,te,ye,ie] = dop853(@Orbite,tspan,y0,options,mu,mustar,y0);
T_dop853 = toc

plot(y(:,1),y(:,2),ye(:,1),ye(:,2),'o');
title ('Restricted three body problem')
ylabel ('y(t)')
xlabel ('x(t)')

options = odeset('RelTol',1e-5,'AbsTol',1e-4,...
                 'Events',@OrbiteEvents);
tic
[t1,y1,te1,ye1,ie1] = ode45(@Orbite,tspan,y0,options,mu,mustar,y0);
T_ode45 = toc

plot(y1(:,1),y1(:,2),ye1(:,1),ye1(:,2),'o');
title ('Restricted three body problem')
ylabel ('y(t)')
xlabel ('x(t)')

% --------------------------------------------------------------
