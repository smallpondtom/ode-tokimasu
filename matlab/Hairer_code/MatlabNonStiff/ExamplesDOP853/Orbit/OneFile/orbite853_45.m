function orbite853_45
%ORBITODE  Restricted three-body problem
clear all
close all
clc
addpath D:\RDPSolver

mu = 1 / 82.45;
mustar = 1 - mu;
y0 = [1.2; 0; 0; -1.04935750983031990726];
tspan = [0 7];
   
options = rdpset('RelTol',1e-5,'AbsTol',1e-4,'Refine',10,  ...
                 'Events',@events);         
tic
[t,y,te,ye,ie] = dop853(@f,tspan,y0,options);
T_dop853 = toc

plot(y(:,1),y(:,2),ye(:,1),ye(:,2),'o');
title ('Restricted three body problem')
ylabel ('y(t)')
xlabel ('x(t)')


options = odeset('RelTol',1e-5,'AbsTol',1e-4,...
                 'Events',@events);
tic
[t1,y1,te1,ye1,ie1] = ode45(@f,tspan,y0,options);
T_ode45 = toc

figure
plot(y1(:,1),y1(:,2),ye1(:,1),ye1(:,2),'o');
title ('Restricted three body problem')
ylabel ('y(t)')
xlabel ('x(t)')


% --------------------------------------------------------------
function dydt = f(t,y)
r13 = ((y(1) + mu)^2 + y(2)^2) ^ 1.5;
r23 = ((y(1) - mustar)^2 + y(2)^2) ^ 1.5;
dydt = [ y(3)
         y(4)
         2*y(4) + y(1) - mustar*((y(1)+mu)/r13) - ...
                         mu*((y(1)-mustar)/r23)
        -2*y(3) + y(2) - mustar*(y(2)/r13) - mu*(y(2)/r23) ];
end   % End nested function f
% --------------------------------------------------------------
function [value,isterminal,direction] = events(t,y)
% Locate the time when the object returns closest to the
% initial point y0 and starts to move away; stop integration.
% Also locate the time when the object is farthest from the 
% initial point y0 and starts to move closer.
% 
% The current distance of the body is
% 
%   DSQ = (y(1)-y0(1))^2 + (y(2)-y0(2))^2 
%       = <y(1:2)-y0(1:2),y(1:2)-y0(1:2)>
% 
% A local minimum of DSQ occurs when d/dt DSQ crosses zero 
% heading in the positive direction. Compute d(DSQ)/dt as
% 
%  d(DSQ)/dt = 2*(y(1:2)-y0(1:2))'*dy(1:2)/dt = ... 
%                 2*(y(1:2)-y0(1:2))'*y(3:4)
% 
dDSQdt = 2 * ((y(1:2)-y0(1:2))' * y(3:4));
value = [dDSQdt; dDSQdt];
isterminal = [1; 0];            % Stop at local minimum
direction  = [1; -1];           % [local minimum, local maximum]
end   % End nested function events
end 