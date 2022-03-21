% Can I use the events function more than once?
% Let's take the ballode demo for example. 
% The events function is used to detect when the ball hits the ground. 
% What if I also wanted to add another events function to detect 
% something else? Maybe I also want to know when the ball is at a height 
% of 0.5. This would mean I would have two different events functions 
% for the same ode solver in the for loop.

function Ballode1
%BALLODE Run a demo of a bouncing ball.

clear all
close all
clc

addpath D:\DOPSolver

tstart = 0;
tfinal = 30;
tspan  = [tstart,tfinal];
y0     = [-0.00000,15];
               
options = rdpset('Events',@events,...
                 'OutputSel',[1,2],'InitialStep',0.1,'MaxStep',0.2);               
               
tout  = []; 
yout  = []; 
teout = [];
yeout = [];
ieout = [];
ie = [];
te = [];
dt = 0.0001;
tic
for i = 1:20 
  [t,y,te,ye,ie] = dop853(@f,tspan,y0,options);     
  % Accumulate output.
  nt    = length(t);
  tout  = [tout; t(1:nt-1)];
  yout  = [yout; y(1:nt-1,:)];
  teout = [teout; te];          % Events at tstart are never reported.
  yeout = [yeout; ye];
  ieout = [ieout; ie];
  
  if ~isempty(ie)
    switch ie(end)
      case 1    
        y0(1) = dt;
        y0(2) = -0.9*ye(2);        
        te    = te+dt;
        Me    = length(te);
        tspan = [te(Me),tfinal];
        tout  = [tout;te(Me)];        
        yout  = [yout;ye(Me,:)];
        tspan = [te(Me),tfinal];
      case 2
        te = te + dt;
        y0 = ye;
        y0(1) = y0(1) + dt*ye(2);
        tspan = [te,tfinal];
    end
  end
  % Set the new initial conditions, with .9 attenuation.    
  if isempty(te)
    break
  end  
end
T_dop853 = toc

plot(teout,yeout(:,1),'ro')
xlabel('time');
ylabel('height');
title('Ball trajectory and the events');
hold off
odeplot([],[],'done');

figure 
plot(tout,yout(:,1),'b')
hold on
plot(teout,yeout(:,1),'ro')


clear all
tstart = 0;
tfinal = 30;
tspan  = [tstart,tfinal];
y0     = [-0.00000,15];
               
options = odeset('Events',@events,...
                 'OutputSel',[1,2],'InitialStep',0.1,'MaxStep',0.2);               
               
tout  = []; 
yout  = []; 
teout = [];
yeout = [];
ieout = [];
ie = [];
te = [];
dt = 0.0001;
tic
for i = 1:20 
  [t,y,te,ye,ie] = ode45(@f,tspan,y0,options);     
  % Accumulate output.
  nt    = length(t);
  tout  = [tout; t(1:nt-1)];
  yout  = [yout; y(1:nt-1,:)];
  teout = [teout; te];          % Events at tstart are never reported.
  yeout = [yeout; ye];
  ieout = [ieout; ie];
  
  if ~isempty(ie)
    switch ie(end)
      case 1    
        y0(1) = dt;
        y0(2) = -0.9*ye(2);        
        te    = te+dt;
        Me    = length(te);
        tspan = [te(Me),tfinal];
        tout  = [tout;te(Me)];        
        yout  = [yout;ye(Me,:)];
        tspan = [te(Me),tfinal];
      case 2
        te = te + dt;
        y0 = ye;
        y0(1) = y0(1) + dt*ye(2);
        tspan = [te,tfinal];
    end
  end
  % Set the new initial conditions, with .9 attenuation.    
  if isempty(te)
    break
  end  
end
T_ode45 = toc

figure
plot(teout,yeout(:,1),'ro')
xlabel('time');
ylabel('height');
title('Ball trajectory and the events');
hold off
odeplot([],[],'done');

figure 
plot(tout,yout(:,1),'b')
hold on
plot(teout,yeout(:,1),'ro')
% --------------------------------------------------------------
function dydt = f(t,y)
dydt = [y(2); -9.8];
% --------------------------------------------------------------
% function [value,isterminal,direction] = events(t,y)
% % Locate the time when height passes through zero in a
% % decreasing direction and stop integration.
% value = y(1); % Detect height = 0
% isterminal = 1; % Stop the integration
% direction = -1; % Negative direction only 

% -------------------------------------------------------------------
% 
% "ill will" <schoolsofthought@gmail.com> wrote in message
% news:gnk9u7$65u$1@fred.mathworks.com...
% > Can I use the events function more than once?
% > Let's take the ballode demo for example. The events function is used to
% > detect when the ball hits the ground. What if I also wanted to add another
% > events function to detect something else? Maybe I also want to know when
% > the ball is at a height of 0.5. This would mean I would have two different
% > events functions for the same ode solver in the for loop.
% 
% You can't do this exactly the way you described, but you can do what you
% want:
% 
% http://www.mathworks.com/access/helpdesk/help/techdoc/ref/odeset.html#f92-1017470
% 
% "value, isterminal, and direction are vectors for which the ith element
% corresponds to the ith event function: "
% 
% So you can have your events function check for multiple events and update
% the appropriate element of the output vectors. Your events function would
% become:
% 
% 
function [value,isterminal,direction] = events(t,y)
% Locate the time when height passes through zero in a
% decreasing direction and stop integration.
value      = [y(1)]; % Detect height = 0 or 0.5
isterminal = [ 1];  % Stop the integration
direction  = [-1]; % 

% value      = [y(1)]; % Detect height = 0 or 0.5
% isterminal = [1]; % Stop the integration
% direction  = [-1]; % 
% 
% 
% Of course, in this example, assuming the ball starts above height 0.5, only
% the second event should trigger.
% 
% -- 
% Steve Lord
% slord@mathworks.com




