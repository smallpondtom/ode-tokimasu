% Driver for the Ballode 1 problem
% BALLODE Run a demo of a bouncing ball.
% Differential equations:
% dy1/dt = y2
% dy2/dt = -g = - 9.81
% Locate the time when height passes through zero in a decreasing direction
% and stop integration.
% ---------------------------
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
%   Denis Bichsel
%   Rue des Deurres 58
%   2000 Neuchâtel
%   Suisse
%   dbichsel@infomaniak.ch
%   End of 2015
% ---------------------------
clear all
close all
clc

addpath D:\RDPSolver

ProblemName = 'Ballode 1';

tstart = 0;
tfinal = 30;
tspan  = [tstart  tfinal];
y0     = [0,15];
    
OdeFcn    = @BallodeFcn;
EventsFcn = @BallodeEventsFcn;

optionsrdp = rdpset('Events',EventsFcn,'Refine',10, ...
                    'OutputSel',[1,2],'InitialStep',0.1,'MaxStep',0.2);
optionsode = odeset('Events',EventsFcn,'Refine',10, ...
                    'OutputSel',[1,2],'InitialStep',0.1,'MaxStep',0.2);               
               
tout  = []; 
yout  = []; 
teout = [];
yeout = [];
ieout = [];
ie    = [];
te    = [];
dt    = 0.00001;
dy0   = dt*y0(2);

tic
for i = 1:20 

  [t,y,te,ye,ie] = dop853(OdeFcn,tspan,y0,optionsrdp);     

  % Accumulate output
  nt    = length(t);
  tout  = [tout; t(1:nt-1)];
  yout  = [yout; y(1:nt-1,:)];    
  teout = [teout; te];          % Events at tstart are never reported.
  yeout = [yeout; ye];
  ieout = [ieout; ie];
  
  % Set the new initial conditions, with .9 attenuation.  
  if ~isempty(ie)
    switch ie
      case 1           
        y0(2) = -0.9*ye(2); 
        y0(1) = dt*y0(2);
        tspan = [te,tfinal];
        tout  = [tout;te];
        te    = te + dt;
        yout  = [yout;ye];
        tspan = [te,tfinal];
      case 2
        te    = te + dt;
        y0    = ye;
        y0(1) = y0(1) + dt*ye(2);
        tspan = [te,tfinal];
    end
  end
  
  if isempty(te)
    break
  end  
end

T_dop853 = toc

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


tstart = 0;
tfinal = 30;
tspan  = [tstart  tfinal];
y0     = [0,15];

tout  = []; 
yout  = []; 
teout = [];
yeout = [];
ieout = [];
ie    = [];
te    = [];
dt    = 0.00001;
dy0   = dt*y0(2);

tic
for i = 1:20 

  [t,y,te,ye,ie] = ode45(OdeFcn,tspan,y0,optionsode);     

  % Accumulate output
  nt    = length(t);
  tout  = [tout; t(1:nt-1)];
  yout  = [yout; y(1:nt-1,:)];    
  teout = [teout; te];          % Events at tstart are never reported.
  yeout = [yeout; ye];
  ieout = [ieout; ie];
  
  % Set the new initial conditions, with .9 attenuation.  
  if ~isempty(ie)
    switch ie
      case 1           
        y0(2) = -0.9*ye(2); 
        y0(1) = dt*y0(2);
        tspan = [te,tfinal];
        tout  = [tout;te];
        te    = te + dt;
        yout  = [yout;ye];
        tspan = [te,tfinal];
      case 2
        te    = te + dt;
        y0    = ye;
        y0(1) = y0(1) + dt*ye(2);
        tspan = [te,tfinal];
    end
  end
  
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

