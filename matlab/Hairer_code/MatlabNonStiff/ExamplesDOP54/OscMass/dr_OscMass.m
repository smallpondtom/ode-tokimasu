% Driver for the harmonic oscillator
% Differential equations:
% y" + y = 0
% or
% du/dt = v
% dv/dt = -u
% ---------------------------
% See
%    E. Hairer S.P. Norsett G. Wanner
%    Solving Ordinary Differential Equations I
%    Nonstiff Problems
%    Springer Verlag
%    ISBN 3-540-17145-2, ISBN 0-387-17145-2
%    
% See also http://www.unige.ch/~hairer/software.html
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

clc
close all
clear all

addpath D:\RDPSolver

ProblemName = 'Oscillator';

% Function and jacobian
FcnDef  = @OscFcn;
MassDef = @OscMass;

% Initial conditions
y0 = [1;0];

Ti    = 0;
Tf    = 2*pi;
tspan = [Ti Tf];

tol        = 1e-12;
RelTolDef  = tol;
AbsTolDef  = tol;

InfoStatDef   = true;
InfoDynDef    = true;

RefineDef    = 10;

options = [];
options = rdpset(options,'RelTol',RelTolDef);
options = rdpset(options,'AbsTol',AbsTolDef);
options = rdpset(options,'Refine',RefineDef);
options = rdpset(options,'MassFcn',MassDef);
options = rdpset(options,'InitialStep',1e-8);

MassValue = 10;

tic
[tout,yout,Stats] = dop54(FcnDef,tspan,y0,options,MassValue);    
T_dop54 = toc

InfoStat = Stats.Stat;
InfoDyn  = Stats.Dyn;

figure
subplot(2,1,1)
plot(tout,real(yout(:,1)))
hold on
plot(tout,real(yout(:,1)),'.r')
grid on
title('Oscillator, first Component') 
subplot(2,1,2)
plot(tout,real(yout(:,1))- cos(tout))
grid on
title('Oscillator, Error on the first Component')

if InfoDynDef
  figure    
  subplot(3,1,1)
  plot(tout,real(yout(:,1)),'b')
  grid on
  hold on
  plot(tout,real(yout(:,1)),'.b')
  plot(tout,real(yout(:,2)),'r')
  plot(tout,real(yout(:,2)),'.r')
  title('Oscillator, u  in blue   v in red')
    
  subplot(3,1,2)
  semilogy(InfoDyn.haccept_t,InfoDyn.haccept);
  hold on
  grid on
  semilogy(InfoDyn.haccept_t,InfoDyn.haccept,'.k');
  plot(InfoDyn.hreject_t,InfoDyn.hreject,'or')
  title('Length of steps in function of time, accept = black, reject = red ')   
  
  subplot(3,1,3)
  semilogy(InfoDyn.haccept_Step,InfoDyn.haccept);
  hold on
  grid on
  semilogy(InfoDyn.haccept_Step,InfoDyn.haccept,'.k');
  plot(InfoDyn.hreject_Step,InfoDyn.hreject,'or')
  title('Length of steps in function of step number, accept = black, reject = red ')     
end 
if InfoStatDef  
  InfoStat = InfoStat
end

options = [];
options = odeset(options,'RelTol',RelTolDef);
options = odeset(options,'AbsTol',AbsTolDef);
options = odeset(options,'Refine',RefineDef);
options = odeset(options,'Mass',MassDef);
options = odeset(options,'InitialStep',1e-8);
tic
[tout,yout,Sta] = ode45(FcnDef,tspan,y0,options,MassValue);    
T_ode45 = toc

figure
subplot(2,1,1)
plot(tout,real(yout(:,1)))
hold on
plot(tout,real(yout(:,1)),'.r')
grid on
title('Oscillator, first Component') 
subplot(2,1,2)
plot(tout,real(yout(:,1))- cos(tout))
grid on
title('Oscillator, Error on the first Component')





