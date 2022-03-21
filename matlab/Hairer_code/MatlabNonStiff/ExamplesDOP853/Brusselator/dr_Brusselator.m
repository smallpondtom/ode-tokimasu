% Driver for the non-diffusive brusselator 
% Differential equations:
% dy1/dt = A + y1.^2.*y2 - (B+1)*y1;
% dy2/dt = B*y1 - y1.^2.*y2;
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

ProblemName = 'Non-diffusive brusselator';

% Function and jacobian
FcnDef = @BrusselatorFcn;

% Parameters
Bru_A = 2; 
Bru_B = 8.533;

% Initial conditions 
y0    = [1.0; 4.2665]; 

Ti    = 0;
Tf    = 20;
tspan = [Ti  Tf];

OdeCoeff{1} = Bru_A;
OdeCoeff{2} = Bru_B;

tol        = 1e-6;
RelTolDef  = tol;
AbsTolDef  = tol;

InfoStatDef    = true;
InfoDynDef     = true;

options = [];
options = rdpset(options,'RelTol',RelTolDef);
options = rdpset(options,'AbsTol',AbsTolDef);

% First, solution via odeplot

tic
[tout1,yout1] =  ode45(FcnDef,tspan,y0,options,OdeCoeff); 
T_ode45 = toc

tic
[tout,yout] = dop853(FcnDef,tspan,y0,options,OdeCoeff);    
T_dop853 = toc

% Second, soltion via normal output
tic
[tout,yout,Stats] = dop853(FcnDef,tspan,y0,options,OdeCoeff);    
T_dop853 = toc


InfoStat = Stats.Stat;
InfoDyn  = Stats.Dyn;

figure
plot(tout,real(yout(:,1)),'b')
hold on
plot(tout,real(yout(:,1)),'.b')
grid on
plot(tout,real(yout(:,2)),'r')
plot(tout,real(yout(:,2)),'.r')
title('Brusselator, u  in blue   v in red')
   
if InfoDynDef
  figure    
  subplot(3,1,1)
  plot(tout,real(yout(:,1)),'b')
  grid on
  hold on
  plot(tout,real(yout(:,1)),'.b')
  plot(tout,real(yout(:,2)),'r')
  plot(tout,real(yout(:,2)),'.r')
  title('Brusselator, u  in blue   v in red')
    
  subplot(3,1,2)
  semilogy(InfoDyn.haccept_t,InfoDyn.haccept);
  hold on
  grid on
  semilogy(InfoDyn.haccept_t,InfoDyn.haccept,'.k');
  plot(InfoDyn.hreject_t,InfoDyn.hreject,'or')
  xlim([Ti Tf])
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

