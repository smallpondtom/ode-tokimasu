% Driver for the Van der Pol problem
% Differential equations:
% dy1/dt = y2
% dy2/dt = ( (1-y1^2)*y2 - y1 )/epsilon
% Parameter epsilon
% Remark : No Mass
% ---------------------------
% See
%    E. Hairer S.P. Norsett G. Wanner
%    Solving Ordinary Differential Equations I
%    Nonstiff Problems
%    Springer Verlag
%    ISBN 3-540-17145-2, ISBN 0-387-17145-2
%    
% See also  http://www.unige.ch/~hairer/software.html
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

ProblemName = 'Van der Pol';

% Function, jacobian, mass
FcnDef      = @VdPFcn;

% Parameter
epsilon   = 1e-2;
MassValue = 1;

% Initial conditions 
y0 = [2; 0.0]; 

Ti    = 0;
Tf    = 2;
tspan = [Ti Tf];

tol        = 1e-4;
RelTolDef  = tol;
AbsTolDef  = tol;

InfoStatDef    = true;
InfoDynDef     = true;


options = [];
options = rdpset(options,'Beta', 0.15);
tic
dop54(FcnDef,tspan,y0,options,epsilon,MassValue);    
toc

tic
[tout,yout,Stats] = dop54(FcnDef,tspan,y0,options,epsilon,MassValue);    
T_dop54 = toc
whos

options = [];
options = dopset(options,'RelTol',RelTolDef);
options = dopset(options,'AbsTol',AbsTolDef);
options = dopset(options,'OutputSel',1);
options = dopset(options,'Mass',@VdPMass);
options = rdpset(options,'Refine',1); 

  
  InfoStat = Stats.Stat;
  InfoDyn  = Stats.Dyn;


if InfoDynDef
  figure  
  subplot(3,1,1)
  plot(tout,real(yout(:,1)))
  grid on
  hold on
  plot(tout,real(yout(:,1)),'.k')
  title('Van der Pol, component no 1')
    
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

