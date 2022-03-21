% Driver for the stiff example, Vol 2 page 2 
% Differential equations:
% dy/dt = -alpha*(y-cos(x))
% Initial condition : y(0) = 0
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

clear all
close all
clc


addpath D:\RDPSolver

NomProbleme = 'Stiff page 2';

% Function and jacobian
FcnDef = @Stiff;


% Parameter
alpha = 20;

% Initial condition
y0 = 0.0;  

Ti    = 0;
Tf    = 2;
dt    = (Tf-Ti)/4;
tspan = [Ti Tf]; 

OdeCoeff    = alpha;

tol         = 1e-8;
RelTolDef   = tol;
AbsTolDef   = tol;

InfoStatDef   = true;
InfoDynDef    = true;

options = [];
options = rdpset(options,'RelTol',RelTolDef);
options = rdpset(options,'AbsTol',AbsTolDef);

tic
[tout,yout,Stats] = dop853(FcnDef,tspan,y0,options,OdeCoeff);    
T_dop853 = toc

tic
[tout1,yout1,Sta] = ode45(FcnDef,tspan,y0,options,OdeCoeff);    
T_ode45 = toc

InfoStat = Stats.Stat;
InfoDyn  = Stats.Dyn;

c        = alpha;
te       = [tspan(1):0.01:tspan(end)];
ExactSol = c^2/(1+c^2)*(- exp(-c*te) + cos(te)) + c/(1+c^2)*sin(te);

figure
plot(tout,real(yout(:,1)))
hold on
plot(tout,real(yout(:,1)),'.r')
grid on
hold on
plot(te,ExactSol,'k')
  
te       = tout;
ExactSol = c^2/(1+c^2)*(- exp(-c*te) + cos(te)) + c/(1+c^2)*sin(te);

figure
plot(tout,real(yout(:,1))- ExactSol)
grid on
 
if InfoDynDef
  figure    
  subplot(3,1,1)
  plot(tout,real(yout(:,1)),'b')
  grid on
  hold on
  plot(tout,real(yout(:,1)),'.b')
  title('Stiff Page 2')
    
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

figure
plot(tout1,real(yout1(:,1)))
hold on
plot(tout1,real(yout1(:,1)),'.r')
grid on
hold on
plot(te,ExactSol,'k')
title('ode45')

te       = tout1;
ExactSol = c^2/(1+c^2)*(- exp(-c*te) + cos(te)) + c/(1+c^2)*sin(te);
figure
plot(tout1,real(yout1(:,1))- ExactSol)
grid on
title('ode45')
