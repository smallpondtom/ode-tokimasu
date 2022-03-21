% Driver for the diffusive brusselator
% Differential equations:
% du/dt = g1 = A + u^2* v -(B+1)*u + alpha*d2u/dx2
% dv/dt = g2 = B*u - U^2*v + a*d2v/dx2
% Parameters  A, B and alpha
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

ProblemName = 'Diffusive brusselator';

% Function and jacobian
FcnDef  = @BrusselatorDiffusifFcn;

% Parameters
A       = 1; 
B       = 3;
alpha   = 1/50;
uBord   = 1;
vBord   = 3;

% Model geometry
Nx    = 10;
xmin  = 0;
xmax  = 1;
dx    = (xmax - xmin)/Nx;
xspan = dx:dx:1-dx;
xspan = xspan(:);       % Boundaries are in x = 0 and x = 1

% Initial conditions
u0 = 1 + sin(2*pi*xspan);
v0 = 3 * ones(size(xspan));

% To get a 5 bands matrix, the equations are interleaved
% y0 = [u0(1);v0(1);u0(2);v0(2); ... ].
N             = length(u0);
y0(1:2:2*N-1) = u0;
y0(2:2:2*N)   = v0;
y0            = y0(:);

Ti         = 0;
Tf         = 10;
NTimeInter = 100;
dt         = (Tf-Ti)/NTimeInter;
tspan      = Ti:dt:Tf;
tspan      = tspan(:);

OdeCoeff{1} = xspan;
OdeCoeff{2} = A;
OdeCoeff{3} = B;
OdeCoeff{4} = alpha;
OdeCoeff{5} = uBord;
OdeCoeff{6} = vBord;

tol       = 1e-6;
RelTolDef = tol;
AbsTolDef = tol;

InfoStatDef    = true;
InfoDynDef     = true;


options = [];
options = rdpset(options,'RelTol',RelTolDef);
options = rdpset(options,'AbsTol',AbsTolDef);

tic
for k = 1 : 1
  [tout,yout,Stats] = dop853(FcnDef,tspan,y0,options,OdeCoeff);    
end
toc

InfoStat = Stats.Stat;
InfoDyn  = Stats.Dyn;

t = tout;
y = yout;

[M,N] = size(y);
u     = real(y(:,1:2:N-1));
v     = real(y(:,2:2:N));
u     = [ones(M,1),u,ones(M,1)];       % We add the boundaries
v     = [3*ones(M,1),v,3*ones(M,1)];
figure;
subplot(1,2,1)
surf(u)
title([ProblemName,',   fonction U'])
subplot(1,2,2)
surf(v)
title([ProblemName,',   fonction V'])
 
   
if InfoDynDef
  figure  
  subplot(3,1,1)
  plot(t,u)
  grid on
  hold on
  plot(t,u,'.r')
  title('Brusselator diffusive, component no 1')
    
  subplot(3,1,2)
  semilogy(InfoDyn.haccept_t,InfoDyn.haccept,'.k');
  hold on
  grid on
  plot(InfoDyn.hreject_t,InfoDyn.hreject,'or')
  xlim([Ti Tf])
  title('Length of steps in function of time, black = accepted  red = rejected')
  
  subplot(3,1,3)
  semilogy(InfoDyn.haccept_Step,InfoDyn.haccept,'.k');
  hold on
  grid on
  plot(InfoDyn.hreject_Step,InfoDyn.hreject,'or')
  title('Length of steps in function of step number, black = accepted  red = rejected')
end
if InfoStatDef
  InfoStat = InfoStat
end


