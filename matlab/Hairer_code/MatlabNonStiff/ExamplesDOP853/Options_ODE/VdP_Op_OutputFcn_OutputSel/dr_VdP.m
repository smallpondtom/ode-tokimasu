% Driver for the Van der Pol problem
% Differential equations:
% dy1/dt = y2
% dy2/dt = ( (1-y1^2)*y2 - y1 )/epsilon
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

% Options:  Gustafsson
% Input:    (VdP,tspan,y0,[],epsilon)
% Output:   odeplot or [tout, yout, InfoStat, InfoDyn]

clc
close all
clear all

addpath D:\RDPSolver

ProblemName = 'Van der Pol';

% Ode Function
FcnDef = @VdPFcn;

% Parameters
epsilon = 1e-2;

% Initial conditions 
y0 = [2; 0.0]; 

Ti    = 0;
Tf    = 2;
tspan = [Ti  Tf]';

options      = [];
OutputFcnDef = @WriteToFile1;
OutputSelDef = 1;
options      = rdpset(options,'OutputFcn',OutputFcnDef);
options      = rdpset(options,'OutputSel',OutputSelDef);
options      = rdpset(options,'Refine',10);

tic
dop853(FcnDef,tspan,y0,options,epsilon);
T_853 = toc

options      = [];
OutputFcnDef = @WriteToFile2;
OutputSelDef = [1 2];
options      = rdpset(options,'OutputFcn',OutputFcnDef);
options      = rdpset(options,'OutputSel',OutputSelDef);
options      = rdpset(options,'Refine',20);

tic
[tout,yout] = dop853(FcnDef,tspan,y0,options,epsilon);
T_dop853 = toc

clear all

FileName1 = 'VdP_Result_1.txt';
FileName2 = 'VdP_Result_2.txt';

Fid1 = fopen(FileName1,'r');
[t_y1,Count1] = fscanf(Fid1, '%g %g', [2 inf]);
t_y1 = t_y1';
fclose(Fid1);

Fid2 = fopen(FileName2,'r');
[t_y2,Count2] = fscanf(Fid1, '%g %g %g',[3 inf]);
t_y2 = t_y2';
fclose(Fid2);

figure(1)
plot(t_y1(:,1),t_y1(:,2))
grid on

figure(2)
subplot(2,1,1)
plot(t_y2(:,1),t_y2(:,2))
grid on
subplot(2,1,2)
plot(t_y2(:,1),t_y2(:,3))
grid on

