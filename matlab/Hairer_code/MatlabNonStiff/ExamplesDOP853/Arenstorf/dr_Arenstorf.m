% ---------------------------
% ARENSTORF ORBIT 
% 
% Denis Bichsel
% Rue des Deurres 58
% 2000 Neuchâtel
% Email : dbichsel@infomaniak.ch
% 2010-01-20
% ---------------------------
clear all
close all
clc

addpath ..\..\

ProblemName = 'Arenstorf';

ti = 0.0;
tf = 17.0652165601579625588917206249;
Nt = 101;
t = linspace(ti,tf,Nt);
 
t = [ti tf];

y0(1) =  0.994;
y0(2) =  0.0;
y0(3) =  0.0;
y0(4) = -2.00158510637908252240537862224;
RPAR1 = 0.012277471;
RPAR2 = 1.0-RPAR1;

warning off

optiondop = rdpset('RelTol',1e-7,'AbsTol',1e-7,'Refine',10);
optionode = odeset('RelTol',1e-7,'AbsTol',1e-7,'Refine',10);

display('dop853')
tic
[t1,y1] = dop853(@Arenstorf,t,y0,optiondop,[RPAR1,RPAR2]);
T_dop853 = toc
NbrPts(1) = length(t1);

display('ode45')
tic
[t2,y2] = ode45(@Arenstorf,t,y0,optionode,[RPAR1,RPAR2]);
T_ode45 = toc
NbrPts(2) = length(t2);

figure(1)
subplot(1,2,1); plot(y1(:,1),y1(:,2)); title(['dop853: Arenstorf, Time = ',num2str(T_dop853),'  Nbr de Pts = ',num2str(NbrPts(1))])
subplot(1,2,2); plot(y2(:,1),y2(:,2)); title(['ode45: Arenstorf, Time = ',num2str(T_ode45),'  Nbr de Pts = ',num2str(NbrPts(2))])
