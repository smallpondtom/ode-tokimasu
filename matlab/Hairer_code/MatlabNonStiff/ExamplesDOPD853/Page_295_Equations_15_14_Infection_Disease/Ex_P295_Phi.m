function [Phi,tLag,tLagTab] = Ex_P295_Phi(t,y,varargin)
% Differential equation with delay
% y(1)' = -y1(t)*y2(t-1) + y2(t-10) 
% y(2)' =  y1(t)*y2(t-1) - y2(t)
% y(3)' =  y2(t) - y2(t-10)
% y1(t) = Phi1(t) = 5    t<=0
% y2(t) = Phi2(t) = 0.1  t<=0
% y3(t) = Phi3(t) = 1    t<=0
% ---------------------------
% See
%    E. Hairer S.P. Norsett G. Wanner
%    Solving Ordinary Differential Equations I
%    Nonstiff Problems
%    Springer Verlag
%    ISBN 3-540-17145-2, ISBN 0-387-17145-2
%     
% See also http://www.unige.ch/~hairer/software.html 
%
%     Matlab version:
%     Denis Bichsel
%     Rue des Deurres 58
%     2000 Neuchâtel
%     Suisse
%     dbichsel@infomaniak.ch
%     End of 2015
% ---------------------------
% Phi contains the values of y before the integration time
% tLag contains the time Lags 
% tLagTab is the table of all time delay
% tLagTab(ny,nt): ny is the indice of y and nt is the indice of the
% corresponding delay tLag(nt)
% In this problem, ny = 3, nt = 2.
Phi          = [5;0.1;1]; 
tLag         = [-1,-10];
tLagTab      = zeros(length(Phi),length(tLag));
tLagTab(2,1) = tLag(1);
tLagTab(2,2) = tLag(2);

