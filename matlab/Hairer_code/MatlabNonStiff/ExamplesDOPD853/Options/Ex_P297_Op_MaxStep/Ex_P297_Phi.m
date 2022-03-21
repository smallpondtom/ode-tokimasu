function [Phi,tLag,tLagTab] = Ex_P297_Phi(t,y,varargin)
% Differential equation with delay
% y(1)' = (h(1) -h(2)*y3)*y1;
% y(2)' = ksi*h(3)*yLag(1,1)*yLag(3,1)-h(5)*(y2-1);
% y(3)' = h(4)*(y2-y3) -h(8)*y3*y1;
% y(4)' = h(6)*y1 -h(7)*y4;
% y1(t) = Phi1(t) =  max(0,1e-6+t)  t <= 0 
% y3(t) = Phi3(t) = 1               t <= 0
% tau   = -0.5;
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
% Phi is the value of y before the integration time
% tLag is the time Lag 
% tLagTab is the table of all time delay
% tLagTab(ny,nt): ny is the indice of y and nt is the indice of the
% corresponding delay tLag(nt)
% In this problem, ny = nt = 1.
Phi          = [max(0,1e-6+t),0,1,0];
tau          = -0.5;
tLag         = tau;
tLagTab      = zeros(length(Phi),length(tLag));
tLagTab(1,1) = tau;
tLagTab(3,1) = tau;

