function [Phi,tLag,tLagTab] = Ex_P296_Phi(t,y,varargin)
% Differential equation with delay
% y(1)' = I - z*-y1(t) 
% y(2)' = z*y1(t) - y2(t)
% y(3)' =   y2(t) - y3(t)
% y(4)' =   y3(t) - 0.5*y4(t)
% z     = 1/(0.0005*y4(t-4)^3)
% y(t) = Phi() = 1, -1<=t<=0
% y0   = [63; 11; 11; 21];
% y4(t) = Phi(t) = y0, t<=0
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
varargin{:}
Phi4         = varargin{2};
Phi          = [0,0,0,21];
tLag         = -4;
tLagTab      = zeros(length(Phi),length(tLag));
tLagTab(4,1) = tLag;

