function [Phi,tLag,tLagTab] = Ex_P287_Phi2(t,varargin)
% Differential equations with delay:
% y'   = -y(t-1), 
% y(0) = 1; 
% y(t) = Phi(t) = 0.8 + t, -1<=t<=0
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
% Phi is the value of y before the integration time
% tLag is the time Lag 
% tLagTab is the table of all time delay
% tLagTab(ny,nt): ny is the indice of y and nt is the indice of the
% corresponding delay tLag(nt)
% In this problem, ny = nt = 1.
Phi     =  0.8+t; 
tLag    = -1;
tLagTab = tLag;
