function yPrime = Ex_P293_Prime(t,y,yLag,a)
% Differential equations with delay:
% y'   = -y(t-1), 
% y(0) = 1; 
% y(t) = Phi(t) = 0, -1<=t<=0
% a is a parameter
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
%   Matlab version:
%   Denis Bichsel
%   Rue des Deurres 58
%   2000 Neuchâtel
%   Suisse
%   dbichsel@infomaniak.ch
%   End of 2015
% ---------------------------
yPrime = (a-yLag)*y;
