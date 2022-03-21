function yPrime = Ex_P286_Prime(t,y,yLag)
% Differential equations with delay:
% y' = -y(t-1), 
% y(0) = 1; 
% y(t) = 1, -1<=t<=0
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
%   Denis Bichsel
%   Rue des Deurres 58
%   2000 Neuchâtel
%   Suisse
%   dbichsel@infomaniak.ch
%   End of 2015
% ---------------------------
yPrime = -yLag;
