function yPrime = Stiff(t,yin,Alpha)
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
yPrime = - Alpha *(yin - cos(t));
