function Mass = VdPMass(t,y,epsilon)
% Differential equations:
% 1 * dy1/dt      = y2
% epsilon *dy2/dt = ( (1-y1^2)*y2 - y1 )
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
%   2000 Neuch�tel
%   Suisse
%   dbichsel@infomaniak.ch
%   End of 2015
% ---------------------------
Mass = sparse([1, 0 ; 0, epsilon]);
