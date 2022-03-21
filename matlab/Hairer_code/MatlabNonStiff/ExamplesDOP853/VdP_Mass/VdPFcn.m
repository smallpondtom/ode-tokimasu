function dydt = VdPFcn(t,y,epsilon,MassValue)
% Differential equations:
% dy1/dt = y2
% dy2/dt = ( (1-y1^2)*y2 - y1 )/epsilon
% Parameter epsilon
% Remark : No Mass
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
dydt = [y(2,:);  ((1-y(1,:).^2).*y(2,:) - y(1,:))/epsilon] * MassValue;
