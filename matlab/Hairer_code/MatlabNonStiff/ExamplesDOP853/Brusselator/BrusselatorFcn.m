function yPrime = BrusselatorFcn(t,yin,OdeCoeff)
% Differential equations:
% dy1/dt = A + y1.^2.*y2 - (B+1)*y1;
% dy2/dt = B*y1 - y1.^2.*y2;
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
A         = OdeCoeff{1};
B         = OdeCoeff{2};
yPrime(1) = A + yin(1).^2.*yin(2) - (B+1)*yin(1);
yPrime(2) = B*yin(1) - yin(1).^2.*yin(2);
yPrime    = yPrime(:);