function yPrime = Ex_P296_Prime(t,y,yLag,varargin)
% Differential equation with delay
% y(1)' = I - z*-y1(t) 
% y(2)' = z*y1(t) - y2(t)
% y(3)' =   y2(t) - y3(t)
% y(4)' =   y3(t) - 0.5*y4(t)
% z     = 1/(0.0005*y4(t-4)^3)
% y(t) = Phi() = 1, -1<=t<=0
% y0   = [63; 11; 11; 21];
% y4(t) = Phi(t) = y0, t<=0

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
I = varargin{1};
z = 1/(1+0.0005*yLag(4,1)^3);
yPrime(1) = I - z * y(1);
yPrime(2) = z* y(1) - y(2);
yPrime(3) = y(2) - y(3);
yPrime(4) = y(3) - 0.5*y(4);
yPrime = yPrime(:);