function yPrime = Ex_P297_Prime(t,y,yLag,varargin)
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
h = varargin{:};
if y(4) <= 0.1
  ksi = 1;
else
  ksi = (1-y(4))*10/9;
end
yPrime(1) = (h(1) -h(2)*y(3))*y(1);
yPrime(2) = ksi*h(3)*yLag(1,1)*yLag(3,1)-h(5)*(y(2)-1);
yPrime(3) = h(4)*(y(2)-y(3)) -h(8)*y(3)*y(1);
yPrime(4) = h(6)*y(1) -h(7)*y(4);
yPrime    = yPrime(:);

