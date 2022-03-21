function yPrime  = Wille_and_Backer_Prime(t,y,yLag,varargin)
% Differential equation with delay
% y(1)' = y1(t-1)  
% y(2)' = y1(t-1) + y_2(t-0.2)
% y(3)' = y2(t)
% y(0)  = [1,1,1]
% y(t) = Phi(t) = [1;1;1]  t <= 0
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
yPrime(1) = yLag(1,1);
yPrime(2) = yLag(1,1) + yLag(2,2);
yPrime(3) = y(2);
yPrime    = yPrime(:);