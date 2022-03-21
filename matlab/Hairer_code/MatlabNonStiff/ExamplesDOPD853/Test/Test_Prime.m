function yPrime = Test_Prime(t,y,yLag)
% Differential equation with delay
% y' = 3*y(t) - 2*y(t-1) - 3*t^2 - 4*t + 7
% y(0) = 1
% y(t) = Phi(t) = 3t^2 - 2t + 1 
% The solution is analytic everywhere and the exact solution is:
% y =  3t^2 - 2t + 1
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
yPrime = 3*y - 2*yLag - 3*t^2 - 4*t + 7;
