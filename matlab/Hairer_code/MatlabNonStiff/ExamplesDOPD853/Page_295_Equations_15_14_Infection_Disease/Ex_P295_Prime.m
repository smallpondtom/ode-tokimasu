function yPrime = Ex_P295_Prime(t,y,yLag)
% Differential equation with delay
% y(1)' = -y1(t)*y2(t-1) + y2(t-10) 
% y(2)' =  y1(t)*y2(t-1) - y2(t)
% y(3)' =  y2(t) - y2(t-10)
% y1(t) = Phi1(t) = 5    t<=0
% y2(t) = Phi2(t) = 0.1  t<=0
% y3(t) = Phi3(t) = 1    t<=0
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
% ------------------------------------------------------------------------
yPrime(1,1)  = -y(1)*yLag(2,1) + yLag(2,2);
yPrime(2,1)  =  y(1)*yLag(2,1) - y(2);
yPrime(3,1)  =  y(2) - yLag(2,2);

