function yPrime = Ex_P287_Prime(t,y,yLag)
% Differential equation with delay
% y' = -y(t-1), y(0) = 1; 
% Three different Phi functions:
% 1) y(t) = Phi1(t) 0.8,    -1<=t<=0
% 2) y(t) = Phi1(t) 0.8+t   -1<=t<=0
% 3) y(t) = Phi1(t) 0.8+2*t -1<=t<=0
% The time delay is 1 for the three cases
% ---------------------------
% See
%    E. Hairer S.P. Norsett G. Wanner
%    Solving Ordinary Differential Equations I
%    Nonstiff Problems
%    Springer Verlag
%    ISBN 3-540-17145-2, ISBN 0-387-17145-2
%     
% See also http://www.unige.ch/~hairer/software.html 
% ------------------------------------------------------------------------
%     Matlab version:
%     Denis Bichsel
%     Rue des Deurres 58
%     2000 Neuchâtel
%     Suisse
%     dbichsel@infomaniak.ch
%     End of 2015
% ---------------------------
yPrime      = -1.4*yLag;
