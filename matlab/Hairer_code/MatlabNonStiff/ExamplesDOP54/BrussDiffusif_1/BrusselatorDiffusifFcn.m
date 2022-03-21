function yPrime = BrusselatorDiffusifFcn(t,uv,OdeCoeff)
% Differential equations:
% du/dt = g1 = A + u^2* v -(B+1)*u + alpha*d2u/dx2
% dv/dt = g2 = B*u - U^2*v + a*d2v/dx2
% Parameters  A, B and alpha
% To get a 5 bands matrix, the equations are interleaved
% Set y0    = [u0(1);v0(1);u0(2);v0(2); ... ].
%     dy/dt = [du(1)/dt; dv(1)/dt;      ... ].
% The jacobian is a sparse matrix (5 bands).
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
xspan    = OdeCoeff{1};
A        = OdeCoeff{2};
B        = OdeCoeff{3};
alpha    = OdeCoeff{4};
uBord    = OdeCoeff{5};
vBord    = OdeCoeff{6};

dx    = xspan(2) - xspan(1);
adx2  = alpha/dx^2;

Nx    = length(xspan);
N     = length(uv);

u     = uv(1:2:N-1);
v     = uv(2:2:N);
f     = zeros(size(u));
g     = zeros(size(v));

% First two lines
f(1)  = A + u(1)^2*v(1) - (B+1)*u(1) + adx2*(uBord - 2*u(1) + u(2));
f(Nx) = A + u(Nx)^2*v(Nx) - (B+1)*u(Nx) + adx2*(u(Nx-1) - 2*u(Nx) + uBord);

for i = 2:Nx-1     
  f(i) = A + u(i)^2*v(i) - (B+1)*u(i) + adx2*(u(i-1) - 2*u(i) + u(i+1));
  g(i) = -u(i)^2*v(i) + B*u(i) + adx2*(v(i-1) -2*v(i) + v(i+1));
end

% Last two lines
g(1)  = -u(1)^2*v(1) + B*u(1) + adx2*(vBord -2*v(1) + v(2));
g(Nx) = -u(Nx)^2*v(Nx) + B*u(Nx) + adx2*(v(Nx-1) -2*v(Nx) + vBord);

yPrime(1:2:2*Nx-1) = f;
yPrime(2:2:2*Nx)   = g;
yPrime = yPrime(:);




