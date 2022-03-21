function fg = BrusselatorDiffusifFcn(t,uv,OdeCoeff)
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
fg    = zeros(size(uv));

% First two lines
fg(1)   = A + uv(1)^2*uv(2) - (B+1)*uv(1) + adx2*(uBord - 2*uv(1) + uv(3));
fg(2)   = -uv(1)^2*uv(2) + B*uv(1) + adx2*(vBord -2*uv(2) + uv(4));

for i = 3:2:N-3     
  fg(i) = A + uv(i)^2*uv(i+1) - (B+1)*uv(i) + adx2*(uv(i-2) - 2*uv(i) + uv(i+2));
end
for i= 4:2:N-2
  fg(i) = -uv(i-1)^2*uv(i) + B*uv(i-1) + adx2*(uv(i-2) -2*uv(i) + uv(i+2));
end

% Last two lines
fg(N-1) = A + uv(N-1)^2*uv(N) - (B+1)*uv(N-1) + adx2*(uv(N-3) - 2*uv(N-1) + uBord);
fg(N)   = -uv(N-1)^2*uv(N) + B*uv(N-1) + adx2*(uv(N-2) -2*uv(N) + vBord);







