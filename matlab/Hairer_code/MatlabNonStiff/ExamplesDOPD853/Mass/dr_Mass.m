% Driver for the harmonic oscillator with constant
% Differential equation with delay
% yPrime(1) = y2;
% yPrime(2) = -y1(t-2*pi); 
% y(0) = [1;0]; 
% y1(t) = Phi1(t) = [cos(t),0];  t <= 0
% 
% Solution: 
% y1(t) = cos(t)
% y2(t) = -sin(t) 
% On one side, the mass function is defined via the option. This mass is 
% used normally and multiplie y'.
% On the other side, the same mass ente in the Mass_Prime function via
% the arguments "varargin". In the Mass_Prime function, the inverse
% of the mass multiply y'.
% Hence, fo every invertible mass, the solution is always:
% y1(t) = cos(t) and y2(t) = -sin(t).
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
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS 
% IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
% TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR 
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% ------------------------------------------------------------------------

clear all
close all
clc

t       = [0 2*pi];
y0      = [1,0];

% Mass = matrice constante

% MassCte = [1,0;0,1];

MaxStepDef = 1;
RefineDef  = 10;


option = rdpset('MaxStep', MaxStepDef);
option = rdpset(option,'Refine', RefineDef);


MassCte       = [2,-3;-1,4];
MassFcnCteDef = @MassFcnCte;
option        = rdpset(option,'MassFcn',MassFcnCteDef);
tic
[t1,y1] = dop853d(@MassCte_Prime,@Mass_Phi,t,y0,option,MassCte);
T_dopd  = toc;
L_dopd  = length(t1);
TL_dopd = [T_dopd,L_dopd]

figure(1)
plot(t1,y1)
title(['dop853d, Mass = cte, time = ',num2str(T_dopd)])


MassCte       = [2,-3;-1,4];
MassFcntDef   = @MassFcn_t;
option        = rdpset(option,'MassFcn',MassFcntDef);
tic
[t1,y1] = dop853d(@Mass_t_Prime,@Mass_Phi,t,y0,option,MassCte);
T_dopd  = toc;
L_dopd  = length(t1);
TL_dopd = [T_dopd,L_dopd]
figure(2)
plot(t1,y1)
title(['dop853d, Mass function of time, time = ',num2str(T_dopd)])


MassCte       = [2,-3;-1,4];
MassFcntyDef  = @MassFcn_ty;
option        = rdpset(option,'MassFcn',MassFcntyDef);
tic
[t1,y1] = dop853d(@Mass_ty_Prime,@Mass_Phi,t,y0,option,MassCte);
T_dopd  = toc;
L_dopd  = length(t1);
TL_dopd = [T_dopd,L_dopd]
figure(3)
plot(t1,y1)
title(['dop853d, Mass function of time, time = ',num2str(T_dopd)])



