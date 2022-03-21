% Driver for the exercise page 296
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

clc
close all
clear all

addpath D:\RDPSolver

ProblemName = 'Ex_P296';

InfoStatDef  = true;
InfoDynDef   = true;

tv  = [0  160];
y0  = [63; 11; 11; 21];
I   = 10.5;

y0equ  = [ I*(1+0.004*I^3);  I; I; 2*I];

OdeDef   = @Ex_P296_Prime;
PhiDef   = @Ex_P296_Phi;

RelTolDef  = 1e-6;
AbsTolDef  = 1e-6;
RefineDef  = 1;

options = rdpset('RelTol',RelTolDef,'AbsTol',AbsTolDef);
options = rdpset(options,'Refine',RefineDef);
% ---------------------------

tic
[t,y,Stats] = dop54d(OdeDef,PhiDef,tv,y0equ,options,I,y0(4)); 
T_dopd  = toc;
L_dopd  = length(t);
TL_dopd = [T_dopd,L_dopd]
plot(t,y)
grid on
title('Equilibrium')

tic
[t,y,Stats] = dop54d(OdeDef,PhiDef,tv,y0,options,I,y0(4)); 
T_dopd  = toc;
L_dopd  = length(t);
TL_dopd = [T_dopd,L_dopd]

figure
plot(t,y)
grid on
title('Non-equilibrium')

tout     = t;
yout     = y;
InfoStat = Stats.Stat;
InfoDyn  = Stats.Dyn;
if InfoDynDef
  figure    
  subplot(3,1,1)
  plot(tout,real(yout(:,1)),'b')
  grid on
  hold on
  plot(tout,real(yout(:,1)),'.b')
 
  title('Ex page 296, y ')
    
  subplot(3,1,2)
  semilogy(InfoDyn.haccept_t,InfoDyn.haccept);
  hold on
  grid on
  semilogy(InfoDyn.haccept_t,InfoDyn.haccept,'.k');
  plot(InfoDyn.hreject_t,InfoDyn.hreject,'or')
  title('Length of steps in function of time, accept = black, reject = red ')   
  
  subplot(3,1,3)
  semilogy(InfoDyn.haccept_Step,InfoDyn.haccept);
  hold on
  grid on
  semilogy(InfoDyn.haccept_Step,InfoDyn.haccept,'.k');
  plot(InfoDyn.hreject_Step,InfoDyn.hreject,'or')
  title('Length of steps in function of step number, accept = black, reject = red ')     
end 

if InfoStatDef  
  InfoStat = InfoStat
end


