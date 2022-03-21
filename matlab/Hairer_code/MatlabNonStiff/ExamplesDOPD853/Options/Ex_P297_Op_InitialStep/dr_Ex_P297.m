% Driver for the exercise page 297
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

ProblemName = 'Ex_P297';

InfoStatDef  = true;
InfoDynDef   = true;

tv  = [0  60];
y0  = [1e-6; 1; 1; 0];

OdeDef   = @Ex_P297_Prime;
PhiDef   = @Ex_P297_Phi;

Coeff_h1 = [2, 0.8, 1e4, 0.17, 0.5, 10, 0.12, 8];

tic
[t,y,Stats] = dop853d(OdeDef,PhiDef,tv,y0,[],Coeff_h1); 
T_dopd  = toc;
L_dopd  = length(t);
TL_dopd = [T_dopd,L_dopd]

y(:,1) = y(:,1)*1e4;
y(:,2) = y(:,2)/2;
y(:,3) = y(:,3);
y(:,4) = y(:,4)*10;
figure
plot(t,y)
grid on
title(['Complete Recovery,  h6 = ',num2str(Coeff_h1(6)),'   InitialStep = Default'])

InitialStepDef    = 0.05;
options = rdpset('InitialStep',InitialStepDef);

tic
[t,y,Stats] = dop853d(OdeDef,PhiDef,tv,y0,options,Coeff_h1); 
T_dopd  = toc;
L_dopd  = length(t);
TL_dopd = [T_dopd,L_dopd]

y(:,1) = y(:,1)*1e4;
y(:,2) = y(:,2)/2;
y(:,3) = y(:,3);
y(:,4) = y(:,4)*10;
figure
plot(t,y)
grid on
title(['Complete Recovery,  h6 = ',num2str(Coeff_h1(6)),'   InitialStep = ',num2str(InitialStepDef)])

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
  title('Ex page 297, y(:,1) ')
    
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


