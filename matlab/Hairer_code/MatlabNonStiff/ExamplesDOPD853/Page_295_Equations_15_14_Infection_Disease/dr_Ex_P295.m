% Driver for the exercise page 295
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

ProblemName = 'Ex_P295';

InfoStatDef = true;
InfoDynDef   = true;

tv  = [0 : 0.1 : 40];
y0 = [5; 0.1;  1];

OdeDef   = @Ex_P295_Prime;
PhiDef   = @Ex_P295_Phi;

RelTolDef  = 1e-6;
AbsTolDef  = 1e-6;
RefineDef  = 1;

options = rdpset('RelTol',RelTolDef,'AbsTol',AbsTolDef);
options = rdpset(options,'Refine',RefineDef);
% ---------------------------

tic
[t,y,Stats] = dop853d(OdeDef,PhiDef,tv,y0,options); 
T_dopd  = toc;
L_dopd  = length(t);
TL_dopd = [T_dopd,L_dopd]

plot(t,y)
grid on

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
 
  title('Ex page 295, y ')
    
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


