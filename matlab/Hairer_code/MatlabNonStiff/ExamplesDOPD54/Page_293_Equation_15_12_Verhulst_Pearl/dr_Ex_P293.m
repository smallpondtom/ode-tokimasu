% Driver for the exercise page 293
% Differential equation with delay
% y'   = -y(t-1), 
% y(0) = 1; 
% y(t) = 1, -1<=t<=0
% The time delay is 1
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

ProblemName = 'Ex_P293';

InfoStatDef  = true;
InfoDynDef   = true;

tv  = [0  60];
y0 = 0.1;

OdeDef   = @Ex_P293_Prime;
PhiDef   = @Ex_P293_Phi;
CoeffVec = [0.35 0.5 1 1.4 1.6];
ColorVec = ['b','r','g','k','m'];

RelTolDef  = 1e-6;
AbsTolDef  = 1e-6;

option = rdpset('RelTol',RelTolDef,'AbsTol',AbsTolDef);

% ---------------------------
figure(1)
tic
for k = 1:length(CoeffVec)
  Coeff = CoeffVec(k);
  [t,y,Stats] = dop54d(OdeDef,PhiDef,tv,y0,option,Coeff); 
  plot(t,y,ColorVec(k))
  hold on
end
T_dopd  = toc;
L_dopd  = length(t);
TL_dopd = [T_dopd,L_dopd]

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
 
  title('Ex page 293, y ')
    
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


