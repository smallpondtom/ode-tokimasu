% Driver for the exercise page 286
% Differential equation with delay
% y' = -y(t-1), 
% y(0) = 1; 
% y(t) = Phi(t) = 1, -1<=t<=0
% 
% Solution: 
% y(t) = 1-t;                                                    0<= t <= 1
% y(t) = 1-t + (t-1)^2/2;                                        1<= t <= 2
% y(t) = 1-t + (t-1)^2/2 - (t-2)^3/6;                            2<= t <= 3
% y(t) = 1-t + (t-1)^2/2 - (t-2)^3/6 + (t-3)^4/24;               3<= t <= 4
% y(t) = 1-t + (t-1)^2/2 - (t-2)^3/6 + (t-3)^4/24 - (t-4)^5/120; 4<= t <= 5
% etc.
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

ProblemName = 'Ex_P286';

InfoStatDef = true;
InfoDynDef  = true;

tv = [0  10];
y0 = 1;

OdeDef = @Ex_P286_Prime;
PhiDef = @Ex_P286_Phi;


RelTolDef      = 1e-6;
AbsTolDef      = 1e-6;
RefineDef      = 1;
InitialStepDef = 0.1;
MaxStepDef     = 0.1;

option = rdpset('RelTol',RelTolDef,'AbsTol',AbsTolDef);
option = rdpset(option,'Refine',RefineDef);
% option = rdpset(option,'InitialStep',InitialStepDef);
% option = rdpset(option,'MaxStep',MaxStepDef);
% ---------------------------
dop853d(OdeDef,PhiDef,tv,y0,option)

tic
[t1,y1,Stats] = dop853d(OdeDef,PhiDef,tv,y0,option);   
T_dopd  = toc;
L_dopd  = length(t1);
TL_dopd = [T_dopd,L_dopd]

for k = 1:L_dopd
  t = t1(k);
  if t <= 1
    yexact(k) = 1-t;
  elseif t <= 2
    yexact(k) = 1-t + (t-1)^2/2;
  elseif t <= 3
    yexact(k) = 1-t + (t-1)^2/2 - (t-2)^3/6;
  elseif t <= 4
    yexact(k) = 1-t + (t-1)^2/2 - (t-2)^3/6 + (t-3)^4/24;
  elseif t <= 5
    yexact(k) = 1-t + (t-1)^2/2 - (t-2)^3/6 + (t-3)^4/24 - (t-4)^5/120; 
  elseif t <= 6
    yexact(k) = 1-t + (t-1)^2/2 - (t-2)^3/6 + (t-3)^4/24 - (t-4)^5/120 ...
                    + (t-5)^6/720; 
  elseif t <= 7
    yexact(k) = 1-t + (t-1)^2/2 - (t-2)^3/6 + (t-3)^4/24 - (t-4)^5/120 ...
                    + (t-5)^6/720 - (t-6)^7/5040; 
  elseif t <= 8
    yexact(k) = 1-t + (t-1)^2/2 - (t-2)^3/6 + (t-3)^4/24 - (t-4)^5/120 ...
                    + (t-5)^6/720 - (t-6)^7/5040 + (t-7)^8/40320; 
  elseif t <= 9
    yexact(k) = 1-t + (t-1)^2/2 - (t-2)^3/6 + (t-3)^4/24 - (t-4)^5/120 ...
                    + (t-5)^6/720 - (t-6)^7/5040 + (t-7)^8/40320 ...
                    - (t-8)^9/362880;                     
  elseif t <= 10
    yexact(k) = 1-t + (t-1)^2/2 - (t-2)^3/6 + (t-3)^4/24 - (t-4)^5/120 ...
                    + (t-5)^6/720 - (t-6)^7/5040 + (t-7)^8/40320 ...
                    - (t-8)^9/362880 + (t-9)^10/3628800;   
  end
end

figure(1)
subplot(2,1,1)
plot(t1,y1);         
title(['dop853d,  Time = ',num2str(T_dopd)])
grid on
subplot(2,1,2)
plot(t1,y1-yexact')
grid on

tout     = t1;
yout     = y1;
InfoStat = Stats.Stat;
InfoDyn  = Stats.Dyn;
if InfoDynDef
  figure    
  subplot(3,1,1)
  plot(tout,real(yout(:,1)),'b')
  grid on
  hold on
  plot(tout,real(yout(:,1)),'.b') 
  title('Ex page 286, y ')
    
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


