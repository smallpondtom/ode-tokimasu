% Driver for the Wille and Backer problem
% Differential equation with delay
% y(1)' = y1(t-1)  
% y(2)' = y1(t-1) + y_2(t-0.2)
% y(3)' = y2(t)
% y(0)  = [1,1,1]
% y(t) = Phi(t) = [1;1;1]  t <= 0
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

addpath D:\RDPSolver

ProblemName = 'Wiile_and_Backer';

InfoStatDef = true;
InfoDynDef  = true;

tv = [0  1];
y0 = [1,1,1];

OdeDef = @Wille_and_Backer_Prime;
PhiDef = @Wille_and_Backer_Phi;

RelTolDef      = 1e-8;
AbsTolDef      = 1e-8;
RefineDef      = 10;
InitialStepDef = 0.1;

option = rdpset('RelTol',RelTolDef,'AbsTol',AbsTolDef);
option = rdpset(option,'Refine',RefineDef);
% option = rdpset(option,'InitialStep',InitialStepDef);

tic
[t1,y1] = dop54d(OdeDef,PhiDef,tv,y0,option); 
T_dopd = toc;
L_dopd  = length(t1);
TL_dopd = [T_dopd,L_dopd]

figure(1)
subplot(3,1,1)
plot(t1,y1(:,1),'b',t1,y1(:,2),'r',t1,y1(:,3),'g')
title('Wille and Backer')
xlabel('time t');
ylabel('solution y');

% Exact solution

for k = 1: length(t1)
  t = t1(k);
  if t <= 1
    ye(k,1) = 1+t;
  elseif t <= 2
    ye(k,1) = 1.5 + 0.5 *t*t;
  elseif t <= 3
    ye(k,1) = 1/6*(t-1).^3 + 1.5*t + 1/3;
  elseif t <= 4
    ye(k,1) = 1/24*(t-2).^4 + 1.5/2*(t-1).^2 +t/3 + 2.125;
  elseif t > 4
    ye(k,1) = -1.13333333333333 + 1/6 *(t-1)^2 + 1.5/6*(t-2)^3 + ...
            1/120*(t-3)^5 + 2.125*t;        
  else
    ye(k,1) = 0;
  end
end

for k = 1: length(t1)
  t = t1(k);
  if t <= 0.2
    ye(k,2) = 1 + 2*t;
    ye(k,3) = 1 + t + t^2;
  elseif t <= 0.4
    ye(k,2) = 1 + 2*t + (t-0.2)^2;
    ye(k,3) = 1 + t + t^2 + (t-0.2)^3/3;
  elseif t <= 0.6
    ye(k,2) = 1 + 2*t + (t-0.2)^2   + (t-0.4)^3/3 ;
    ye(k,3) = 1 + t + t^2 + (t-0.2)^3/3 + (t-0.4)^4/12;
  elseif t <= 0.8
    ye(k,2) = 1 + 2*t + (t-0.2)^2 + (t-0.4)^3/3 + (t-0.6)^4/12;
    ye(k,3) = 1 + t + t^2 + (t-0.2)^3/3 + (t-0.4)^4/12 + (t-0.6)^5/60;    
  elseif t <= 1
    ye(k,2) = 1 + 2*t + (t-0.2)^2 + (t-0.4)^3/3 + (t-0.6)^4/12 + (t-0.8)^5 /60; 
    ye(k,3) = 1 + t + t^2 + (t-0.2)^3/3 + (t-0.4)^4/12 + (t-0.6)^5/60 + (t-0.8)^6/360;      
  else
    ye(k,2) = 0;
    ye(k,3) = 0;
  end  
end
   

subplot(3,1,2)
plot(t1,ye(:,1),'b',t1,ye(:,2),'r',t1,ye(:,3),'g')

subplot(3,1,3)
plot(t1,y1(:,1)-ye(:,1),'b',t1,y1(:,2)-ye(:,2),'r',t1,y1(:,3)-ye(:,3),'g')




