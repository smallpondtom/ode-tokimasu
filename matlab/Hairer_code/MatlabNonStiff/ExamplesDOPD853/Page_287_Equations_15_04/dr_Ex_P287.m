% Driver for the exercise page 287
% Differential equation with delay
% y' = -y(t-1), y(0) = 1; 
% Three different Phi functions:
% 1) y(t) = Phi1(t) 0.8,    -1<=t<=0
% 2) y(t) = Phi2(t) 0.8+t   -1<=t<=0
% 3) y(t) = Phi3(t) 0.8+2*t -1<=t<=0
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

ProblemName = 'Ex_P287';

InfoStatDef  = true;
InfoDynDef   = true;

tv = [0  7];
y0 = 0.8;

OdeDef  = @Ex_P287_Prime;
PhiDef1 = @Ex_P287_Phi1;
PhiDef2 = @Ex_P287_Phi2;
PhiDef3 = @Ex_P287_Phi3;

RelTolDef      = 1e-6;
AbsTolDef      = 1e-6;
RefineDef      = 10;
InitialStepDef = 0.1;
MaxStepDef     = 0.1;

option = rdpset('RelTol',RelTolDef,'AbsTol',AbsTolDef,'MaxStep',MaxStepDef);
% option = rdpset(option,'Refine',RefineDef);
option = rdpset(option,'InitialStep',InitialStepDef);

% ---------------------------

tic
[t1,y1] = dop853d(OdeDef,PhiDef1,tv,y0,option);  
[t2,y2] = dop853d(OdeDef,PhiDef2,tv,y0,option);  
[t3,y3] = dop853d(OdeDef,PhiDef3,tv,y0,option);
T_dopd  = toc;
L_dopd  = length(t1);
TL_dopd = [T_dopd,L_dopd]

t    = [-1:0.1:0]';
for k = 1 : length(t)
  Phi1(k) = Ex_P287_Phi1(t(k));
end
t1   = [t;t1];
y1   = [Phi1';y1];

for k = 1 : length(t)
  Phi2(k) = Ex_P287_Phi2(t(k));
end
t2 = [t;t2];
y2 = [Phi2';y2];

for k = 1 : length(t)
  Phi3(k) = Ex_P287_Phi3(t(k));
end
t3 = [t;t3];
y3 = [Phi3';y3];

figure(1)
plot(t1,y1,'b'); 
hold on
plot(t2,y2,'r');
plot(t3,y3,'k');
grid on
title(['dop853d,  Time = ',num2str(T_dopd)])



  
