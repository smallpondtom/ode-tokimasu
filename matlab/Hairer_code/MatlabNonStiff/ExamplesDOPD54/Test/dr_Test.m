% Driver for a Test
% Differential equation with delay
% y' = 3*y(t) - 2*y(t-1) - 3*t^2 - 4*t + 7
% y(0) = 1
% y(t) = Phi(t) = 3t^2 - 2t + 1 
% The solution is analytic everywhere and the exact solution is:
% y =  3t^2 - 2t + 1
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

t  = [0,2];
y0 = 1;

AbsTolDef  = 1e-6;
RelTolDef  = 1e-3;
MaxStepDef = 0.1;
InitialStepDef = 0.1;

option = rdpset('RelTol',RelTolDef,'AbsTol',AbsTolDef);
% option = rdpset(option,'MaxStep',MaxStepDef);
% option = rdpset(option,'InitialStep',InitialStepDef);
% ---------------------------
tic
[t1,y1]    = dop54d(@Test_Prime,@Test_Phi,t,y0,option);   
T_dop54dN  = toc;
L_dop54dN  = length(t1);
TL_dop54dN = [T_dop54dN,L_dop54dN]
yExact1    = 3*t1.^2 - 2*t1 + 1;

figure(1)
subplot(3,1,1); plot(t1,y1);         
title(['dop54dN,  Time = ',num2str(T_dop54dN)])
grid on
subplot(3,1,2); plot(t1,yExact1); 
title('Exact sol')
grid on
subplot(3,1,3); plot(t1,y1-yExact1); 
title('Exact sol - dop54d sol')
grid on
