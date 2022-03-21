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

tv  = [0  60];
y0  = [1e-6; 1; 1; 0];

OdeDef   = @Ex_P297_Prime;
PhiDef   = @Ex_P297_Phi;

Coeff_h1 = [2, 0.8, 1e4, 0.17, 0.5, 10, 0.12, 8];

OutputFcnDef1 = @WriteToFile1;
OutputSelDef1 = 1;
options       = rdpset('OutputFcn',OutputFcnDef1);
options       = rdpset(options,'OutputSel',OutputSelDef1);
dop54d(OdeDef,PhiDef,tv,y0,options,Coeff_h1); 

OutputFcnDef2 = @WriteToFile2;
OutputSelDef2 = [1,2];
options       = rdpset('OutputFcn',OutputFcnDef2);
options = rdpset(options,'OutputSel',OutputSelDef2);
dop54d(OdeDef,PhiDef,tv,y0,options,Coeff_h1); 


clear all

FileName1 = 'Ex_P297_Result_1.txt';
FileName2 = 'Ex_P297_Result_2.txt';

Fid1 = fopen(FileName1,'r');
[t_y1,Count1] = fscanf(Fid1, '%g %g', [2 inf]);
t_y1 = t_y1';
fclose(Fid1);

Fid2 = fopen(FileName2,'r');
[t_y2,Count2] = fscanf(Fid1, '%g %g %g',[3 inf]);
t_y2 = t_y2';
fclose(Fid2);

figure(1)
plot(t_y1(:,1),t_y1(:,2))
grid on

figure(2)
subplot(2,1,1)
plot(t_y2(:,1),t_y2(:,2))
grid on
subplot(2,1,2)
plot(t_y2(:,1),t_y2(:,3))
grid on





