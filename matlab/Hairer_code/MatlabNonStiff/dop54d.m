function [varargout] = dop54d(OdeFcn,PhiFcn,tspan,y0,options,varargin)
% ------------------------------------------------------------------------
% Solver Name :  dop54d.
% -----------
%
%     Numerical solution of a non-stiff system of first order ordinary 
%     differential equations with delay.
%     This is an explicit Runge-Kutta method of order (4)5 due to
%     Dormand & Prince (with step size control and dense output).
%
%
%     AUTHORS: E. HAIRER AND G. WANNER
%              UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
%              CH-1211 GENEVE 24, SWITZERLAND 
%              E-MAIL:  Ernst.Hairer@math.unige.chfbeta
%                       Gerhard.Wanner@math.unige.ch
%     THIS CODE IS DESCRIBED IN:
%         E. HAIRER, S.P. NORSETT AND G. WANNER, SOLVING ORDINARY
%         DIFFERENTIAL EQUATIONS I. NONSTIFF PROBLEMS. 2ND EDITION.
%         SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
%         SPRINGER-VERLAG (1993)       
%
%     Matlab version:
%     Denis Bichsel
%     Rue des Deurres 58
%     2000 Neuch�tel
%     Suisse
%     dbichsel@infomaniak.ch
%     Version of end 2015
% 
% Case 1)
%    dop54d(OdeFcn,PhiFcn,tspan,y0) or
%    dop54d(OdeFcn,PhiFcn,tspan,y0,options) 
%    dop54d(OdeFcn,PhiFcn,tspan,y0,options,Parameters)
%    solve the system of first order delay differential equations 
%    y'  = OdeFcn(t,y,varargin)
%    Phi = PhiFcn(t,y,varargin) defines the phi (y) value between
%    tLag and tspan(1)
%  Input:
%    tspan, a vector with 2 or more components. The components must be
%    monotonic ordered.
%    If tspan is a two components vector, the output value are only
%    returned at the calculated points if 'Refine' is <= 1. If refine is
%    bigger than 1, say 10, 9 more values are interpolated between two
%    succesive evaluated points.
%    If length(tspan) > 2, dop54d returns the solutions at these t-values
%    only.
%    y0, the initial condition must be a column vector.
%    Options overwrite the default integration parameters (see rdpset).
%    Parameters may be passed to OdeFcn, to MassFcn, and to EventsFcn.
%  Output:
%    In this cas, the default output is odeplot. The option 'outputFcn' may
%    overwrite the default and, for exmaple, write the results in a file.
%   
% Case 2)
%    [tout,yout] = dop54d(OdeFcn,PhiFcn,tspan,y0) or
%    [tout,yout] = dop54d(OdeFcn,PhiFcn,tspan,y0,options) 
%    [tout,yout] = dop54d(OdeFcn,PhiFcn,tspan,y0,options,Parameters)
%  Input:
%    Same as case 1).
%  Output:
%   The solution is returned in [tout,yout].
%  
% Case 3)
%    [tout,yout,Stats] = dop54d(OdeFcn,PhiFcn,tspan,y0) or
%    [tout,yout,Stats] = dop54d(OdeFcn,PhiFcn,tspan,y0,options) 
%    [tout,yout,Stats] = dop54d(OdeFcn,PhiFcn,tspan,y0,options,Parameters)
%  Input:
%    Same as case 1).
%  Output:
%    The solution is returned in [tout,yout].
%    Stats contains some informations on the calculation.
%      Stats.Stat gives the the following global (static) informations. 
%      Stats.Stat.FcnNbr:     The call number to OdeFcn. 
%      Stats.Stat.StepNbr:    The number of main steps. 
%      Stats.Stat.AccptNbr:   The number of accepted steps
%      Stats.Stat.StepRejNbr: The number of rejected steps
%
%    Stats.Dyn gives the following dynamical informations. 
%      Stats.Dyn.haccept_t:      Times of the accepted steps
%      Stats.Dyn.haccepted_Step: Accepted steps' number
%      Stats.Dyn.haccept:        Values of the accepted step sizes
%      Stats.Dyn.hreject_t:      Times of the rejected steps  
%      Stats.Dyn.hreject_Step:   Rejected steps' number
%      Stats.Dyn.hreject:        Rejected steps'value
%
% case 4)
%    The option Events is set, (EventsFcn ~= [] ) and there are five 
%    outputs.     
%    [tout,yout,te,ye,ie] = dop54d(OdeFcn,PhiFcn,tspan,y0) or
%    [tout,yout,te,ye,ie] = dop54d(OdeFcn,PhiFcn,tspan,y0,options) 
%    [tout,yout,te,ye,ie] = dop54d(OdeFcn,PhiFcn,tspan,y0,options,Parameters)
%  Input:
%    Same as case 1).
%  Output:
%    The solution is returned in [tout,yout]. 
%    te is a the t-vector of the events
%    ye contains the value of the events
%    ie is the indice of the events function corresponding to te, ye
%
% case 5)
%    The option Events is set, (EventsFcn ~= [] ) and there are five 
%    outputs.     
%    [tout,yout,te,ye,ie,Stats] = dop54d(OdeFcn,PhiFcn,tspan,y0) or
%    [tout,yout,te,ye,ie,Stats] = dop54d(OdeFcn,PhiFcn,tspan,y0,options) 
%    [tout,yout,te,ye,ie,Stats] = dop54d(OdeFcn,PhiFcn,tspan,y0,options,Parameters)
%  Input:
%    Same as case 1).
%  Output: same as case 4) plus Stats as in case 3)
%
% -------------------------------------------------------------------------          
% See DOP54 DOP853 DOP853D RADAU RDPGET RDPSET
% ------------------------------------------------------------------------
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS �AS 
% IS� AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
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
Solver_Name = 'dop54d';

Ny  = length(y0);

% ---------------------------
% Default options values
AbsTolDef           = 1e-6;
RelTolDef           = 1e-3;
InitialStepDef      = 1e-2;                  % WORK(7)
MaxStepDef          = tspan(end) - tspan(1); % WORK(6)
MaxNbrStepDef       = 1e5;                   % IWORK(1)
RefineDef           = 1;
OutputFcnDef        = [];
OutputSelDef        = 1:Ny;
FacLDef             = 0.2;                   % WORK(3)
FacRDef             = 10.0;                  % WORK(4)
SafeDef             = 0.9;                   % WORK(2)
NormControlDef      = false;
MassFcnDef          = [];
EventsFcnDef        = [];
StiffTestDef        = 1000;                  % IWORK(5)
BetaDef             = 0.04;                  % WORK(5)

OpDefault = {AbsTolDef;       RelTolDef;        InitialStepDef;  ...
             MaxStepDef;      MaxNbrStepDef;    ...
             RefineDef;       OutputFcnDef;     OutputSelDef;    ...
             FacLDef;         FacRDef;          SafeDef;         ...
             NormControlDef;  MassFcnDef;       EventsFcnDef;     ...
             StiffTestDef;    BetaDef}; 

OpNames = ['AbsTol          ';'RelTol          ';'InitialStep     '; ...
           'MaxStep         ';'MaxNbrStep      ';                    ...
           'Refine          ';'OutputFcn       ';'OutputSel       '; ...
           'FacL            ';'FacR            ';'Safe            '; ...
           'NormControl     ';'MassFcn         ';'EventsFcn       ';
           'StiffTest       ';'Beta            ']; ...        

% ---------------------------
% Tests on inputs
% ---------------
if (nargin < 4 )    %  dopxx('odefile',phifile',tspan,y0) 
  error ([Solver_Name,': Number of input arguments must be at least equal to 4.']);
end   

% nargin >= 4
if ~isa (OdeFcn, 'function_handle') && ~(exist(OdeFcn,'file') == 2)
  error ([Solver_Name,': First input argument must be a valid function handle or name']);        
end

if ~isa (PhiFcn, 'function_handle') && ~(exist(OdeFcn,'file') == 2)
  error ([Solver_Name,': First input argument must be a valid function handle or name']);        
end

if (~isvector (tspan) || length (tspan) < 2)
  error ([Solver_Name,': Second input argument must be a valid vector']);
end
% Test of tspan monotony
tspan  = tspan(:);                     % Column vector
PosNeg = sign(tspan(end) - tspan(1));  % Monotony check 
if any(PosNeg*diff(tspan)) <= 0
  error([Solver_Name,': Time vector must be strictly monotonic']);
end

if ~isvector (y0)
  error ([Solver_Name,': Initial conditions argument must be a valid vector or scalar']);
end
y0 = y0(:);           % Column vector

if nargin < 5
  options = [];
end
  
Arg.In = nargin > 5;

% ---------------------------
% Tests on options
% ----------------
Op  = [];
for n = 1:size(OpNames,1)
  Op.(deblank(OpNames(n,:))) =  ...
             rdpget(options,deblank(OpNames(n,:)),OpDefault{n});  
end 

% ------- ABSTOL
if ~isnumeric(Op.AbsTol)
   error([Solver_Name, ': Wrong input "AbsTol must be a positive number" ']) 
end
if any(Op.AbsTol <= 0 )
  error([Solver_Name,': Absolute tolerance are too small.']);
end
if (length(Op.AbsTol) ~= Ny) && (length(Op.AbsTol) ~=1)
  error([Solver_Name, ': AbsTol vector of length 1 or %d.',num2str(Ny)]);
end
Op.AbsTol = Op.AbsTol + zeros(size(y0));
% ------- RELTOL
if ~isnumeric(Op.RelTol)
  error([Solver_Name, ': Wrong input "RelTol" must be a positve number" ']) 
end
if any(Op.RelTol < 10*eps)  
  error([Solver_Name,': Relative tolerance are too small.']);
end
if (length(Op.RelTol) ~= Ny) && (length(Op.RelTol) ~=1)
  error([Solver_Name, ': RelTol vector of length 1 or %d.',num2str(Ny)]);
end
Op.RelTol = Op.RelTol + zeros(size(y0));
% ------- INITIAL STEP SIZE
if ~isnumeric( Op.InitialStep)
  error([Solver_Name, ': Wrong input "InitialStep" must be a number']) 
end 
% ------- MAXIMAL STEP SIZE
if ~isnumeric(Op.MaxStep)
  error([Solver_Name, ': Wrong input "MaxStep" must be a number '])
end
% ------- MAXIMAL NUMBER OF STEPS
if ~isnumeric(Op.MaxNbrStep)
  error([Solver_Name,': Wrong input "MaxNbrStep" must be a positive number']);
elseif Op.MaxNbrStep <= 0
  error([Solver_Name,': Wrong input "MaxNbrStep" ',num2str(Op.MaxNbrStep),', must be > 0'])
end
% ------- REFINE
if ~isempty(Op.Refine)
  if ~isnumeric(Op.Refine)
    error([Solver_Name,': Wrong input "Refine" must empty be or a positive number']);
  end
end
% ------- OUTPUTFCN 
if ~isempty(Op.OutputFcn)
  if isa (Op.OutputFcn,'function_handle')
    Op.OutputFcn = func2str(Op.OutputFcn);
  end
  if ~(exist(Op.OutputFcn,'file') == 2)
    error ([Solver_Name,': OutputFcn must be a valid function handle']);  
  end  
end
% ------- OUTPUTSEL  
if ~isempty(Op.OutputSel) 
  if ~isvector (Op.OutputSel)
    error ([Solver_Name,': OutputSel must be a scalar or a vector of positive integer']);  
  end 
  IndComp = 1:1:Ny;
  for n = 1:length(Op.OutputSel)
    if ~ismember(Op.OutputSel(n),IndComp);
      error ([Solver_Name,': OutputSel must be an integer in 1 .. ',num2str(Ny)])
    end
  end    
end
% -------  FACL,FACR     PARAMETERS FOR STEP SIZE SELECTION
if ~isnumeric(Op.FacL)
  error([Solver_Name, ': Wrong input "Facl" must be numeric default 0.2 '])
end
if Op.FacL > 1.0 
   error([Solver_Name, ': Curious input for "FacL" default 0.2 '])
end  
if ~isnumeric(Op.FacR)
  error([Solver_Name, ': Wrong input "FacR" must be numeric default 8 '])  
end
if Op.FacR < 1.0 
   error([Solver_Name, ': Curious input for "FacR" default 8 '])
end  
% ------- SAFE SAFETY FACTOR IN STEP SIZE PREDICTION   
if ~isnumeric(Op.Safe)
   error([Solver_Name, ': Wrong input "Safe" must be a positive number '])  
elseif Op.Safe <= 0.001 || Op.Safe >= 1 
  error([Solver_Name,': Curious input for Safe,' num2str(Op.Safe),' must be in ]0.001 .. 1[ '])
end
% ------- NORMCONTROL
if ~islogical(Op.NormControl) 
  error ([Solver_Name,': NormControl must be logical']);
end
% ------- MASS  

% ------- EVENTS  
EventsExist = ~isempty(Op.EventsFcn);
% ------- STIFFTEST
if ~isnumeric(Op.StiffTest)
   error([Solver_Name, ': Wrong input "StiffTest" must be a positive number ']) 
   if Op.StiffTest < 1
     error([Solver_Name, ': Wrong input "StiffTest" must be a positive number '])
   end
end
% ------- BETA
if ~isnumeric(Op.Beta)
   error([Solver_Name, ': Wrong input "Beta" must be a positive number 0 - 0.2']) 
elseif Op.Beta > 0.2
   error([Solver_Name, ': Curious input "Beta" must be a positive number 0 - 0.2']) 
elseif Op.Beta <= 0.0 
  Op.Beta = 0; 
end

% ---------------------------
% Tests on outputs
% ---------------
OutputNbr = abs(nargout());
if OutputNbr == 0  
  Op.Stats = false;
  if isempty(Op.OutputFcn)
    Op.OutputFcn = rdpget(options,'OutputFcn','odeplot');       
  end  
elseif OutputNbr == 2 
  Op.Stats = false;
elseif OutputNbr == 3
  Op.Stats = true;
elseif OutputNbr == 5 
  Op.Stats = false;
  if ~EventsExist
     error([Solver_Name,': Events not set, too much output']);
  end
elseif OutputNbr == 6
  Op.Stats = true;
  if ~EventsExist 
    error([Solver_Name,': Events and Stats must be set for 6 outputs']);
  end
else
  error([Solver_Name,': Outputs number not correct']);
end
    
solver54d = {'dop54dsolver',OdeFcn,PhiFcn,tspan,y0,Op};
if Arg.In
  solver54d = [solver54d, {varargin{:}}];
end 

switch OutputNbr
  case 0
    feval(solver54d{:});
  case 2
    [tout,yout] = feval(solver54d{:});
    varargout = {tout,yout};
  case 3
    [tout,yout,Stats] = feval(solver54d{:});
    varargout = {tout,yout,Stats};
  case 5
    [tout,yout,teout,yeout,ieout] = feval(solver54d{:});
    varargout = {tout,yout,teout,yeout,ieout};
  case 6
    [tout,yout,teout,yeout,ieout,Stats] = feval(solver54d{:});
    varargout = {tout,yout,teout,yeout,ieout,Stats};
end
return
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
function [varargout] = dop54dsolver(OdeFcn,PhiFcn,tspan,y,Op,varargin)

global tcont54Vec ycont54Vec hcont54Vec

Solver_Name = 'dop54dsolver';

% ------- INPUT PARAMETERS
% Time properties
tspan  = tspan(:);                     % Column vector
ntspan = length(tspan);
t      = tspan(1);
tfinal = tspan(end);
PosNeg = sign(tfinal-t);

% Number of equations, y is a column vector
Ny = length(y);

% General options
AbsTol      = Op.AbsTol;
RelTol      = Op.RelTol;
h           = Op.InitialStep;         % h may be positive or negative
hmax        = abs(Op.MaxStep);        % hmax is positive
MaxNbrStep  = Op.MaxNbrStep;
Refine      = Op.Refine;
OutputFcn   = Op.OutputFcn;
OutputSel   = Op.OutputSel;
FacL        = 1/Op.FacL;
FacR        = 1/Op.FacR;
Safe        = Op.Safe;                % WORK(2)
NormControl = Op.NormControl;
MassFcn     = Op.MassFcn;
EventsFcn   = Op.EventsFcn;
StiffTest   = Op.StiffTest;
Beta        = Op.Beta;
% Initialisation of Stat parameters
Stat.FcnNbr     = 0;
Stat.StepNbr    = 0;
Stat.StepAccNbr = 0;
Stat.StepRejNbr = 0;

% Initialisation of Dyn parameters
StatsExist = false;
if nargout == 3 || nargout == 6
  StatsExist       = true;
  Dyn.haccept_t    = [];
  Dyn.haccept_Step = [];
  Dyn.haccept      = [];
  Dyn.hreject_t    = [];
  Dyn.hreject_Step = [];
  Dyn.hreject      = [];  
end
% -------  ARGUMENTS
Arg.In = nargin > 4;

% ------- ODEFCN arguments or not
if isa (OdeFcn, 'function_handle')
  OdeError  = func2str(OdeFcn);
else
  OdeError  = OdeFcn;
end 
Arg.Ode = abs(nargin(OdeFcn)) > 2;
if Arg.Ode && ~Arg.In
  error([Solver_Name,':  ',OdeError,', parameters are missing'])
end

if isa (PhiFcn, 'function_handle')
  PhiError  = func2str(PhiFcn);
else
  PhiError  = PhiFcn;
end 
Arg.Phi = abs(nargin(PhiFcn)) >= 2;
if Arg.Phi && ~Arg.In
  error([Solver_Name,':  ',PhiError,', parameters are missing'])
end

% ------- MASSFCN arguments or not
MassExist  = false;
MassArgNbr = -1;
if ~isempty(MassFcn) 
  if isa (MassFcn, 'function_handle')
    MassError  = func2str(MassFcn);
  else
    MassError  = MassFcn;
  end    
  MassExist  = true;
  MassArgNbr = abs(nargin(MassFcn));
  switch MassArgNbr
    case 0      
      Mass    = MassFcn(); 
      detMass = det(Mass);
      if detMass == 0
        error([Solver_Name,':  Mass is singular '])
      end    
    otherwise
      if MassArgNbr >= 3 && ~Arg.In
        error([Solver_Name,':  ',MassError,', parameters are missing'])
      end     
  end  
end

% Set the output flag and output buffer
if ntspan == 2                                               
  if Refine <= 1                       % Computed points only   
    OutFlag = 1; 
    nBuffer = 100;
    nout    = 0;
    tout    = zeros(nBuffer,1);
    yout    = zeros(nBuffer,Ny);
  else
    OutFlag = 2;    
    nBuffer = 10*Refine;
	nout    = 0;
    tout    = zeros(nBuffer,1);
    yout    = zeros(nBuffer,Ny);% Computed + refined points    
  end
else  % ntspan > 2
  OutFlag = 3;                         % Computed points  
  nout    = 0;
  nout3   = 0;
  tout    = zeros(ntspan,1);
  yout    = zeros(ntspan,Ny); 
  if Refine > 1
    Refine = 1;
    warning([Solver_Name,': Refine set equal 1, because length(tspan) > 2 '])
  end
end

OutputFcnExist = false;
if ~isempty(OutputFcn) 
  OutputFcnExist = true;   
  % Initialize the OutputFcn
  OutputFcnArg = {OutputFcn,[t tfinal],y(OutputSel),'init'};
  feval(OutputFcnArg{:}); 
end

% Initialiation of internal parameters
FacOld = 1e-4;
Expo1  = 0.2 - Beta*0.75;             
hLamb  = 0.0;
IaSti  = 0;
K      = zeros(Ny,7);  % For DOP54 a is a matrix(7,7)
y1     = zeros(Ny,7);
% --------------
% Integration step

% First, tLag tLagTab
if Arg.Phi
  [Phi,tLag,tLagTab] = feval(PhiFcn,t,y,varargin{:});
else
  [Phi,tLag,tLagTab] = feval(PhiFcn,t,y);
end
% Second Phi value
NtLag = length(tLag);
yLag = zeros(Ny,NtLag);
for ny = 1:Ny
  for nt = 1:NtLag
    if tLagTab(ny,nt)
      yLag(ny,nt) = feval('dop54d_Lag',PhiFcn,ny,t+tLag(nt),y,varargin{:});
    end
  end
end
if any(isnan(yLag))
  error([Solver_Name,':  Function ',PhiError,' Some components of are NaN'])
end

if Arg.Ode
  f0 = feval(OdeFcn,t,y,yLag,varargin{:});
else
  f0 = feval(OdeFcn,t,y,yLag);
end
if any(isnan(f0))
  error([Solver_Name,':  Function ',OdeError,' Some components of are NaN'])
end
Stat.FcnNbr = Stat.FcnNbr+1;

[m,n]      = size(f0);
if n > 1
  error([Solver_Name,':  Function ',OdeError,'(t,y) must return a column vector.'])
elseif m ~= Ny
  error([Solver_Name,':  Vector ',OdeError, '(t,y) must be same length as initial conditions.']);
end

if MassExist                % There is a Mass
  if MassArgNbr > 0
    switch MassArgNbr     
      case 1
        MassFcnArg = {MassFcn,t};
      case 2
        MassFcnArg = {MassFcn,t,y};
      otherwise
        MassFcnArg = {MassFcn,t,y,varargin{:}};
    end
    Mass = feval(MassFcnArg{:}); 
    if any(isnan(Mass))
      error([Solver_Name, ': Some components of Mass are NAN'])       
    end
  end
  K(:,1) = Mass \ f0;
else
  K(:,1) = f0;
end

% hmax must be <= min(abs(tLag)) 
tLagMin = min(abs(tLag));
tLagMax = max(abs(tLag));
hmax    = min(hmax,tLagMin);

if abs(h) <= 10*eps
  h = hInitFcn(OdeFcn,MassFcn,MassArgNbr,Arg,PhiFcn,tLag,tLagTab, ...
               t,y,PosNeg,hmax,K(:,1),RelTol,AbsTol,varargin{:});
  Stat.FcnNbr = Stat.FcnNbr+1;
end
h      = PosNeg * min(abs(h),hmax);
Reject = false;

% Important
% ---------
% The results are much better when the steps don't cross the time lags.
% Construction of the next time's lags, t0 + k * tLag
% For best precision, the end of the step must be as near as possible 
% of a time LAG.
% It's very important when the derivative are not continuous. 

tStopVec = [];
for k = 1:length(tLag)
  ts       = tspan(1):-tLag(k):tspan(end);
  tStopVec = [tStopVec,ts];  
end
tStopVec = unique([tStopVec,tspan(end)]);

EventsExist = false;
if ~isempty(EventsFcn) 
  EventsExist = true;
  if nargin(EventsFcn) == 0
    error([Solver_Name,'  EventsFcn needs input arguments'])
  end 
  teout = [];
  yeout = [];
  ieout = [];
  [teout,yeout,ieout] = EventZeroFcn(EventsFcn,t,h,y,[],f0,'init',varargin{:});
end

% ---------------------------
% Coefficients values for dop54d
% See Hairer Universit� de Gen�ve
% ---------------------------
% c: time increment coefficients
c     = zeros(1,8); 
c(1)  = 0;
c(2)  = 0.2;
c(3)  = 0.3; 
c(4)  = 0.8;
c(5)  = 8/9;
c(6)  = 1;
c(7)  = 1;
% ---------------------------
% a: slope calculation coefficients 
a      =  zeros(7,7);
a(2,1) =  0.2;
a(3,1) =  3/40;
a(3,2) =  9/40;
a(4,1) =  44/45;
a(4,2) = -56/15;
a(4,3) =  32/9;
a(5,1) =  19372/6561;
a(5,2) = -25360/2187;
a(5,3) =  64448/6561;
a(5,4) = -212/729;
a(6,1) =  9017/3168;
a(6,2) = -355/33;
a(6,3) =  46732/5247;
a(6,4) =  49/176;
a(6,5) = -5103/18656;
a(7,1) =  35/384;
a(7,3) =  500/1113;
a(7,4) =  125/192;
a(7,5) = -2187/6784;
a(7,6) =  11/84;
% ---------------------------
% Dense output coefficients
d = zeros(1,7);
d(1) = -12715105075/11282082432;
d(3) =  87487479700/32700410799;
d(4) = -10690763975/1880347072;
d(5) =  701980252875/199316789632;
d(6) = -1453857185/822651844;
d(7) =  69997945/29380423;
% ---------------------------
% Error calculation coefficients
er(1) =  71/57600;
er(2) =  0.0;
er(3) = -71/16695;
er(4) =  71/1920;
er(5) = -17253/339200;
er(6) =  22/525;
er(7) = -1/40;
% ---------------------------

% ------------
% --- BASIC INTEGRATION STEP  
% ------------

Done   = false;
tcont54Vec = [];

while ~Done 

  % -------------------------
  % THE FIRST 6 STAGES
  % ------------------------- 
  Stat.StepNbr = Stat.StepNbr+1; 
 
  if Stat.StepNbr > MaxNbrStep
    warning([Solver_Name,':  NbrSteps = ',num2str(Stat.StepNbr),' > MaxNbrStep']);
    break
  end

  if (0.1*abs(h) <= abs(t)*eps)
    warning([Solver_Name,' : Too small step size']);    
    Done = true;
  end
 
  if t == tStopVec(1) && length(tStopVec) >= 2 
    tStopVec = tStopVec(2:end);
  end
  if (t + 1.01*h - tStopVec(1) ) >= 0
    h = tStopVec(1)-t;
  elseif t+1.8*h - tStopVec(1) >= 0
    h = ( tStopVec(1) - t )* 0.5;
  end  
  
  ch = h*c;  
  ah = h*a';     % Needed for matrix calculation  

  for j = 2:7  
    time    = t+ch(j);
    y1(:,j) = y+K*ah(:,j);    
    
    for ny = 1:Ny
      for nt = 1:NtLag
        if tLagTab(ny,nt)
          yLag(ny,nt) = feval('dop54d_Lag',PhiFcn,ny,time+tLag(nt),y1(:,j),varargin{:});
        end
      end
    end
    if any(isnan(yLag))
      error([Solver_Name, ': Some components of yLag = PhiFcn are NAN'])       
    end
      
    OdeFcnArg = {OdeFcn,time,y1(:,j),yLag};
    if (Arg.Ode)
      OdeFcnArg = [OdeFcnArg,{varargin{:}}];
    end 
    f0 = feval(OdeFcnArg{:});
    if any(isnan(f0))
      error([Solver_Name, ': Some components of f0 = OdeFcn are NAN'])       
    end
    if MassExist                % There is a Mass
      if MassArgNbr > 0
        switch MassArgNbr
          case 1
            MassFcnArg = {MassFcn,time};
          case 2
            MassFcnArg = {MassFcn,time,y1(:,j)};
          otherwise
            MassFcnArg = {MassFcn,time,y1(:,j),varargin{:}};
        end
        Mass = feval(MassFcnArg{:}); 
        if any(isnan(Mass))
          error([Solver_Name, ': Some components of Mass are NAN'])       
        end
      end             
      K(:,j) = Mass \ f0; 
    else        
      K(:,j) = f0;
    end              
  end;
 
  Stat.FcnNbr = Stat.FcnNbr+6;
  ySti = y1(:,6);  
  % K2 in Hairer fortran -->  K(:,7)
  K4   = h * K * er';      %  K4 ~= K(:,4)
  tph  = t + h;
  % --- ERROR ESTIMATION   (450)
  if NormControl
    %  norm(e) <= max(RelTol*norm(y),AbsTol)
    Sk = max(AbsTol) + max(RelTol)* max(norm(y),norm(y1(:,7)));
    Err = norm(K4)/Sk;
  else 
    Sk   = AbsTol + RelTol .* max( abs(y),abs(y1(:,7)));         
    Err  = sqrt( sum((K4./ Sk).^2)/Ny );
  end
  % --- COMPUTATION OF HNEW -----> 662 Hairer
  Fac11 = Err^Expo1;
  % --- LUND-STABILIZATION
  Fac   = Fac11/FacOld^Beta;
  % --- WE REQUIRE  FAC1 <= HNEW/H <= FAC2
  Fac   = max(FacR,min(FacL,Fac/Safe));
  hNew  = h/Fac;  
  
  if(Err < 1.D0)
    % --- STEP IS ACCEPTED          (470 Hairer)
    Stat.StepAccNbr = Stat.StepAccNbr+ 1;
    if StatsExist
      Dyn.haccept_t    = [Dyn.haccept_t,t];
      Dyn.haccept_Step = [Dyn.haccept_Step,Stat.StepNbr];
      Dyn.haccept      = [Dyn.haccept,h];
    end
    FacOld = max(Err,1e-4);
    % ------- STIFFNESS DETECTION                     675
    if mod(Stat.StepAccNbr,StiffTest)== 0 || IaSti > 0        
      StNum = sum( (K(:,7) - K(:,6)).^2 );
      StDen = sum( (ySti   - y1(:,7)).^2 ); 
      if StDen > 0 
        hLamb = abs(h) * sqrt(StNum/StDen);
      end         
      if hLamb > 3.25      
        StiffTest = 0;
        IaSti    = IaSti + 1;           
        if IaSti == 15
          warning([' The problem seems to become stiff at t = ', num2str(t)]);                  
        end
      else
        StiffTest = StiffTest + 1;
        if StiffTest == 6 
          IaSti = 0;
        end
      end
    end        
    
    % ------- FINAL PREPARATION FOR DENSE OUTPUT     495       
    YDiff      = y1(:,7)  - y;
    Bspl       = h*K(:,1) - YDiff;
    ycont(:,1) = y;
    ycont(:,2) = YDiff;
    ycont(:,3) = Bspl;
    ycont(:,4) = -h*K(:,7)+YDiff-Bspl;
    ycont(:,5) = h* K*d';
    
    if EventsExist
      [te,ye,ie,Stop] = EventZeroFcn(EventsFcn,t,h,y1(:,7),ycont,[],'',varargin{:});
      if ie > 0
        teout = [teout;te];
        yeout = [yeout;ye];
        ieout = [ieout;ie];
        if Stop 
          t = te;
          y = ye; 
          break 
        end
      end      
    end  
    
    if nargout > 0
      switch OutFlag        
        case 1          % Computed points, no Refinement 
          nout = nout + 1;
          if nout > length(tout)
            tout = [tout;zeros(nBuffer,1)];
            yout = [yout;zeros(nBuffer,Ny)];  
          end
          tout(nout)   = t;
          yout(nout,:) = y';                                    
        case 2          % Computed points, with refinement   
          nout         = nout + 1;
          tout(nout)   = t;
          yout(nout,:) = y';
          oldnout      = nout;
          nout         = nout + Refine - 1;
          S            = (1:Refine-1)' / Refine;
          if nout > length(tout)
            tout         = [tout; zeros(10*Refine,1)]; 	
            yout         = [yout; zeros(10*Refine,Ny)];
          end
          ii           = oldnout+1:nout;
          tinterp      = t+h*S;
          yinterp      = ntrprad(tinterp,t,h,ycont);
          tout(ii)     = tinterp;
          yout(ii,:)   = yinterp;
        case 3          % Output only at tspan points
          while ( PosNeg > 0 && t<= tspan(nout+1) && tspan(nout+1) < tph || ...
                  PosNeg < 0 && t>= tspan(nout+1) && tspan(nout+1) > tph)
            nout         = nout + 1;         
            yinterp      = ntrprad(tspan(nout),t,h,ycont);     
            tout(nout)   = tspan(nout);  
            yout(nout,:) = yinterp';  % Column output
          end      
      end
    end 

    if ~isempty(OutputFcn)
      switch OutFlag      % Output function required
        case 1            % Computed points, no Refinement 
          feval(OutputFcn,t,y(OutputSel)','');
        case 2            % Computed points, with refinement 
          [tout2,yout2] = OutFcnSolout2(t,h,y,ycont,OutputSel,Refine);
          for k = 1 : length(tout2)
            feval(OutputFcn,tout2(k),yout2(k,:),'');                    
          end       
        case 3            % Output only at tspan points
          [nout3,tout3,yout3] = OutFcnSolout3(nout3,t,h,ycont,OutputSel,tspan);          
          for k = 1 : length(tout3)
            feval(OutputFcn,tout3(k),yout3(k,:),'');                    
          end
      end  
    end
    
    % ---------------
    % We must save ycont to be able to evaluate y(t-tau)
    % even if tspan = [t0,tfinal] because we need to 
    % evaluate y at points like t0+h-tau  (time lag) at each step. 
    IndMax               = length(tcont54Vec) + 1;
    tcont54Vec(IndMax)     = t; 
    ycont54Vec(IndMax,:,:) = ycont;
    hcont54Vec(IndMax)     = h;   
    % Free no more useful tcondVec ycont54Vec and hcondVec components    
    Ltcont54Vec = length(tcont54Vec);
    Ind     = 1;      
    while t - tLagMax > tcont54Vec(Ind) 
      Ind = Ind + 1;
    end
    if Ind > 1
      Ind      = Ind - 1;
      ycont54Vec = ycont54Vec(Ind:Ltcont54Vec,:,:);
      tcont54Vec = tcont54Vec(Ind:Ltcont54Vec);
      hcont54Vec = hcont54Vec(Ind:Ltcont54Vec);
    end
    
    K(:,1) = K(:,7);
    t      = tph;
    y      = y1(:,7);  
    
    if t == tfinal  
      Done    = true;
    end           
    
    if abs(hNew) > hmax
      hNew = PosNeg*hmax;
    end
    if Reject
      hNew = PosNeg*min(abs(hNew),abs(h));
    end
	     
    Reject = false;
	    
  else % --- STEP IS REJECTED      depuis 457    (769 Hairer)
 
    hNew = h/min(FacL,Fac11/Safe);
    Reject = true;
    if Stat.StepAccNbr > 1 
      Stat.StepRejNbr = Stat.StepRejNbr + 1;
    end       
    if StatsExist
      Dyn.hreject_t    = [Dyn.hreject_t,t];
      Dyn.hreject_Step = [Dyn.hreject_Step,Stat.StepNbr];
      Dyn.hreject      = [Dyn.hreject,h];
    end 
  end
  h = hNew;

end  % while

if StatsExist
  Stats.Stat = Stat;
  Stats.Dyn  = Dyn;
end

OutputNbr = abs(nargout);
if OutputNbr > 0 
  nout         = nout + 1;
  tout(nout)   = t;
  yout(nout,:) = y(OutputSel)';
  tout         = tout(1:nout);
  yout         = yout(1:nout,:);
  varargout    = {tout,yout};
  if EventsExist    
    varargout = [varargout,{teout,yeout,ieout}];
  end
  if StatsExist
    varargout = [varargout,{Stats}];
  end      
end  
 
if OutputFcnExist   % Close the OutputFcn        
  OutputFcnArg = {OutputFcn,t,y(OutputSel)',''};
  feval(OutputFcnArg{:}); 
  OutputFcnArg = {OutputFcn,t,y(OutputSel)','done'};
  feval(OutputFcnArg{:});
end
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

% Computed points, with Refinement
function [tout,yout] = OutFcnSolout2(t,h,y,ycont,OutputSel,Refine)
tout       = zeros(Refine,1); 	
yout       = zeros(Refine,length(OutputSel));
tout(1)    = t;
yout(1,:)  = y(OutputSel)';
ii         = 2:Refine;
S          = (1:Refine-1)' / Refine;
tinterp    = t+h*S;
yinterp    = ntrprad(tinterp,t,h,ycont);
tout(ii)   = tinterp;
yout(ii,:) = yinterp(:,OutputSel);

% Output only at tspan points 
function [nout,tout,yout] = OutFcnSolout3(nout,t,h,ycont,OutputSel,tspan)
PosNeg = sign(tspan(end) - tspan(1));
tph    = t+h;
tout   = [];
yout   = [];
Compt  = 0;
while ( PosNeg > 0 && t <= tspan(nout+1) && tspan(nout+1) < tph || ...
        PosNeg < 0 && t>= tspan(nout+1) && tspan(nout+1) > tph)
  Compt         = Compt + 1 ;
  nout          = nout + 1;
  yinterp       = ntrprad(tspan(nout),t,h,ycont);     
  tout(Compt)   = tspan(nout); 
  yout(Compt,:) = yinterp(OutputSel)';  % Column output
end

function yinterp = ntrprad(tinterp,t,h,ycont)
S       = (tinterp-t)/h;   % S is theta in the book
S1      = 1-S;
for k = 1:length(S)
  yinterp(k,:) = ycont(:,1) + S(k)*(ycont(:,2) + S1(k)*(ycont(:,3) + ...
                 S(k)*(ycont(:,4) + S1(k)*ycont(:,5))));
end
return

function h = hInitFcn(OdeFcn,MassFcn,MassArgNbr,Arg,PhiFcn,tLag,tLagTab, ...
             t,y,PosNeg,hmax,f0,RelTol,AbsTol,varargin)                          
% ----------------------------------------------------------
% ----  COMPUTATION OF AN INITIAL STEP SIZE GUESS
% ----------------------------------------------------------      
% ---- COMPUTE A FIRST GUESS FOR EXPLICIT EULER AS
% ----   H = 0.01 * NORM (Y0) / NORM (F0)
% ---- THE INCREMENT FOR EXPLICIT EULER IS SMALL
% ---- COMPARED TO THE SOLUTION


% ------- MASSFCN arguments or not
Function_Name = 'hInitFcn';

Ordre = 5;

Sk  = AbsTol + RelTol.*abs(y);
Dnf = sum( (f0./Sk).^2 );
Dny = sum( (y./Sk).^2 );
if (Dnf < 1e-10 || Dny < 1e-10)
  h = 1e-6;
else
  h = sqrt(Dny/Dnf) * 0.01; 
end
h = min(h,hmax) * PosNeg;

% ---- PERFORM AN EXPLICIT EULER STEP
y1   = y + h*f0;
time = t+h;

Ny = length(y);
NtLag = length(tLag);
yLag = zeros(Ny,NtLag);
for ny = 1:Ny
  for nt = 1:NtLag
    if tLagTab(ny,nt)
      yLag(ny,nt) = feval('dop54d_Lag',PhiFcn,ny,time+tLag(nt),y1,varargin{:});
    end
  end
end

OdeFcnArg = {OdeFcn,t+h,y1,yLag};
if ~isempty(Arg.Ode)
  OdeFcnArg = [OdeFcnArg, {varargin{:}}];
end 
f1 = feval(OdeFcnArg{:}); 
if any(isnan(f1))
  error([Solver_Name, ': Some components of f0 = OdeFcn are NAN'])       
end
if ~isempty(MassFcn) 
  switch MassArgNbr
    case 0      
      MassFcnArg = {MassFcn};    
    case 1
      MassFcnArg = {MassFcn,t+h};
    case 2
      MassFcnArg = {MassFcn,t+h,y1};
    otherwise
      MassFcnArg = {MassFcn,t+h,y1,varargin{:}};
  end
  Mass = feval(MassFcnArg{:});
  if any(isnan(Mass))
    error([Solver_Name, ': Some components of Mass are NAN'])       
  end
  f1 = Mass \ f1;  
end

% ---- ESTIMATE THE SECOND DERIVATIVE OF THE SOLUTION
Sk   = AbsTol + RelTol .* abs(y);
Der2 = sum ( ((f1-f0)./Sk).^2 );   
Der2 = sqrt(Der2)/h;
% ---- STEP SIZE IS COMPUTED SUCH THAT
% ----  H**IORD * MAX ( NORM (F0), NORM (DER2)) = 0.01
Der12 = max(abs(Der2),sqrt(Dnf));
if Der12 <= 1e-15 
  h1 = max(1e-6,abs(h)*1e-3);
else
  h1 = (0.01/Der12)^(1/Ordre);
end
h = min([100*abs(h),h1,hmax])*PosNeg;
return
% -------------------------------------------------------------------
% -------------------------------------------------------------------


function [tout,yout,iout,Stop] = EventZeroFcn(EventFcn,t,y,T0,h, ...
                                 ycont,Flag,varargin)
% Cette fonction calcule la position d'un z�ro sur l'intervalle [tOld,tNew]
% connaissant yOld, yNew et ayant � disposition yCont qui permet le calcul
% de y(t) sur tout l'intervalle [tOld,tNew]. (1 zero par intervalle)
% La m�thode utilis�e est la Regula Falsi.
% La fonction est appel�e apr�s chaque nouvelle �valuation de (t,y) par
% la fonction dop54d. Elle utilise yCont, qui est sp�cifique � dop54d.

persistent t1  E1v

tout = [];
yout = [];
iout = [];
Stop = 0;

EventFcnArg = abs(nargin(EventFcn));
EventFcnVar = {EventFcn,t};
if EventFcnArg > 2
  EventFcnVar = {EventFcnVar{:},y,varargin{:}};
elseif EventFcnArg == 2
  EventFcnVar = {EventFcnVar{:},y};
end
  
if strcmp(Flag,'init')
  [E1v,Stopv,Slopev] = feval(EventFcnVar{:});   
  t1 = t; 
  return
end

[E2v,Stopv,Slopev] = feval(EventFcnVar{:});

t2 = t;
ky = length(E1v);

% Y a-t-il un z�ro sur l'intervalle?

IterMax = 50;
tAbsTol = 1e-6;
tRelTol = 1e-6;
EAbsTol = 1e-6;

Indk = 0;
for k = 1: ky 
  E1    = E1v(k);
  E2    = E2v(k);
  E12   = E1*E2;
  p12   = (E2-E1)/(t2-t1);
    
  if (E12 < 0) && (p12*Slopev(k) >= 0)
    
    Indk = Indk + 1;   % Indice de stockage   
    Done = false;
    Iter = 0;
    
    tNew = t2; % Pour �viter les ennuis dus � de tr�s petits intervalles
    yNew = y;
    ENew = E2;   
    while ~Done      
      Iter = Iter + 1;
      if Iter >= IterMax 
        fprintf(1,'EventZeroFcn: iteration number > maximal iteration number \n')    
        break
      end
      tRel = abs(t1-t2)*tRelTol < max(abs(t1),abs(t2));
      tAbs = abs(t1-t2) < tAbsTol;
      if (abs(ENew) < EAbsTol) && (tRel & tAbs)  % On a trouv�        
        break
      else    
        % Regula falsi  si possible, sinon dichotomie
        if abs(E1) < 200*EAbsTol || abs(E2) < 100*EAbsTol
          tNew = 0.5*(t1+t2);          
        else          
          tNew = (t1*E2-t2*E1)/(E2-E1);
        end
        S    = (tNew-T0)/h;
        S1   = 1-S;
        yNew = ycont(:,1) + S*(ycont(:,2) + S1*(ycont(:,3) + ...
                            S*(ycont(:,4) + S1*ycont(:,5))));
                                                                                                       
        EventFcnVar = {EventFcn,tNew};
        if EventFcnArg > 2
          EventFcnVar = {EventFcnVar{:},yNew,varargin{:}};
        elseif EventFcnArg == 2
          EventFcnVar = {EventFcnVar{:},yNew};
        end                        
        ENew = feval(EventFcnVar{:});                        
       
        ENew = ENew(k);
        if ENew * E1 > 0 
          t1 = tNew;
          E1 = ENew;
        else
          t2 = tNew;
          E2 = ENew;
        end
      end
    end
    ioutk(Indk)   = k;   
    toutk(Indk)   = tNew;
    youtk(Indk,:) = yNew;
    Stopk(Indk)   = Stopv(k);
  end
end

if exist('toutk')
  if t1 < t2
    [mt,Ind] = min(toutk);
  else
    [mt,Ind] = max(toutk);
  end
  iout = ioutk(Ind(1));
  tout = mt(1);
  yout = youtk(Ind(1),:);
  Stop = Stopk(Ind(1));
end
t1  = t2;
E1v = E2v;
% -------------------------------------------------------------------
% -------------------------------------------------------------------

function yLag = dop54d_Lag(PhiFcn,ny,t,y,varargin)

global tcont54Vec ycont54Vec hcont54Vec
% Search for the interval which contains t  


IndMax = length(tcont54Vec);
if IndMax == 0
  yLag = feval(PhiFcn,t,y,varargin{:}); 
  yLag = yLag(ny);
else
  Ind    = 1;
  while Ind <= IndMax && t > tcont54Vec(Ind) 
    Ind = Ind + 1;
  end
  if Ind > 1  
    Ind = Ind -1;
  end
  if t < tcont54Vec(Ind)  % --> phi
    yLag = feval(PhiFcn,t,y,varargin{:});
    yLag = yLag(ny);
  else
    t1  = tcont54Vec(Ind);
    h   = hcont54Vec(Ind);
    S   = (t-t1)/h;    % S is theta in the book
    S1  = 1-S;
    yLag  = ycont54Vec(Ind,ny,1) + S*(ycont54Vec(Ind,ny,2) + ...
            S1*(ycont54Vec(Ind,ny,3) + S*(ycont54Vec(Ind,ny,4) + ... 
            S1*ycont54Vec(Ind,ny,5))));        
  end
end
