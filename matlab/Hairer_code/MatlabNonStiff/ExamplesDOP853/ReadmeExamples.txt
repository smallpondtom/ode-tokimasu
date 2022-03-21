% --------------------------------------------------------- %
% Readme.txt                                                %
%                                                           %
% List of the examples used to test dop54 and dop853        %
%                                                           %
% Denis Bichsel                                             %
% 58 Rue des Deurres                                        %
% 2000 Neuchâtel                                            %
% Tel. 41 (0) 32 730 10 16                                  %
% email: dbichsel@infomaniak.ch                             %
% --------------------------------------------------------- %
The functions dop54 and dop853 are translation in Matlab from
the Fortran code written by :
E. Hairer and G. Wanner
Université de Genève.
% --------------------------------------------------------- %
The options for dop54 and dop853 can be set with dopset
and read with dopget.
These function are very similar to Matalb function ode45, 
odeset and odeget.

The allowed options are essentially the same as in Matlab.
The options
  BDF
  Jacobian
  JConstant
  JPattern
  Refine
  Vectorized
  MvPattern
  MassSingular
  InitialSlope
  MaxOrder
are not implemented.

Remark : The following options are not implemented for ode45
in Matlab   
  BDF
  Jacobian
  JConstant
  JPattern
  Vectorized
  MvPattern
  InitialSlope
  MassSingular
  InitialSlope
  MaxOrder
  
The option "Refine" is just useless with dop54 and dop 853,
because it's possible to choose all the output value for 
the t variable without any effect on the step size calculation.

DenseOuputFcn,  DenseOutputSel.

Like the "OutputFcn" and "OutputSel" options, "DenseOutputFcn" 
and "DenseOutputSel" allow to get the t and y values after each 
step calculation. More, with this option the user gets also the
coefficients which allow dense calculation on the intervalle
[t, t+h], h = step size,  5 coefficients for dop54 and 
8 coefficients for dop853 (see example below).
These options allows for example, the calculation with high 
precision of the extrema of a function 
% --------------------------------------------------------- %


List of examples:
------------------

The name of the main function of the examples end always
with   "_dop54_dop853.m". The other functions with a name 
which begin on the same way, are used for the calculation.
For example :  Arenstorf_dop54_dop853.m uses the function 
"Arenstorf.m"

Arenstorf_dop54_dop853.m
  This function tests the following options:
    RelTol
    AbsTol
    MaxIter
    MaxStep
  and compare the execution time of dop54, dop853 and 
  Matlab ode45
	  
Ballode1_dop54_dop853.m	  
Ballode2_dop54_dop853.m
  These functions test the following options:
    Events
	OutputFcn
	OutputSel
	InitialStep
	MaxStep
  
Evol_Soleilplanete_dop54_dop853.m
  This function tests the following options:
    RelTol
    AbsTol
	InitialStep
	MaxStep
	Stats
	DenseOutputFcn
	DenseOutputSel
  	
  This function shows very clearly the difference of
  dop54 and dop853 in solving the Mercury's perihelion 
  problem. dop54 like ode45 don't succeed but dop853 result is
  531.88 arc second (true value is 531.54 in the Newton's
  theory)     

Mass_dop54_dop853.m
  This function tests the following options:
    Mass 
	MStateDependence

NonNegativ_dop54_dop853.m
  This function tests the following option:
    NonNegative

Osci_dop54_dop853.m
  This very simple function compare the time needed in solving
  the same problem (harmonic oscillator) by dop54, dop853 and
  the Matlab function ode45.

VanderPol_dop54_dop853.m
  This function tests the following options:
    OutputSel
	DenseOutputFcn
	DenseOutputSel
	InitialStep
	MaxStep
% --------------------------------------------------------- %