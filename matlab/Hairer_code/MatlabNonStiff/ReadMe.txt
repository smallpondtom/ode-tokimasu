dop54.m      Numerical solution of a non-stiff system of first order ordinary 
             differential equations:
             This is an explicit Runge-Kutta method of order (4)5 due to
             Dormand & Prince (with step size control and dens output).
			 Good for mid precision (1e-4 - 1e-7)
			            
dop853.m     Numerical solution of a non-stiff system of first order ordinary 
             differential equations:
             This is an explicit Runge-Kutta method of order (8)5,3 due to
             Dormand & Prince (with step size control and dens output).
             Good for high precision (1e-6 - 1e-12)

dop54d.m     Numerical solution of a non-stiff system of first order ordinary 
			 differential equations with delay.
			 This is an explicit Runge-Kutta method of order (4)5 due to
			 Dormand & Prince (with step size control and dense output).
			 Good for mid precision (1e-4 - 1e-7)
			 
dop853D.m    Numerical solution of a non-stiff system of first order ordinary 
             differential equations with delay.
             This is an explicit Runge-Kutta method of order (8) 5,3 due to
             Dormand & Prince (with step size control and dense output).			 
			 Good for high precision (1e-6 - 1e-12)
			 
			 AUTHORS OF THE FORTRAN CODE:
			 E. HAIRER AND G. WANNER
             UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
             CH-1211 GENEVE 24, SWITZERLAND 
             E-MAIL:  Ernst.Hairer@math.unige.ch
                      Gerhard.Wanner@math.unige.ch
			 
			 THESE CODES ARE DESCRIBED IN:
             E. HAIRER, S.P. NORSETT AND G. WANNER, SOLVING ORDINARY
             DIFFERENTIAL EQUATIONS I. NONSTIFF PROBLEMS. 2ND EDITION.
             SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
			 SPRINGER-VERLAG (1993)
			 
			 
radau.m      Numerical solution of a stiff (or differential algebraic) system of
             first order ordinary differential equations:
             Mass*y' = OdeFcn(t,y).
             The system can be (linearly) implicit (mass-matrix Mass ~= I)
             or explicit (Mass = I)
             The code is based on implicit Runge-Kutta methods (Radau IIa)
             with variable order (1, 5, 9, 13), with step size control and
             continuous output.

             AUTHORS: E. HAIRER AND G. WANNER
             UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
             CH-1211 GENEVE 24, SWITZERLAND 
             E-MAIL:  Ernst.Hairer@math.unige.ch
                      Gerhard.Wanner@math.unige.ch
					  
 			 THIS CODE IS DESCRIBED IN:                  
             E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
             EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
             SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS 14,
             SPRINGER-VERLAG 1991, SECOND EDITION 1996.

    Matlab version:
    Denis Bichsel
    Rue des Deurres 58
    2000 Neuchâtel
    Suisse
    dbichsel@infomaniak.ch
	Version of end 2015
