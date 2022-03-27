#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
___________                   _________            .___      
\__    ___/___   _____   ____ \_   ___ \  ____   __| _/____  
  |    | /  _ \ /     \ /  _ \/    \  \/ /  _ \ / __ |/ __ \ 
  |    |(  <_> )  Y Y  (  <_> )     \___(  <_> ) /_/ \  ___/ 
  |____| \____/|__|_|  /\____/ \______  /\____/\____ |\___  >
                     \/               \/            \/    \/ 
   
    Author: Tomoki Koike
    Contact: tkoike3@gatech.edu
    Last Edited: 03-27-2022
    Description: ODE solvers for numerical integration with following methods
    (1) Euler (with half and full step)
    (2) RK4: Runge-Kutta 4th order (fixedstep)
    (3) RALS4: Ralston's (Runge-Kutta) 4th order (fixedstep)
    (4) RK42: Runge-Kutta 4th order (with half and full step)
    (5) RKF45: Embedded Runge-Kutta Fehlberg 4(5) 
    (6) CKRK45: Embedded Cash-Karp 4(5) 
    (7) DVERK65: Embedded Verner's 6(5) (unstable)
    (8) RKF78: Embeded Runge-Kutta Fehlberg 7(8) (unstable)

    References:
        [1] Solving Ordinary Differential Equations I, vol. 8. Berlin, 
        Heidelberg: Springer Berlin Heidelberg, 1993. 
        doi: 10.1007/978-3-540-78862-1.

        [2] T. Ritschel, “Numerical Methods For Solution of Differential 
        Equations,” p. 224.

"""


import inspect
import numpy as np
from numpy.typing import NDArray
from typing import Callable, Optional, Any, Union, List
from ode_solver import *
from stepsize_controllers import StepSizeControl, Ho110SSC


ODE = Union[
    Callable[[float, NDArray[Any], Optional[Any]],
             Union[NDArray[Any], NDArray[np.float64]]],
    Callable[[float, NDArray[Any]],
             Union[NDArray[Any], NDArray[np.float64]]]
]


METHOD = {
    "CKRK45": CKRK45,
    "DOP54": DOP54,
    "DVERK65": DVERK65,
    "Euler": Euler,
    "RALS4": RALS4,
    "RK4": RK4,
    "RK42": RK42,
    "RKF45": RKF45,
    "RKF78": RKF78
}

ORDER = {
    "CKRK45": 4,
    "DOP54": 5,
    "DVERK65": 6,
    "Euler": 1,
    "RALS4": 4,
    "RK4": 4,
    "RK42": 4,
    "RKF45": 4,
    "RKF78": 11
}

SSC = [
    "H211SSC",
    "H321SSC",
    "Ho110SSC",
    "PISSC",
    "PIDSSC",
    "PPIDSSC"
]


class NumInt:

    def __init__(self, derivative: ODE, h: float, 
                 tspan: Union[NDArray[Any], List[Any]],
                 IC: Union[NDArray[Any], List[Any]],
                 Atol: float = 1e-8, Rtol: float = 1e-5,
                 method: str = "RKF45", stepsizeControl: str = "PISSC"):
        self.h = h
        self.tspan = tspan
        self.IC = IC
        self.Atol = Atol
        self.Rtol = Rtol
        
        # Fixed step or adaptive step size depending on time span provided
        self.fixedstepsize = len(tspan) > 2
        if self.fixedstepsize:
            self.h = (tspan[-1] - tspan[0]) / (len(tspan) - 1)

        # Euler can only handle 1st order ODEs
        if len(IC) > 1 and method == "Euler":
            raise ValueError("`Euler` is only for 1st order ODEs.")

        # Check solver method and stepsize controller
        if method not in METHOD or not inspect.isclass(METHOD[method]):
            raise ValueError("`method` must be one of {}".format(METHOD))
        else:
            if method == "Euler":
                self.solver = METHOD[method](derivative)
            else:
                self.solver = METHOD[method](derivative, IC)
            self.p = ORDER[method]
            self.embedded = True
            # If the method is either `RK4` or `RALS4` override the fixedstepsize
            if method == "RK4" or method == "RALS4":
                self.fixedstepsize = True
                self.embedded = False
        
        if stepsizeControl not in SSC:
            raise ValueError("`stepsizeControl` must be one of {}".format(SSC))
        else:
            self.ssc = StepSizeControl(stepsizeControl, self.p)

        # Data holder for time and solution
        self.time = [tspan[0]]
        self.soln = [IC]

        # Store normalized errors and stepsizes 
        self.Es = [1.0, 1.0, 1.0]
        self.Hs = [1.0, 1.0, h]

        # States and time
        self.x = IC
        self.t = tspan[0]

        # Flag to check if it is the second step
        self.firststep = True
        self.secondstep = False


    
    def __call__(self, *args: Optional[Any]) -> bool:
        """Integration with stepsize control.

        Args:
             (Optional[Any]): Any optional parameters for ODE.

        Returns:
            (bool): step again flag -> True or False.
        """
        # Make sure the end point is included
        dt = self.tspan[-1] - self.t
        if self.h > dt:
            self.h = dt

        if self.fixedstepsize:  # fixed stepsize
            if self.embedded:
                xnp1, _ = self.solver(self.t, self.h, self.x, *args)
            else:
                xnp1 = self.solver(self.t, self.h, self.x, *args)

            # Update time and states
            self.t += self.h
            self.x = xnp1
            
            # Append values to the data storage
            self.soln.append(xnp1)
            self.time.append(self.t)
            return False  # Do not have to rerun step since fixed stepsize
        else:  # adaptive stepsize
            # Run the numerical integration 
            xnp1, err = self.solver(self.t, self.h, self.x, *args)

            # Normalized error estimate
            Enp1 = max(
                1e-10,
                np.max(np.abs(err) / (self.Atol + np.abs(xnp1)*self.Rtol))
            )

            if Enp1 <= 1:  # accpetable
                # Update normalized errors
                self.Es[0] = self.Es[1]
                self.Es[1] = self.Es[2]
                self.Es[2] = Enp1
                
                # Update states and time
                self.x = xnp1
                self.t += self.h

                # Append values to the data storage
                self.soln.append(xnp1)
                self.time.append(self.t)
                
                # Update stepsizes
                if self.firststep or self.secondstep:
                    # Do this for the first two steps
                    eta = Ho110SSC(Enp1, self.p)

                    self.firststep = False
                    self.secondstep = self.firststep or not self.secondstep
                else:
                    eta = self.ssc(self.Es, self.Hs)

                # Update stepsize 
                self.h *= eta
                self.Hs[0] = self.Hs[1]
                self.Hs[1] = self.Hs[2]
                self.Hs[2] = self.h
                return False  # Do not have to rerun step if E <= 1
            else:  # unacceptable
                eta = Ho110SSC(Enp1, self.p)
                self.h *= eta
                return True  # Rerun step with different stepsize since E > 1



if __name__=="__main__":

    # Test ------------------------------------------------------------------------
    def test1(t: float, y: NDArray[np.float64]) -> NDArray[Any]:
        return np.array([
            y[1],
            -4*t*y[1] - (2 + 4*t**2)*y[0]
        ])

    def test2(t: float, y: NDArray[np.float64]) -> NDArray[Any]:
        return -y[0] * np.cos(t)

    def test3(t: float, y: NDArray[np.float64]) -> NDArray[Any]:
        return np.array([
            (np.cos(t) - np.sin(t)*y[0])/y[1],
            np.sin(t)
        ])

    def test_exactsol1(t):
        # IC: y(0)=1, ydot(0)=0.5
        return np.array([
            0.5*np.exp(-t**2)*(t+2),
            (-t**2 - 2*t + 0.5)*np.exp(-t**2)
        ])

    def test_exactsol2(t):
        # IC: 0.5
        return 0.5 * np.exp(-np.sin(t))

    def test_exactsol3(t):
        # IC: y(0)=2, ydot(0)=1
        return np.array([
            (np.sin(t)+2)/(-np.cos(t)+2),
            -np.cos(t) + 2
        ])

    import matplotlib.pyplot as plt
    import time
    
    # Setup values 
    tf = 10  # final time
    t = 0   # initial time
    ti = t
    x = np.array([1, 0.5])  # initial conditions of states
    h = 0.001  # initial stepsize
    p = 4  # order of the numerical integration
    atol = 1e-10  # absolute tolerance
    rtol = 1e-8  # relative tolerance

    # Store data
    exsol = []  # exact solutions
    numint = NumInt(test1, h, [ti, tf], x, atol, rtol,
                    method="RKF45", stepsizeControl="PISSC")
    # numint = NumInt(test3, h, np.linspace(ti,tf,1000), x, atol, rtol, method="DOP54")


    # Numerical integration loop
    start = time.perf_counter()
    while numint.t < tf:

        exsol.append(test_exactsol1(numint.t))

        # Step again flag
        step_again = True
        
        while step_again:
            step_again = numint()
    end = time.perf_counter()
    
    # Append values for final time to exact solution and time
    exsol.append(test_exactsol1(tf))
    
    # Convert lists to np arrays 
    T = np.array(numint.time)
    simsol = np.array(numint.soln)
    exsol = np.array(exsol)

    print("Total data points: ", len(T)) 
    print(f"Total time elapsed: {end - start:0.6f} seconds")

    # Plot results
    plt.scatter(T, simsol[:,0], label="sim", c='r')
    plt.plot(T, exsol[:,0], label="exact")
    plt.plot(T, np.abs(simsol[:,0] - exsol[:,0]), label="error")
    plt.xlabel("t")
    plt.ylabel("x1")
    plt.legend()
    plt.grid(True, which='major', linestyle='--')
    plt.minorticks_on()
    plt.grid(True, which='minor', linestyle=':')
    plt.show()

    plt.scatter(T, simsol[:,1], label="sim", c='r')
    plt.plot(T, exsol[:,1], label="exact")
    plt.plot(T, np.abs(simsol[:,1] - exsol[:,1]), label="error")
    plt.xlabel("t")
    plt.ylabel("x2")
    plt.legend()
    plt.grid(True, which='major', linestyle='--')
    plt.minorticks_on()
    plt.grid(True, which='minor', linestyle=':')
    plt.show()

