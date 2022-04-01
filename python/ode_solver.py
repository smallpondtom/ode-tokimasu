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


import numpy as np
from numpy.typing import NDArray
from typing import Optional, Any, Union, Callable, Tuple, List


ODE = Union[
    Callable[[float, NDArray[Any], Optional[Any]],
             Union[NDArray[Any], NDArray[np.float64]]],
    Callable[[float, NDArray[Any]],
             Union[NDArray[Any], NDArray[np.float64]]],
    Callable[[float, NDArray[Any], NDArray[Any], Optional[Any]],
             Union[NDArray[Any], NDArray[np.float64]]],
    Callable[[float, NDArray[Any], NDArray[Any]],
             Union[NDArray[Any], NDArray[np.float64]]]
]

class OdeSolver:
    """Base ODE solver class.

    Attributes: 
        func: Derivative function.
        IC: Initial conditions.
        N: State dimension.
        s: Integrator order plus 1.
    """
    A: NDArray = NotImplemented
    B: NDArray = NotImplemented
    C: NDArray = NotImplemented
    D: NDArray = NotImplemented
    K: NDArray = NotImplemented
    Bhat: NDArray = NotImplemented
    p: int = NotImplemented


    def __init__(self, derivative: ODE, IC: Union[NDArray[Any], List[Any]]):
        self.func = derivative
        self.IC = np.array(IC).reshape(-1)
        self.N = len(IC)  # order of the ODE
        self.s = self.p + 1


    def butcher_tableau(self):
        # Print out the Butcher Tableau for verifications
        print("A:\n")
        print('\n'.join([''.join(['{:4}'.format(item) for item in row]) for row in self.A]))
        print("B:\n")
        print('\n'.join([''.join(['{:4}'.format(item) for item in row]) for row in self.B]))
        print("Bhat:\n")
        print('\n'.join([''.join(['{:4}'.format(item) for item in row]) for row in self.Bhat]))
        print("C:\n")
        print('\n'.join([''.join(['{:4}'.format(item) for item in row]) for row in self.C]))
        print("D:\n")
        print('\n'.join([''.join(['{:4}'.format(item) for item in row]) for row in self.D]))
        print("K:\n")
        print('\n'.join([''.join(['{:4}'.format(item) for item in row]) for row in self.K]))
   

    def solve(self, t: float, h: float, x: NDArray[Any],
              *params: Optional[Any]) -> NDArray[np.float64]:
        # Each stages
        for i in range(self.p):
            self.K[i,:] = self.func(t + h*self.C[i], x + h*np.dot(self.A[i], self.K), *params)

        xnew = x + h * np.dot(self.B, self.K)  

        return xnew.reshape(-1)


    def embedded_solve(self, t: float, h: float, x: NDArray[Any],
              *params: Optional[Any]) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
        # Each stages
        for i in range(self.s+1):
            self.K[i,:] = self.func(t + h*self.C[i], x + h*np.dot(self.A[i], self.K), *params)

        xnew = x + h * np.dot(self.B, self.K)  
        err = h * np.dot(self.D, self.K)   

        return xnew.reshape(-1), err


class Euler:

    def __init__(self, derivative: ODE):
        self.func = derivative

    def __call__(self, t: float, h: float, x: NDArray[Any],
                 *params: Optional[Any]) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Euler method for 1st order ODEs with full and half setp.

        Args:
            derivative (Callable): ODE function. Should return NDArray of states.
            t (float): time
            h (float): stepsize
            x (NDArray[np.float64]): (n)th states
            (Optional[Any]): optional arguments

        Returns:
            (NDArray[np.float64]): new states
            (NDArray[np.float64]): error
        """
        # Full step
        fn = self.func(t, x, *params).reshape(-1)
        xfull = x + h * fn

        # Double step
        xhalf = x + (h/2) * fn
        fn = self.func(t + h/2, xhalf, *params).reshape(-1)
        xdouble = xhalf + (h/2) * fn

        # Error estimate
        err = xfull - xdouble

        return xdouble, err
        

class RK4(OdeSolver):
    p = 4
    A = np.array([
        [0, 0, 0, 0],
        [1/2, 0, 0, 0],
        [0, 1/2, 0, 0],
        [0, 0, 1, 0]
    ])
    B = np.array([1/6, 1/3, 1/3, 1/6])     
    C = np.array([0, 1/2, 1/2, 1])
    
    def __init__(self, derivative: ODE, IC: Union[NDArray[Any], List[Any]]):
        super().__init__(derivative, IC)
        self.K = np.zeros((self.p, self.N))
    
    def __call__(self, t: float, h: float, x: NDArray[Any],
                 *params: Optional[Any]) -> NDArray[np.float64]:
        return self.solve(t, h, x, *params)


class RALS4(OdeSolver):
    p = 4  # order of solver
    A = np.array([
        [0, 0, 0, 0],
        [0.4, 0, 0, 0],
        [0.29697761, 0.15875964, 0, 0],
        [0.21810040, -3.05096516, 3.83286476, 0]
    ])
    B = np.array([0.17476028, -0.55148066, 1.20553560, 0.17118478])
    C = np.array([0, 0.4, 0.45573725, 1])

    def __init__(self, derivative: ODE, IC: Union[NDArray[Any], List[Any]]):
        super().__init__(derivative, IC)
        self.K = np.zeros((self.p, self.N))

    def __call__(self, t: float, h: float, x: NDArray[Any],
                 *params: Optional[Any]) -> NDArray[np.float64]:
        return self.solve(t, h, x, *params)


class RK42(RK4):

    def __init__(self, derivative: ODE, IC: Union[NDArray[Any], List[Any]]):
        super().__init__(derivative, IC)
        self.K = np.zeros((self.p, self.N))
    
    def RK42_step(self, fn, t, h, x, *params):                                      
        # Each stages (full step)
        self.K[0,:] = fn
        for i in range(1,self.p):
            self.K[i,:] = self.func(t + h*self.C[i], x + h*np.dot(self.A[i], self.K), *params)

        xnew = x + h * np.dot(self.B, self.K)
        return xnew.reshape(-1)

    def __call__(self, t: float, h: float, x: NDArray[Any],
                 *params: Optional[Any]
                ) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Classic Runge-Kutta 4th order (explicit) method with full and half step.

        Args:                                                                           
            derivative (function): ODE function. Should return NDArray of states.       
            t (float): time                                                             
            h (float): stepsize                                                         
            x (np.ndarray[float]): (n)th states                                         
            *params (any): optional arguments.                                          
                                                                                        
        Returns:
            (np.ndarray[float]): (n+1)th states
            (NDArray[np.float64]): error
        """

        # Each stages (full step)
        fn = self.func(t, x, *params)
        xfull = self.RK42_step(fn, t, h, x, *params)

        # Each stages (half step)
        xhalf = self.RK42_step(fn, t, h/2, x, *params)
        fn = self.func(t + h/2, xhalf, *params)
        xdouble = self.RK42_step(fn, t + h/2, h/2, xhalf, *params)
        err = xfull - xdouble

        return xdouble.reshape(-1), err


class RKF45(OdeSolver):
    p = 4  # order of solver
    s = p+1
    A = np.array([
        [0, 0, 0, 0, 0, 0],
        [1/4, 0, 0, 0, 0, 0],
        [3/32, 9/32, 0, 0, 0, 0],
        [1932/2197, -7200/2197, 7296/2197, 0, 0, 0],
        [439/216, -8, 3680/513, -845/4104, 0, 0],
        [-8/27, 2, -3544/2565, 1859/4104, -11/40, 0]
    ])
    B = np.array([25/216, 0, 1408/2565, 2197/4104, -1/5, 0])
    Bhat = np.array([16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55])
    C = np.array([0, 1/4, 3/8, 12/13, 1, 1/2])
    D = np.subtract(B, Bhat)

    def __init__(self, derivative: ODE, IC: Union[NDArray[Any], List[Any]]):
        super().__init__(derivative, IC)
        self.K = np.zeros((self.s+1, self.N))
    
    def __call__(self, t: float, h: float, x: NDArray[Any],
                 *params: Optional[Any]
                 ) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
        return self.embedded_solve(t, h, x, *params)


class CKRK45(RKF45):
    # Butcher Tableau                                                                 
    p = 4  # order of solver                                                          
    s = p+1                                                                           
    A = np.array([                                                                    
        [0, 0, 0, 0, 0, 0],                                                           
        [1/5, 0, 0, 0, 0, 0],
        [3/40, 9/40, 0, 0, 0, 0],
        [3/10, -9/10, 6/5, 0, 0, 0],
        [-11/54, 5/2, -70/27, 35/27, 0, 0],
        [1631/55296, 175/512, 575/13824, 44275/110592, 253/4096, 0]
    ])
    B = np.array([37/378, 0, 250/621, 125/594, 0, 512/1771])
    Bhat = np.array([2825/27648, 0, 18575/48384, 13525/55296, 277/14336, 1/4])
    C = np.array([0, 1/5, 3/10, 3/5, 1, 7/8])
    D = np.subtract(B, Bhat)


class DOP54(RKF45):
    # Butcher Tableau
    p = 5  # order of solver
    s = p+1
    A = np.array([                                                                    
        [0, 0, 0, 0, 0, 0, 0],                                                        
        [1/5, 0, 0, 0, 0, 0, 0],                                                      
        [3/40, 9/40, 0, 0, 0, 0, 0],                                                  
        [44/45, -56/15, 32/9, 0, 0, 0, 0],                                            
        [19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0, 0],                     
        [9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0, 0],                  
        [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0]
    ])
    Bhat = np.array([5179/57600, 0, 7571/16695, 393 /
                 640, -92097/339200, 187/2100, 1/40])
    B = A[-1]
    C = np.array([0, 1/5, 3/10, 4/5, 8/9, 1, 1])
    D = np.subtract(B, Bhat)


class RKF78(RKF45):
    # Butcher Tableau
    p = 11
    s = p+1
    A = np.zeros((s+1, s+1))
    B = np.zeros(s+1)
    Bhat = np.zeros(s+1)
    C = np.zeros(s+1)

    A[1,0] = 2/27
    A[2,0] = 1/36                                                                   
    A[2,1] = 1/12                                                                   
    A[3,0] = 1/24                                                                   
    A[3,2] = 1/8                                                                    
    A[4,0] = 5/12                                                                   
    A[4,2] = -25/16                                                                 
    A[4,3] = 25/16                                                                  
    A[5,0] = 1/20
    A[5,3] = 1/4
    A[5,4] = 1/5
    A[6,0] = -25/108
    A[6,3] = 125/108
    A[6,4] = -65/27
    A[6,5] = 125/54
    A[7,0] = 31/300
    A[7,4] = 61/225
    A[7,5] = -2/9
    A[7,6] = 13/900
    A[8,0] = 2
    A[8,3] = -53/6
    A[8,4] = 704/45
    A[8,5] = -107/9
    A[8,6] = 67/90
    A[8,7] = 3
    A[9,0] = -91/108
    A[9,3] = 23/108
    A[9,4] = -976/135
    A[9,5] = 311/54
    A[9,6] = -19/60
    A[9,7] = 17/6
    A[9,8] = -1/12
    A[10,0] = 2383/4100
    A[10,3] = -341/164
    A[10,4] = 4496/1025
    A[10,5] = -301/82
    A[10,6] = 2133/4100
    A[10,7] = 45/82
    A[10,8] = 45/164
    A[10,9] = 18/41
    A[11,0] = 3/205
    A[11,5] = -6/41
    A[11,6] = -3/205
    A[11,7] = -3/41
    A[11,8] = 3/41
    A[11,9] = 6/41
    A[12,0] = -1777/4100
    A[12,3] = -341/164
    A[12,4] = 4496/1025
    A[12,5] = -289/82            
    A[12,6] = 2193/4100          
    A[12,7] = 51/82              
    A[12,8] = 33/164             
    A[12,9] = 19/41              
    A[12,11] = 1                 
                                 
    B[0] = 41/840
    B[5] = 34/105
    B[6] = 9/35
    B[7] = 9/35
    B[8] = 9/280
    B[9] = 9/280
    B[10] = 41/840

    Bhat[5] = 34/105
    Bhat[6] = 9/35
    Bhat[7] = 9/35
    Bhat[8] = 9/280
    Bhat[9] = 9/280
    Bhat[11] = 41/840
    Bhat[12] = 41/840

    D = np.subtract(B, Bhat)


class DVERK65(RKF45):
    # Butcher Tableau
    p = 6  # order of solver
    s = p+1
    A = np.array([
        [0, 0, 0, 0, 0, 0, 0, 0],
        [1/6, 0, 0, 0, 0, 0, 0, 0],
        [4/75, 16/75, 0, 0, 0, 0, 0, 0],
        [5/6, -8/3, 5/2, 0, 0, 0, 0, 0],
        [-165/64, 55/6, -425/64, 85/96, 0, 0, 0, 0],
        [12/5, -8, 4015/612, -11/36, 88/255, 0, 0, 0],
        [-8263/15000, 124/75, -643/680, -81/250, 2484/10625, 0, 0, 0],
        [3501/1720, -300/43, 297275/52632, -319/2322, 24068/84065, 0, 3850/26703, 0]
    ])
    B = np.array([3/40, 0, 875/2244, 23/72, 264/1955, 0, 125/11592, 43/616])
    Bhat = np.array([13/160, 0, 2375/5984, 5/16, 12/85, 3/44, 0, 0])
    C = np.array([0, 1/6, 4/15, 2/3, 5/6, 1, 1/15, 1])
    D = np.subtract(B, Bhat)
