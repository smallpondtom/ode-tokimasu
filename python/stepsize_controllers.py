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
    Description: Stepsize controller for numerical integration.

    References:
        [1] C. A. Kennedy and M. H. Carpenter, “Diagonally implicit Runge–Kutta 
        methods for stiff ODEs,” Applied Numerical Mathematics, vol. 146, 
        pp. 221–244, Dec. 2019, doi: 10.1016/j.apnum.2019.07.008.

        [2] T. Ritschel, “Numerical Methods For Solution of Differential 
        Equations,” p. 224.

"""


from typing import List
import numba as nb

METHODS = [
    "Ho110SSC",
    "PISSC",
    "H211SSC",
    "PIDSSC",
    "PPIDSSC",
    "H321SSC"
]


def ssc(E: List[float], H: List[float], a: float, b: float, c: float, d: float,
        e: float, epsilon: float = 0.8) -> float:
    sf = 0.95  # safety factor
    eta = sf * ((epsilon/E[2])**a * (E[1]/epsilon)**b * (epsilon/E[0])**c
           * (H[2]/H[1])**d * (H[1]/H[0])**e)

    # Set bounds for this factor of eta
    if eta < 0.01:
        eta = 0.01
    elif eta > 10:
        eta = 10

    return eta


class StepSizeControl:
    """Stepsize controller class.

    Attributes: 
        epsilon: 
    """
    def __init__(self, method: str, p: float, epsilon: float = 0.8):
        self.epsilon = epsilon

        if method not in METHODS:
            raise ValueError("`method` must be one of {}.".format(METHODS))
        else:
            self.method = method

        # Set the coefficients for each method
        self.a, self.b, self.c, self.d, self.e, self.f = 0, 0, 0, 0, 0, 0
        if method == "Ho110SSC":  # asymptotic method
            self.a = 1/(p+1)
        elif method == "PISSC":  # PI controller
            self.a = 0.7/(p+1)
            self.b = 0.3/(p+1)
        elif method == "H211SSC":  # H211 controller
            self.a = 1/4/(p+1)
            self.b = -1/4/(p+1)
            self.d = -1/4
        elif method == "PIDSSC":  # PID controller
            self.a = 1/18/(p+1)
            self.b = -1/9/(p+1)
            self.c = 1/18/(p+1)
        elif method == "PPIDSSC":  # PPID controller
            self.a = 6/20/(p+1)
            self.b = -1/20/(p+1)
            self.c = -5/20/(p+1)
            self.d = 1
        elif method == "H321SSC":  # H321 controller
            self.a = 1/3/(p+1)
            self.b = -1/18/(p+1)
            self.c = -5/18/(p+1)
            self.d = 5/6
            self.e = 1/6

    def __call__(self, E: List[float], H: List[float]):
        return ssc(E, H, self.a, self.b, self.c, self.d, self.e, self.epsilon)


@nb.njit
def Ho110SSC(err2: float, p: int, epsilon: float = 0.8) -> float:
    """Compute adaptive asymptotic stepsize controller (SSC).

    Args:
        err2 (float): normalized error between (n)th and (n+1)th states
        p (int): order of the integration algorithm
        epsilon (float): tolerance. Should be between [0.8, 0.9].

    Returns:
        (float): optimized stepsize scale
    """
    sf = 0.95  # safety factor
    eta = sf * (epsilon / err2)**(1/(p+1))

    return eta

