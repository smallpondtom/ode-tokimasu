import numpy as np
from numpy.typing import NDArray 
from typing import Callable, Optional, Any


ODE = Callable[[float, NDArray[np.float64],
                Optional[Any]], NDArray[np.float64]]


def rk42(derivative: ODE, t: float, h: float, x: NDArray[np.float64],
        *params: Optional[Any]) -> tuple[float, NDArray[np.float64]]:
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
    fn = derivative(t, x, *params)
    xfull = rk42_step(derivative, fn, t, h, x, *params)

    # Each stages (half step)
    xhalf = rk42_step(derivative, fn, t, h/2, x, *params)
    fn = derivative(t + h/2, xhalf, *params)
    xdouble = rk42_step(derivative, fn, t + h/2, h/2, xhalf, *params)

    err = xfull - xdouble

    return xdouble.reshape(-1), err


def rk42_step(func, fn, t, h, x, *params):
    # Butcher Tableau 
    p = 4  # order of solver
    A = np.array([
        [0, 0, 0, 0],
        [1/2, 0, 0, 0],
        [0, 1/2, 0, 0],
        [0, 0, 1, 0]
    ])
    B = np.array([1/6, 1/3, 1/3, 1/6])
    C = np.array([0, 1/2, 1/2, 1])
    K = np.zeros((len(x), p))

    # Each stages (full step)
    K[:,0] = fn.reshape(-1)
    for i in range(1,p):
        K[:,i] = func(t + h*C[i], x + h*np.dot(A[i], K.T), *params).reshape(-1)

    xnew = x + h * np.dot(B, K.T) 
    return xnew
