import numpy as np
from numpy.typing import NDArray
from typing import Callable, Optional, Any


ODE = Callable[[float, NDArray[np.float64],
                Optional[Any]], NDArray[np.float64]]


def rkf45(derivative: ODE, t: float, h: float, x: NDArray[np.float64],
          *params: Optional[Any]) -> tuple[float, NDArray[np.float64]]:
    """Runge-Kutta-Fehlburn 4(5)th order (embedded) method.

    Args:
        derivative (function): ODE function. Should return NDArray of states.
        t (float): time
        h (float): stepsize
        x (np.ndarray[float]): (n)th states
        *params (any): optional arguments.

    Returns:
        (np.ndarray[float]): (n+1)th states
        (np.ndarray[float]): (n+1)th error
    """

    # Butcher Tableau
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
    D = B - Bhat
    K = np.zeros((len(x), s+1))

    # Each stages
    for i in range(s+1):
        K[:, i] = derivative(
            t + h*C[i], x + h*np.dot(A[i], K.T), *params).reshape(-1)

    xnew = x + h * np.dot(B, K.T)
    err = h * np.dot(D, K.T)

    return xnew.reshape(-1), err
