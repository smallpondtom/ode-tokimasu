import numpy as np
from numpy.typing import NDArray
from typing import Callable, Optional, Any


ODE = Callable[[float, NDArray[np.float64],
                Optional[Any]], NDArray[np.float64]]


def dverk65(derivative: ODE, t: float, h: float, x: NDArray[np.float64],
            *params: Optional[Any]) -> tuple[float, NDArray[np.float64]]:
    """Verner's Runge-Kutta 6(5)th order (embedded) method. A higher order 
    solution.

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
    D = B - Bhat
    K = np.zeros((len(x), s+1))

    # Each stages
    for i in range(s+1):
        K[:, i] = derivative(
            t + h*C[i], x + h*np.dot(A[i], K.T), *params).reshape(-1)

    xnew = K[:, s+1]
    err = h * np.dot(D, K.T)

    return xnew.reshape(-1), err
