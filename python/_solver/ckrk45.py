import numpy as np
from numpy.typing import NDArray
from typing import Callable, Optional, Any


ODE = Callable[[float, NDArray[np.float64],
                Optional[Any]], NDArray[np.float64]]


def ckrk45(derivative: ODE, t: float, h: float, x: NDArray[np.float64],
           *params: Optional[Any]) -> tuple[float, NDArray[np.float64]]:
    """Cash-Karp Runge-Kutta 4(5)th order (embedded) method.

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
        [1/5, 0, 0, 0, 0, 0],
        [3/40, 9/40, 0, 0, 0, 0],
        [3/10, -9/10, 6/5, 0, 0, 0],
        [-11/54, 5/2, -70/27, 35/27, 0, 0],
        [1631/55296, 175/512, 575/13824, 44275/110592, 253/4096, 0]
    ])
    B = np.array([37/378, 0, 250/621, 125/594, 0, 512/1771])
    Bhat = np.array([2825/27648, 0, 18575/48384, 13525/55296, 277/14336, 1/4])
    C = np.array([0, 1/5, 3/10, 3/5, 1, 7/8])
    D = B - Bhat
    K = np.zeros((len(x), s+1))

    # Each stages
    for i in range(s+1):
        K[:, i] = derivative(
            t + h*C[i], x + h*np.dot(A[i], K.T), *params).reshape(-1)

    xnew = x + h * np.dot(B, K.T)
    err = h * np.dot(D, K.T)

    return xnew.reshape(-1), err
