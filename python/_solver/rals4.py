import numpy as np
from numpy.typing import NDArray
from typing import Callable, Optional, Any


ODE = Callable[[float, NDArray[np.float64],
                Optional[Any]], NDArray[np.float64]]


def rals4(derivative: ODE, t: float, h: float, x: NDArray[np.float64],
          *params: Optional[Any]) -> tuple[float, NDArray[np.float64]]:
    """Ralston's 4th order (explicit) method which has minimum truncation error.

    Args:
        derivative (function): ODE function. Should return NDArray of states.
        t (float): time
        h (float): stepsize
        x (np.ndarray[float]): (n)th states
        *params (any): optional arguments.

    Returns:
        (float): new time
        (np.ndarray[float]): (n+1)th states
    """

    # Butcher Tableau
    p = 4  # order of solver
    A = np.array([
        [0, 0, 0, 0],
        [0.4, 0, 0, 0],
        [0.29697761, 0.15875964, 0, 0],
        [0.21810040, -3.05096516, 3.83286476, 0]
    ])
    B = np.array([0.17476028, -0.55148066, 1.20553560, 0.17118478])
    C = np.array([0, 0.4, 0.45573725, 1])
    K = np.zeros((len(x), p))

    # Each stages
    for i in range(p):
        K[:, i] = derivative(
            t + h*C[i], x + h*np.dot(A[i], K.T), *params).reshape(-1)

    tnew = t + h
    xnew = x + h * np.dot(B, K.T)

    return tnew, xnew.reshape(-1)
