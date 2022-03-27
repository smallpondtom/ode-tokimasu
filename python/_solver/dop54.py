import numpy as np
from numpy.typing import NDArray
from typing import Callable, Optional, Any


ODE = Callable[[float, NDArray[np.float64],
                Optional[Any]], NDArray[np.float64]]


def dop54(derivative: ODE, t: float, h: float, x: NDArray[np.float64],
          *params: Optional[Any]) -> tuple[float, NDArray[np.float64]]:
    """Dormand-Prince 5(4)th order (embedded) method.

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
    B = np.array([5179/57600, 0, 7571/16695, 393 /
                 640, -92097/339200, 187/2100, 1/40])
    Bhat = np.array([35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0])
    C = np.array([0, 1/5, 3/10, 4/5, 8/9, 1, 1])
    D = B - Bhat
    K = np.zeros((len(x), s+1))

    # Each stages
    for i in range(s+1):
        K[:, i] = derivative(
            t + h*C[i], x + h*np.dot(A[i], K.T), *params).reshape(-1)

    xnew = K[:, s+1]
    err = h * np.dot(D, K.T)

    return xnew.reshape(-1), err
