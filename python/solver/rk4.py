import numpy as np
from numpy.typing import NDArray 
from typing import Optional, Any


def rk4(derivative, t: float, h: float, x: NDArray[np.float64],
        *params: Optional[Any]) -> tuple[float, NDArray[np.float64]]:
    """Classic Runge-Kutta 4th order method.

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
        [1/2, 0, 0, 0],
        [0, 1/2, 0, 0],
        [0, 0, 1, 0]
    ])
    B = np.array([1/6, 1/3, 1/3, 1/6])
    C = np.array([0, 1/2, 1/2, 1])
    K = np.zeros((len(x), p))
    
    # Each stages
    for i in range(p):
        K[:,i] = derivative(t + h*C[i], x + h*np.dot(A[i], K.T), *params).reshape(-1)

    tnew = t + h
    xnew = x + h * np.dot(B, K.T)

    return tnew, xnew.reshape(-1)

