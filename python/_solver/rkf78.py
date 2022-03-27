import numpy as np
from numpy.typing import NDArray 
from typing import Callable, Optional, Any


ODE = Callable[[float, NDArray[np.float64],
                Optional[Any]], NDArray[np.float64]]


def rkf78(derivative: ODE, t: float, h: float, x: NDArray[np.float64],
        *params: Optional[Any]) -> tuple[float, NDArray[np.float64]]:
    """Runge-Kutta-Fehlburn 7(8)th order (embedded) method for more precision.

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
    s = 13
    A = np.zeros((s, s))
    B = np.zeros(s)
    Bhat = np.zeros(s)
    C = np.zeros(s)

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

    B[0,0] = 41/840
    B[0,5] = 34/105
    B[0,6] = 9/35
    B[0,7] = 9/35
    B[0,8] = 9/280
    B[0,9] = 9/280
    B[0,10] = 41/840

    Bhat[0,5] = 34/105
    Bhat[0,6] = 9/35
    Bhat[0,7] = 9/35
    Bhat[0,8] = 9/280
    Bhat[0,9] = 9/280
    Bhat[0,11] = 41/840
    Bhat[0,12] = 41/840

    D = B - Bhat
    K = np.zeros((len(x), s))
    
    # Each stages
    for i in range(s+1):
        K[:,i] = derivative(t + h*C[i], x + h*np.dot(A[i], K.T), *params).reshape(-1)

    xnew = x + h * np.dot(B, K.T)
    err = h * np.dot(D, K.T)

    return xnew.reshape(-1), err


