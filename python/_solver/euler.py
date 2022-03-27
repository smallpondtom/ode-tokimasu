import numpy as np
from numpy.typing import NDArray
from typing import Callable, Optional, Any


ODE = Callable[[float, NDArray[np.float64],
                Optional[Any]], NDArray[np.float64]]


def euler(derivative: ODE, t: float, h: float, x: NDArray[np.float64],
          *params: Optional[Any]) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
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
    fn = derivative(t, x, *params).reshape(-1)
    xfull = x + h * fn

    # Double step
    xhalf = x + (h/2) * fn
    fn = derivative(t + h/2, xhalf, *params).reshape(-1)
    xdouble = xhalf + (h/2) * fn

    # Error estimate
    err = xfull - xdouble

    return xdouble, err
