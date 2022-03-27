import numba as nb


@nb.njit
def pissc(err1: float, err2: float, p: int, epsilon: float = 0.8) -> float:
    """PI controller based adaptive stepsize controller (SSC).

    Args:
        err1 (float): normalized error between (n-1)th and (n)th states
        err2 (float): normalized error between (n)th and (n+1)th states
        p (int): order of the integration algorithm
        epsilon (float): tolerance [default 0.8]

    Returns:
        (float): optimized stepsize scale
    """
    sf = 0.9  # safety factor

    # Error controller coefficients
    ki = 0.4/(p+1)
    kp = 0.3/(p+1)

    eta = sf * (epsilon / err2)**ki * (err1 / err2)**kp
    return eta
