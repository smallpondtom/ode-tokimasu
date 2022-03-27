import numba as nb


@nb.njit
def h321ssc(hn1: float, h0: float, h1: float, err0: float, err1: float,
             err2: float, p: int, epsilon: float = 0.8) -> float:
    """H321 adaptive stepsize controller (SSC).

    Args:
        hn1 (float): (n-2)th stepsize (h negative 1)
        h0 (float): (n-1)th stepsize
        h1 (float): (n)th stepsize
        err0 (float): normalized error between (n-2)th and (n-1)th states
        err1 (float): normalized error between (n-1)th and (n)th states
        err2 (float): normalized error between (n)th and (n+1)th states
        p (int): order of the integration algorithm
        epsilon (float): tolerance. Should be between [0.8, 0.9].

    Returns:
        (float): optimized stepsize scale
    """
    sf = 0.9  # safety factor

    # Error controller coefficients
    alpha = 1/3/p
    beta = -1/18/p
    gamma = -5/18/p
    delta = 5/6
    epsilon = 1/6

    eta = (sf * (epsilon / err2)**alpha * (err1 / epsilon) **
            beta * (epsilon / err0)**gamma * (h1 / h0)**delta * (h0 / hn1)**epsilon)
    return eta
