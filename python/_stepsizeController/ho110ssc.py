import numba as nb


@nb.njit
def ho110ssc(err2: float, p: int, epsilon: float = 0.8) -> float:
    """Compute adaptive asymptotic stepsize controller (SSC).

    Args:
        err2 (float): normalized error between (n)th and (n+1)th states
        p (int): order of the integration algorithm
        epsilon (float): tolerance. Should be between [0.8, 0.9].

    Returns:
        (float): optimized stepsize scale
    """
    sf = 0.9  # safety factor
    eta = sf * (epsilon / err2)**(1/(p+1))

    return eta

