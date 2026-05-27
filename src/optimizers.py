"""General-purpose optimisation methods for maximisation.

Provides three classic algorithms:
  - Golden-section search for unimodal scalar functions.
  - Newton's method for optimisation (root of the first derivative).
  - Gradient ascent for multivariate functions.
"""

import numpy as np
from typing import Callable, List, Tuple


def golden_section_max(
    f: Callable[[float], float],
    a: float,
    b: float,
    tol: float = 1e-6,
    max_iter: int = 100,
) -> float:
    """
    Maximise a unimodal scalar function *f* on the interval [*a*, *b*].

    Uses the golden-section search, which reduces the bracket size by the
    golden-ratio factor (≈ 0.618) at each iteration.  The function is assumed
    to be unimodal with a single maximum inside the given interval.

    Parameters
    ----------
    f : callable(x) -> float
        Unimodal function to maximise.
    a : float
        Left endpoint of the search interval.
    b : float
        Right endpoint of the search interval.  Must satisfy a < b.
    tol : float, optional
        Convergence tolerance on the bracket width (default 1e-6).
    max_iter : int, optional
        Maximum number of iterations (default 100).

    Returns
    -------
    x_opt : float
        Location of the maximum (midpoint of the final bracket).
    """
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    inv_phi = 1.0 / phi

    c = b - (b - a) * inv_phi
    d = a + (b - a) * inv_phi

    for _ in range(max_iter):
        if abs(b - a) < tol:
            return (a + b) / 2.0

        if f(c) > f(d):
            b = d
        else:
            a = c

        c = b - (b - a) * inv_phi
        d = a + (b - a) * inv_phi

    return (a + b) / 2.0


def newton_opt_max(
    f: Callable[[float], float],
    df: Callable[[float], float],
    d2f: Callable[[float], float],
    x0: float,
    tol: float = 1e-8,
    max_iter: int = 100,
) -> float:
    """
    Maximise *f* via Newton's method for optimisation.

    Finds a stationary point of *f* by iteratively solving  df(x) = 0  using
    the update  x_{k+1} = x_k - f'(x_k) / f''(x_k).  Convergence is quadratic
    near the optimum provided *f* is twice-differentiable and the starting
    guess *x0* is sufficiently close.

    Parameters
    ----------
    f : callable(x) -> float
        Objective function (required for the interface; not evaluated).
    df : callable(x) -> float
        First derivative of *f*.
    d2f : callable(x) -> float
        Second derivative of *f*.
    x0 : float
        Initial guess.
    tol : float, optional
        Convergence tolerance on |df(x)| and step size (default 1e-8).
    max_iter : int, optional
        Maximum number of iterations (default 100).

    Returns
    -------
    x_opt : float
        Estimated location of the maximum.

    Raises
    ------
    ValueError
        If |d2f(x)| < 1e-15, indicating the second derivative is near-zero
        and the Newton step would be ill-defined.
    """
    x = x0

    for _ in range(max_iter):
        fprime = df(x)
        fdouble = d2f(x)

        if abs(fprime) < tol:
            return x

        if abs(fdouble) < 1e-15:
            raise ValueError(f"Second derivative too small at x = {x}")

        x_new = x - fprime / fdouble

        if abs(x_new - x) < tol:
            return x_new

        x = x_new

    return x


def gradient_ascent(
    f: Callable[[np.ndarray], float],
    grad_f: Callable[[np.ndarray], np.ndarray],
    x0: np.ndarray,
    lr: float = 0.01,
    tol: float = 1e-6,
    max_iter: int = 1000,
) -> Tuple[np.ndarray, List[np.ndarray]]:
    """
    Maximise a multivariate function *f* using gradient ascent.

    Iteratively moves in the direction of the gradient:

        x_{k+1} = x_k + lr * ∇f(x_k)

    The iteration stops when ‖∇f(x_k)‖ < tol or *max_iter* is reached.

    Parameters
    ----------
    f : callable(x) -> float
        Objective function of a vector argument; required for the interface.
    grad_f : callable(x) -> ndarray
        Gradient of *f*; returns an ndarray of the same shape as *x*.
    x0 : ndarray
        Initial point.
    lr : float, optional
        Learning rate (step-size multiplier, default 0.01).
    tol : float, optional
        Convergence tolerance on the gradient norm (default 1e-6).
    max_iter : int, optional
        Maximum number of iterations (default 1000).

    Returns
    -------
    x_opt : ndarray
        Final point when the iteration terminated.
    history : list of ndarray
        Sequence of points visited, starting with *x0* and ending with
        the final iterate.
    """
    x = np.array(x0, dtype=float)
    history: List[np.ndarray] = [x.copy()]

    for _ in range(max_iter):
        g = grad_f(x)
        g_norm = np.linalg.norm(g)

        if g_norm < tol:
            break

        x = x + lr * g
        history.append(x.copy())

    return x, history
