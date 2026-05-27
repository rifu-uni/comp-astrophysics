"""Tridiagonal matrix solver using the Thomas algorithm (TDMA)."""

import numpy as np


def thomas_solve(
    a: np.ndarray,
    b: np.ndarray,
    c: np.ndarray,
    d: np.ndarray,
) -> np.ndarray:
    """
    Solve a tridiagonal linear system using the Thomas algorithm.

    The system is of the form::

        b[0] x[0] + c[0] x[1]                           = d[0]
        a[i-1] x[i-1] + b[i] x[i] + c[i] x[i+1]        = d[i]
        a[N-2] x[N-2] + b[N-1] x[N-1]                   = d[N-1]

    where the middle equation holds for ``1 <= i <= N-2``.

    Parameters
    ----------
    a : ndarray, shape (N-1,)
        Sub-diagonal entries, ``a[i] = A[i+1, i]``.
    b : ndarray, shape (N,)
        Main diagonal entries.
    c : ndarray, shape (N-1,)
        Super-diagonal entries, ``c[i] = A[i, i+1]``.
    d : ndarray, shape (N,)
        Right-hand side vector.

    Returns
    -------
    x : ndarray, shape (N,)
        Solution vector.

    Raises
    ------
    ValueError
        If the matrix is singular (near-zero pivot encountered).
    """
    n = len(d)

    if n == 1:
        if abs(b[0]) < 1e-15:
            raise ValueError("Thomas algorithm: singular matrix (b[0] ≈ 0)")
        return np.array([d[0] / b[0]])

    if abs(b[0]) < 1e-15:
        raise ValueError("Thomas algorithm: singular matrix (b[0] ≈ 0)")

    c_prime = np.zeros(n - 1)
    d_prime = np.zeros(n)

    c_prime[0] = c[0] / b[0]
    d_prime[0] = d[0] / b[0]

    for i in range(1, n):
        pivot = b[i] - a[i - 1] * c_prime[i - 1]
        if abs(pivot) < 1e-15:
            raise ValueError(f"Thomas algorithm: singular matrix at row i={i}")
        inv_pivot = 1.0 / pivot
        if i < n - 1:
            c_prime[i] = c[i] * inv_pivot
        d_prime[i] = (d[i] - a[i - 1] * d_prime[i - 1]) * inv_pivot

    x = np.zeros(n)
    x[-1] = d_prime[-1]
    for i in range(n - 2, -1, -1):
        x[i] = d_prime[i] - c_prime[i] * x[i + 1]

    return x
