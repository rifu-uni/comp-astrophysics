from typing import Callable, Tuple, Optional

RootStepFn = Callable[
    [Callable[[float], float], Optional[Callable[[float], float]], tuple],
    tuple,
]

def find_root(
    f: Callable[[float], float],
    df: Optional[Callable[[float], float]],
    state0: tuple,
    step_fn: RootStepFn,
    tol: float = 1e-8,
    max_iter: int = 200,
) -> float:
    """Generic scalar root-finder."""
    state = state0
    for _ in range(max_iter):
        if len(state) == 1:
            x = state[0]
        elif len(state) == 2:
            x = state[-1] if 'secant' in step_fn.__name__ else 0.5 * (state[0] + state[1])
        else:
            x = state[-1]

        if abs(f(x)) < tol:
            return x

        if len(state) == 2 and 'bisection' in step_fn.__name__ and abs(state[1] - state[0]) < tol:
            return x

        state = step_fn(f, df, state)

    raise RuntimeError(
        f"find_root did not converge after {max_iter} iterations "
        f"(last |f(x)| = {abs(f(x)):.3e})"
    )

def newton_step(
    f: Callable[[float], float],
    df: Callable[[float], float],
    state: Tuple[float],
) -> Tuple[float]:
    """Newton update."""
    x = state[0]
    dfx = df(x)
    if dfx == 0.0:
        raise ZeroDivisionError(f'Derivative is zero at x = {x}')
    return (x - f(x) / dfx,)

def bisection_step(
    f: Callable[[float], float],
    df: Optional[Callable[[float], float]],
    state: Tuple[float, float],
) -> Tuple[float, float]:
    """Bisection update."""
    a, b = state
    m = 0.5 * (a + b)
    return (a, m) if f(a) * f(m) <= 0.0 else (m, b)

def secant_step(
    f: Callable[[float], float],
    df: Optional[Callable[[float], float]],
    state: Tuple[float, float],
) -> Tuple[float, float]:
    """Secant method update."""
    x_prev, x_curr = state
    f_curr = f(x_curr)
    f_prev = f(x_prev)
    if f_curr - f_prev == 0.0:
        raise ZeroDivisionError(f"Denominator is zero in secant step at x={x_curr}")
    x_next = x_curr - f_curr * (x_curr - x_prev) / (f_curr - f_prev)
    return (x_curr, x_next)
