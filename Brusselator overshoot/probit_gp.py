"""Small deterministic probit Gaussian-process classifier.

The implementation follows the standard Laplace approximation for binary GP
classification. It is intentionally limited to the small two-dimensional
boundary experiment and avoids introducing an unavailable ML dependency.
"""

from dataclasses import dataclass
from typing import Iterable

import numpy as np
from scipy.linalg import cholesky, cho_solve, solve_triangular
from scipy.optimize import minimize
from scipy.special import log_ndtr, ndtr


@dataclass
class LaplaceState:
    X: np.ndarray
    y: np.ndarray
    lengthscales: np.ndarray
    signal_std: float
    jitter: float
    K: np.ndarray
    f: np.ndarray
    a: np.ndarray
    sqrt_W: np.ndarray
    L: np.ndarray
    log_marginal_likelihood: float
    iterations: int
    converged: bool


def ard_rbf_kernel(
    X1: np.ndarray,
    X2: np.ndarray,
    lengthscales: np.ndarray,
    signal_std: float,
) -> np.ndarray:
    X1 = np.asarray(X1, dtype=float)
    X2 = np.asarray(X2, dtype=float)
    lengthscales = np.asarray(lengthscales, dtype=float)
    scaled = (X1[:, None, :] - X2[None, :, :]) / lengthscales
    return signal_std**2 * np.exp(-0.5 * np.sum(scaled**2, axis=2))


def _probit_terms(y: np.ndarray, f: np.ndarray):
    z = y * f
    log_likelihood = log_ndtr(z)
    log_density = -0.5 * z**2 - 0.5 * np.log(2.0 * np.pi)
    ratio = np.exp(np.clip(log_density - log_likelihood, -700.0, 50.0))
    gradient = y * ratio
    W = np.maximum(ratio * (ratio + z), 1e-14)
    return log_likelihood, gradient, W


def laplace_posterior(
    X: np.ndarray,
    y: np.ndarray,
    lengthscales: np.ndarray,
    signal_std: float,
    *,
    jitter: float = 1e-8,
    maximum_iterations: int = 100,
    tolerance: float = 1e-10,
) -> LaplaceState:
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    if X.ndim != 2 or y.shape != (X.shape[0],):
        raise ValueError("X and y dimensions are inconsistent")
    if not np.all(np.isin(y, [-1.0, 1.0])):
        raise ValueError("Probit labels must be encoded as -1 and +1")

    K = ard_rbf_kernel(X, X, lengthscales, signal_std)
    K = K + jitter * np.eye(X.shape[0])
    f = np.zeros(X.shape[0], dtype=float)
    converged = False
    iterations = 0
    for iterations in range(1, maximum_iterations + 1):
        _, gradient, W = _probit_terms(y, f)
        sqrt_W = np.sqrt(W)
        B = np.eye(X.shape[0]) + (sqrt_W[:, None] * K) * sqrt_W[None, :]
        L = cholesky(B, lower=True, check_finite=False)
        b = W * f + gradient
        Kb = K @ b
        correction = sqrt_W * cho_solve(
            (L, True), sqrt_W * Kb, check_finite=False
        )
        a = b - correction
        updated = K @ a
        if np.max(np.abs(updated - f)) < tolerance:
            f = updated
            converged = True
            break
        f = updated

    log_likelihood, gradient, W = _probit_terms(y, f)
    sqrt_W = np.sqrt(W)
    B = np.eye(X.shape[0]) + (sqrt_W[:, None] * K) * sqrt_W[None, :]
    L = cholesky(B, lower=True, check_finite=False)
    b = W * f + gradient
    Kb = K @ b
    a = b - sqrt_W * cho_solve((L, True), sqrt_W * Kb, check_finite=False)
    log_marginal = (
        -0.5 * float(a @ f)
        + float(np.sum(log_likelihood))
        - float(np.sum(np.log(np.diag(L))))
    )
    if not all(
        np.all(np.isfinite(value)) for value in [f, a, sqrt_W, L]
    ) or not np.isfinite(log_marginal):
        raise FloatingPointError("Nonfinite Laplace posterior state")
    return LaplaceState(
        X=X,
        y=y,
        lengthscales=np.asarray(lengthscales, dtype=float),
        signal_std=float(signal_std),
        jitter=float(jitter),
        K=K,
        f=f,
        a=a,
        sqrt_W=sqrt_W,
        L=L,
        log_marginal_likelihood=log_marginal,
        iterations=iterations,
        converged=converged,
    )


def predict_proba(state: LaplaceState, X_star: np.ndarray):
    X_star = np.asarray(X_star, dtype=float)
    K_star = ard_rbf_kernel(
        state.X, X_star, state.lengthscales, state.signal_std
    )
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        latent_mean = K_star.T @ state.a
    projected = solve_triangular(
        state.L,
        state.sqrt_W[:, None] * K_star,
        lower=True,
        check_finite=False,
    )
    prior_variance = np.full(X_star.shape[0], state.signal_std**2)
    latent_variance = np.maximum(
        prior_variance - np.sum(projected**2, axis=0), 0.0
    )
    probability_positive = ndtr(
        latent_mean / np.sqrt(1.0 + latent_variance)
    )
    if not all(
        np.all(np.isfinite(value))
        for value in [latent_mean, latent_variance, probability_positive]
    ):
        raise FloatingPointError("Nonfinite GP prediction")
    return probability_positive, latent_mean, latent_variance


def predict_latent_joint(state: LaplaceState, X_star: np.ndarray):
    """Return the Laplace posterior mean and joint covariance of latent f.

    Unlike :func:`predict_proba`, this retains cross-covariances between test
    states.  Those correlations are required when estimating whether two
    nearby states lie on opposite sides of a learned separator.
    """
    X_star = np.asarray(X_star, dtype=float)
    K_star = ard_rbf_kernel(
        state.X, X_star, state.lengthscales, state.signal_std
    )
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        latent_mean = K_star.T @ state.a
        projected = solve_triangular(
            state.L,
            state.sqrt_W[:, None] * K_star,
            lower=True,
            check_finite=False,
        )
        prior_covariance = ard_rbf_kernel(
            X_star, X_star, state.lengthscales, state.signal_std
        )
        covariance = prior_covariance - projected.T @ projected
    covariance = 0.5 * (covariance + covariance.T)
    diagonal = np.diag(covariance)
    if np.min(diagonal) < -1e-8:
        raise FloatingPointError("Joint latent covariance has a negative variance")
    covariance[np.diag_indices_from(covariance)] = np.maximum(diagonal, 0.0)
    if not np.all(np.isfinite(latent_mean)) or not np.all(np.isfinite(covariance)):
        raise FloatingPointError("Nonfinite joint GP prediction")
    return latent_mean, covariance


def optimize_hyperparameters(
    X: np.ndarray,
    y: np.ndarray,
    starts: Iterable[Iterable[float]],
    *,
    lengthscale_bounds=(0.12, 3.0),
    signal_std_bounds=(0.2, 10.0),
    jitter=1e-8,
    maximum_laplace_iterations=100,
    laplace_tolerance=1e-10,
    maximum_optimizer_iterations=250,
    function_tolerance=1e-10,
    log_prior_mean=None,
    log_prior_std=None,
):
    bounds = [
        (np.log(lengthscale_bounds[0]), np.log(lengthscale_bounds[1]))
        for _ in range(X.shape[1])
    ] + [(np.log(signal_std_bounds[0]), np.log(signal_std_bounds[1]))]

    if (log_prior_mean is None) != (log_prior_std is None):
        raise ValueError("log_prior_mean and log_prior_std must be supplied together")
    if log_prior_mean is not None:
        log_prior_mean = np.asarray(log_prior_mean, dtype=float)
        log_prior_std = np.asarray(log_prior_std, dtype=float)
        expected_shape = (X.shape[1] + 1,)
        if log_prior_mean.shape != expected_shape or log_prior_std.shape != expected_shape:
            raise ValueError("Log-parameter prior dimensions are inconsistent")
        if np.any(log_prior_std <= 0):
            raise ValueError("Log-parameter prior standard deviations must be positive")

    def objective(log_parameters):
        lengthscales = np.exp(log_parameters[:-1])
        signal_std = float(np.exp(log_parameters[-1]))
        try:
            state = laplace_posterior(
                X,
                y,
                lengthscales,
                signal_std,
                jitter=jitter,
                maximum_iterations=maximum_laplace_iterations,
                tolerance=laplace_tolerance,
            )
            if not np.isfinite(state.log_marginal_likelihood):
                return 1e30
            value = -state.log_marginal_likelihood
            if log_prior_mean is not None:
                standardized = (log_parameters - log_prior_mean) / log_prior_std
                value += 0.5 * float(standardized @ standardized)
            return value
        except (np.linalg.LinAlgError, FloatingPointError, ValueError):
            return 1e30

    optimization_records = []
    best = None
    for start_index, start in enumerate(starts, start=1):
        start = np.asarray(start, dtype=float)
        if start.shape != (X.shape[1] + 1,):
            raise ValueError("Each optimization start must contain one lengthscale per feature and one signal standard deviation")
        result = minimize(
            objective,
            np.log(start),
            method="L-BFGS-B",
            bounds=bounds,
            options={
                "maxiter": maximum_optimizer_iterations,
                "ftol": function_tolerance,
            },
        )
        record = {
            "start_index": start_index,
            "success": bool(result.success),
            "message": str(result.message),
            "objective": float(result.fun),
            "signal_std": float(np.exp(result.x[-1])),
            "iterations": int(result.nit),
            "function_evaluations": int(result.nfev),
        }
        for feature_index, value in enumerate(np.exp(result.x[:-1]), start=1):
            record[f"lengthscale_{feature_index}"] = float(value)
        optimization_records.append(record)
        if np.isfinite(result.fun) and (best is None or result.fun < best.fun):
            best = result
    if best is None:
        raise RuntimeError("Every GP hyperparameter optimization start failed")

    lengthscales = np.exp(best.x[:-1])
    signal_std = float(np.exp(best.x[-1]))
    state = laplace_posterior(
        X,
        y,
        lengthscales,
        signal_std,
        jitter=jitter,
        maximum_iterations=maximum_laplace_iterations,
        tolerance=laplace_tolerance,
    )
    return state, optimization_records
