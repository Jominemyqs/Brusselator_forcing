"""Regression tests for the kernel-generic probit GP used by replay robustness."""

import numpy as np

from probit_gp import laplace_posterior as original_laplace_posterior
from probit_gp import predict_proba as original_predict_proba
from probit_gp_kernels import SUPPORTED_KERNELS, ard_kernel, laplace_posterior, predict_latent_joint, predict_proba


def synthetic_data():
    X = np.array([[-1.0, -0.5], [-0.5, 0.3], [0.2, -0.4], [0.6, 0.4], [1.0, 0.8]])
    y = np.array([-1.0, -1.0, -1.0, 1.0, 1.0])
    X_star = np.array([[-0.75, 0.0], [0.0, 0.0], [0.75, 0.5]])
    return X, y, X_star


def test_each_kernel_is_symmetric_and_positive_semidefinite():
    X, _, _ = synthetic_data()
    for family in SUPPORTED_KERNELS:
        K = ard_kernel(X, X, np.array([0.7, 1.2]), 1.4, family)
        np.testing.assert_allclose(K, K.T, atol=1e-13)
        assert np.min(np.linalg.eigvalsh(K)) > -1e-11


def test_squared_exponential_matches_frozen_primary_implementation():
    X, y, X_star = synthetic_data()
    lengthscales = np.array([0.8, 1.1])
    signal_std = 1.7
    original = original_laplace_posterior(X, y, lengthscales, signal_std)
    extended = laplace_posterior(
        X,
        y,
        lengthscales,
        signal_std,
        kernel_family="squared_exponential",
    )
    np.testing.assert_allclose(extended.K, original.K, rtol=0.0, atol=1e-13)
    np.testing.assert_allclose(extended.f, original.f, rtol=1e-11, atol=1e-12)
    for left, right in zip(predict_proba(extended, X_star), original_predict_proba(original, X_star)):
        np.testing.assert_allclose(left, right, rtol=1e-11, atol=1e-12)


def test_predictions_and_joint_covariance_are_finite_for_all_kernels():
    X, y, X_star = synthetic_data()
    for family in SUPPORTED_KERNELS:
        state = laplace_posterior(
            X,
            y,
            np.array([0.8, 1.1]),
            1.7,
            kernel_family=family,
        )
        probability, mean, variance = predict_proba(state, X_star)
        joint_mean, covariance = predict_latent_joint(state, X_star)
        assert np.all(np.isfinite(probability))
        assert np.all((probability > 0.0) & (probability < 1.0))
        assert np.all(np.isfinite(mean))
        assert np.all(variance >= 0.0)
        np.testing.assert_allclose(joint_mean, mean, rtol=1e-11, atol=1e-12)
        np.testing.assert_allclose(covariance, covariance.T, atol=1e-13)
        assert np.min(np.linalg.eigvalsh(covariance)) > -1e-9
