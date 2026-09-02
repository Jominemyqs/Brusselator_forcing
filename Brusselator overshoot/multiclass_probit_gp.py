"""Small one-vs-rest multiclass wrapper around the local probit-GP code."""

from dataclasses import dataclass

import numpy as np

from probit_gp import optimize_hyperparameters, predict_latent_joint, predict_proba


@dataclass
class OneVsRestModel:
    classes: tuple
    states: dict
    optimization_records: dict


def fit_one_vs_rest(X, labels, classes, optimization):
    labels = np.asarray(labels)
    states = {}
    records = {}
    for label in classes:
        targets = np.where(labels == label, 1.0, -1.0)
        state, class_records = optimize_hyperparameters(
            X,
            targets,
            optimization["initial_starts"],
            lengthscale_bounds=tuple(optimization["lengthscale_bounds"]),
            signal_std_bounds=tuple(optimization["signal_std_bounds"]),
            jitter=float(optimization["jitter"]),
            maximum_laplace_iterations=int(optimization["maximum_laplace_iterations"]),
            laplace_tolerance=float(optimization["laplace_tolerance"]),
            maximum_optimizer_iterations=int(optimization["maximum_optimizer_iterations"]),
            function_tolerance=float(optimization["function_tolerance"]),
            log_prior_mean=np.log(np.asarray(optimization["prior_mean_parameters"], dtype=float)),
            log_prior_std=np.asarray(optimization["prior_log_std"], dtype=float),
        )
        if not state.converged:
            raise RuntimeError(f"Laplace inference did not converge for class {label}")
        states[label] = state
        records[label] = class_records
    return OneVsRestModel(tuple(classes), states, records)


def predict_normalized_probabilities(model, X):
    raw = np.column_stack([predict_proba(model.states[label], X)[0] for label in model.classes])
    return raw / np.maximum(raw.sum(axis=1, keepdims=True), 1e-15)


def sample_joint_class_indices(model, X, sample_count, random_seed, eigenvalue_floor=1e-10):
    rng = np.random.default_rng(int(random_seed))
    latent = np.empty((int(sample_count), len(X), len(model.classes)))
    diagnostics = {}
    for class_index, label in enumerate(model.classes):
        mean, covariance = predict_latent_joint(model.states[label], X)
        eigenvalues, eigenvectors = np.linalg.eigh(covariance)
        clipped = np.maximum(eigenvalues, float(eigenvalue_floor))
        standard = rng.standard_normal((int(sample_count), len(X)))
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            latent[:, :, class_index] = mean + standard @ (eigenvectors * np.sqrt(clipped)).T
        if not np.all(np.isfinite(latent[:, :, class_index])):
            raise FloatingPointError(f"Nonfinite joint latent samples for class {label}")
        diagnostics[label] = {
            "minimum_covariance_eigenvalue": float(eigenvalues.min()),
            "negative_eigenvalue_count": int(np.sum(eigenvalues < 0.0)),
        }
    return np.argmax(latent, axis=2), diagnostics
