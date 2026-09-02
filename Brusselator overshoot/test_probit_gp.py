"""Minimal numerical tests for the local probit GP implementation."""

import numpy as np

from probit_gp import optimize_hyperparameters, predict_latent_joint, predict_proba


def main():
    X = np.array(
        [
            [-1.0, -1.0],
            [-1.0, 0.0],
            [0.0, -1.0],
            [1.0, 1.0],
            [1.0, 0.0],
            [0.0, 1.0],
        ]
    )
    y = np.array([-1.0, -1.0, -1.0, 1.0, 1.0, 1.0])
    state, records = optimize_hyperparameters(
        X,
        y,
        starts=[[0.5, 0.5, 1.0], [1.0, 1.0, 1.0]],
        maximum_optimizer_iterations=80,
    )
    probabilities, means, variances = predict_proba(
        state, np.array([[-0.8, -0.8], [0.8, 0.8], [0.0, 0.0]])
    )
    assert state.converged
    assert len(records) == 2
    assert np.all(np.isfinite(probabilities))
    assert np.all(variances >= 0.0)
    assert probabilities[0] < 0.5 < probabilities[1]
    assert abs(probabilities[2] - 0.5) < 0.15
    assert means[0] < means[1]
    joint_points = np.array([[-0.1, -0.1], [0.1, 0.1], [0.8, 0.8]])
    joint_mean, joint_covariance = predict_latent_joint(state, joint_points)
    assert joint_mean.shape == (3,)
    assert joint_covariance.shape == (3, 3)
    assert np.max(np.abs(joint_covariance - joint_covariance.T)) < 1e-10
    assert np.min(np.linalg.eigvalsh(joint_covariance)) > -1e-8
    assert joint_covariance[0, 1] > joint_covariance[0, 2]

    prior_state, prior_records = optimize_hyperparameters(
        X,
        y,
        starts=[[0.5, 0.5, 1.0], [1.0, 1.0, 1.0]],
        maximum_optimizer_iterations=80,
        log_prior_mean=np.log([1.0, 1.0, 2.0]),
        log_prior_std=np.array([1.0, 1.0, 1.0]),
    )
    assert prior_state.converged and len(prior_records) == 2
    print("probit_gp synthetic test passed")


if __name__ == "__main__":
    main()
