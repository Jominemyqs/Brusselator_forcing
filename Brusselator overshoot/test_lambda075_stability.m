clear; clc; close all;

%% ============================================================
% PARAMETERS
% ============================================================

sigma = 0.45;
b0    = 10.0;

% First run long enough to get a candidate state
Tfinal_long = 5000;

% Then rerun from a perturbed version of that state
Tfinal_pert = 3000;

flag = 0;
lambda = 0.75;

% Perturbation size
pert_amp = 1e-2;

%% ============================================================
% FILE PATHS
% ============================================================

control_file    = fullfile('Thold_bmax_sweep_fast_outputs', 'control_results.mat');
transition_file = fullfile('Thold_bmax_sweep_fast_outputs', 'run_Thold_80_bmax_11.14.mat');

%% ============================================================
% LOAD A AND B
% ============================================================

C = load(control_file);
T = load(transition_file);

x  = C.ic_control(:,1);
uA = C.ic_control(:,2);
vA = C.ic_control(:,3);

uB = T.run_data.ic_forced(:,2);
vB = T.run_data.ic_forced(:,3);

N = length(x);

%% ============================================================
% BASELINE PARAMETER
% ============================================================

par.sigma = sigma;
par.Bfun  = @(t) b0;

%% ============================================================
% BUILD lambda = 0.75 INITIAL CONDITION
% ============================================================

u0 = (1-lambda)*uA + lambda*uB;
v0 = (1-lambda)*vA + lambda*vB;
ic0 = [x, u0, v0];

%% ============================================================
% RUN 1: LONG BASELINE EVOLUTION
% ============================================================

fprintf('\n====================================\n');
fprintf('Running long lambda = %.2f trajectory\n', lambda);
fprintf('====================================\n');

[ic_star, S1, V1] = solve_brusselator_1d_forced(ic0, par, Tfinal_long, flag);

u_star = ic_star(:,2);
v_star = ic_star(:,3);

%% ============================================================
% BUILD PERTURBED STATE
% ============================================================

rng(1);
u_pert = u_star + pert_amp * randn(size(u_star));
v_pert = v_star + pert_amp * randn(size(v_star));

ic_pert = [x, u_pert, v_pert];

%% ============================================================
% RUN 2: EVOLVE FROM PERTURBED STATE
% ============================================================

fprintf('\n====================================\n');
fprintf('Running perturbed restart from candidate state\n');
fprintf('====================================\n');

[ic_pert_final, S2, V2] = solve_brusselator_1d_forced(ic_pert, par, Tfinal_pert, flag);

u_pert_final = ic_pert_final(:,2);
v_pert_final = ic_pert_final(:,3);

%% ============================================================
% METRICS
% ============================================================

% Compare perturbed-final to candidate state
rel_diff_star = norm(v_pert_final - v_star, 2) / max(norm(v_star,2),1e-12);

% Compare to A and B too
rel_diff_A = norm(v_pert_final - vA, 2) / max(norm(vA,2),1e-12);
rel_diff_B = norm(v_pert_final - vB, 2) / max(norm(vB,2),1e-12);

fprintf('\n=== Final comparison from perturbed restart ===\n');
fprintf('Distance to candidate state: %.6e\n', rel_diff_star);
fprintf('Distance to A: %.6e\n', rel_diff_A);
fprintf('Distance to B: %.6e\n', rel_diff_B);

%% ============================================================
% PLOT 1: SPACE-TIME OF LONG lambda=0.75 RUN
% ============================================================

figure;
imagesc(x, S1, V1);
axis xy;
xlabel('x');
ylabel('t');
title(sprintf('Space-time plot of v(x,t): lambda = %.2f, long run', lambda));
colorbar;

%% ============================================================
% PLOT 2: SPACE-TIME OF PERTURBED RESTART
% ============================================================

figure;
imagesc(x, S2, V2);
axis xy;
xlabel('x');
ylabel('t');
title('Space-time plot after perturbed restart');
colorbar;

%% ============================================================
% PLOT 3: PROFILE COMPARISON
% ============================================================

figure;
plot(x, vA, 'k-', 'LineWidth', 2.0, 'DisplayName', 'A');
hold on;
plot(x, vB, 'k--', 'LineWidth', 2.0, 'DisplayName', 'B');
plot(x, v_star, 'LineWidth', 2.0, 'DisplayName', 'candidate state');
plot(x, v_pert_final, 'LineWidth', 2.0, 'DisplayName', 'perturbed-final');
xlabel('x');
ylabel('v(x)');
title('Profile comparison');
legend('Location','eastoutside');
grid on;
hold off;

%% ============================================================
% PLOT 4: TIME SERIES AT SELECTED POINTS FOR PERTURBED RUN
% ============================================================

idx1 = round(N/4);
idx2 = round(N/2);
idx3 = round(3*N/4);

figure;
plot(S2, V2(:,idx1), 'LineWidth', 1.4, 'DisplayName', sprintf('x index %d', idx1));
hold on;
plot(S2, V2(:,idx2), 'LineWidth', 1.4, 'DisplayName', sprintf('x index %d', idx2));
plot(S2, V2(:,idx3), 'LineWidth', 1.4, 'DisplayName', sprintf('x index %d', idx3));
xlabel('t');
ylabel('v(t,x_*)');
title('Time series after perturbed restart');
legend('Location','best');
grid on;
hold off;

%% ============================================================
% SAVE OUTPUT
% ============================================================

save('test_lambda075_stability_results.mat', ...
    'x', 'uA', 'vA', 'uB', 'vB', ...
    'lambda', 'pert_amp', ...
    'S1', 'V1', 'u_star', 'v_star', ...
    'S2', 'V2', 'u_pert_final', 'v_pert_final', ...
    'rel_diff_star', 'rel_diff_A', 'rel_diff_B', ...
    'Tfinal_long', 'Tfinal_pert');

disp(' ');
disp('Done. Saved results to test_lambda075_stability_results.mat');