clear; clc; close all;

%% ============================================================
% PARAMETERS
% ============================================================

sigma  = 0.45;
b0     = 10.0;
Tfinal = 300;
flag   = 0;

lambda_list = 0:0.05:1;

recovery_thresh   = 1e-2;
transition_thresh = 1e-1;

%% ============================================================
% FILE PATHS
% ============================================================

control_file    = fullfile('Thold_bmax_sweep_fast_outputs', 'control_results.mat');
transition_file = fullfile('Thold_bmax_sweep_fast_outputs', 'run_Thold_80_bmax_11.14.mat');

%% ============================================================
% LOAD STATE A: CONTROL FINAL STATE
% ============================================================

C = load(control_file);

uA = C.ic_control(:,2);
vA = C.ic_control(:,3);
x  = C.ic_control(:,1);

%% ============================================================
% LOAD STATE B: TRANSITION FINAL STATE
% ============================================================

T = load(transition_file);

uB = T.run_data.ic_forced(:,2);
vB = T.run_data.ic_forced(:,3);

%% ============================================================
% BASELINE SOLVER
% ============================================================

par.sigma = sigma;
par.Bfun  = @(t) b0;

%% ============================================================
% STORAGE
% ============================================================

nL = length(lambda_list);
rel_diff_to_A = zeros(nL,1);
rel_diff_to_B = zeros(nL,1);
class_label   = strings(nL,1);

%% ============================================================
% LOOP OVER LAMBDA
% ============================================================

for i = 1:nL
    lambda = lambda_list(i);

    fprintf('\n====================================\n');
    fprintf('Running lambda = %.3f\n', lambda);

    % Interpolated IC
    u0 = (1-lambda)*uA + lambda*uB;
    v0 = (1-lambda)*vA + lambda*vB;
    ic0 = [x, u0, v0];

    % Constant-b0 simulation
    [ic_final, ~, ~] = solve_brusselator_1d_forced(ic0, par, Tfinal, flag);

    u_final = ic_final(:,2);
    v_final = ic_final(:,3);

    % Compare to both candidate states
    rel_diff_to_A(i) = norm(v_final - vA,2) / max(norm(vA,2),1e-12);
    rel_diff_to_B(i) = norm(v_final - vB,2) / max(norm(vB,2),1e-12);

    % Classification by whichever attractor is closer
    if rel_diff_to_A(i) < recovery_thresh
        class_label(i) = "ReturnsToA";
    elseif rel_diff_to_B(i) < recovery_thresh
        class_label(i) = "ReturnsToB";
    else
        if rel_diff_to_A(i) < rel_diff_to_B(i)
            class_label(i) = "CloserToA";
        else
            class_label(i) = "CloserToB";
        end
    end

    fprintf('diff to A = %.3e, diff to B = %.3e → %s\n', ...
        rel_diff_to_A(i), rel_diff_to_B(i), class_label(i));
end

%% ============================================================
% SUMMARY TABLE
% ============================================================

summary_table = table(lambda_list(:), rel_diff_to_A, rel_diff_to_B, class_label, ...
    'VariableNames', {'lambda','rel_diff_to_A','rel_diff_to_B','class_label'});

disp(' ');
disp(summary_table);

%% ============================================================
% PLOTS
% ============================================================

figure;
plot(lambda_list, rel_diff_to_A, '-o', 'LineWidth', 1.8, 'DisplayName', 'distance to A');
hold on;
plot(lambda_list, rel_diff_to_B, '-s', 'LineWidth', 1.8, 'DisplayName', 'distance to B');
yline(recovery_thresh, '--', 'threshold');
xlabel('\lambda');
ylabel('Relative difference');
title('Basin-slice experiment');
legend('Location','best');
grid on;