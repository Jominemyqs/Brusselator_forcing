clear; clc; close all;

%% ============================================================
% USER SETTINGS
% ============================================================

sigma = 0.45;

% Baseline values to test
b0_list = [9.8, 10.0, 10.2];

% Overshoot protocol
Tup    = 40;
Thold  = 80;
Tdown  = 40;
Tfinal_forced = 300;

% Long baseline relaxation for lambda classification
Tfinal_lambda = 2000;

% Choose a supercritical overshoot relative to b0
delta_b = 1.2;   % you may tune this if needed

% Lambda scan
lambda_list = 0.65:0.01:0.80;

% Threshold for "returns to A"
tol_A = 1e-2;

flag = 0;

%% ============================================================
% INITIAL GRID / IC TEMPLATE
% ============================================================

N  = 400;
Lx = 40;
x  = linspace(0, Lx, N)';

a = 2.5;

%% ============================================================
% STORAGE
% ============================================================

nB0 = length(b0_list);
nL  = length(lambda_list);

lambda_threshold_low  = nan(nB0,1);
lambda_threshold_high = nan(nB0,1);

return_matrix = false(nB0, nL);

%% ============================================================
% MAIN LOOP OVER b0
% ============================================================

for ib = 1:nB0
    b0 = b0_list(ib);
    bmax = b0 + delta_b;

    fprintf('\n=================================================\n');
    fprintf('Processing b0 = %.3f, bmax = %.3f\n', b0, bmax);
    fprintf('=================================================\n');

    %% --------------------------------------------------------
    % Build initial condition near homogeneous equilibrium at b0
    %% --------------------------------------------------------
    u1_ss = a;
    u2_ss = b0 / a;

    rng(1);
    u1 = u1_ss + 0.01 * randn(N,1);
    u2 = u2_ss + 0.01 * randn(N,1);
    ic0 = [x, u1, u2];

    %% --------------------------------------------------------
    % Compute control final state A at constant b0
    %% --------------------------------------------------------
    par_control.sigma = sigma;
    par_control.Bfun  = @(t) b0;

    [icA, ~, ~] = solve_brusselator_1d_forced(ic0, par_control, Tfinal_lambda, flag);
    uA = icA(:,2);
    vA = icA(:,3);

    %% --------------------------------------------------------
    % Compute overshoot final state B
    %% --------------------------------------------------------
    par_forced.sigma = sigma;
    par_forced.Bfun  = @(t) overshoot_B(t, b0, bmax, Tup, Thold, Tdown);

    [icB, ~, ~] = solve_brusselator_1d_forced(ic0, par_forced, Tfinal_forced, flag);
    uB = icB(:,2);
    vB = icB(:,3);

    %% --------------------------------------------------------
    % Scan lambda along slice from A to B
    %% --------------------------------------------------------
    for il = 1:nL
        lam = lambda_list(il);

        u0 = (1-lam)*uA + lam*uB;
        v0 = (1-lam)*vA + lam*vB;
        ic_lam = [x, u0, v0];

        [ic_final, ~, ~] = solve_brusselator_1d_forced(ic_lam, par_control, Tfinal_lambda, flag);
        v_final = ic_final(:,3);

        distA = norm(v_final - vA, 2) / max(norm(vA,2), 1e-12);
        return_matrix(ib, il) = (distA < tol_A);

        fprintf('b0 = %.3f, lambda = %.2f, distA = %.3e, return = %d\n', ...
            b0, lam, distA, return_matrix(ib, il));
    end

    %% --------------------------------------------------------
    % Extract threshold interval
    %% --------------------------------------------------------
    idx_return = find(return_matrix(ib,:)==1, 1, 'last');
    idx_leave  = find(return_matrix(ib,:)==0, 1, 'first');

    if ~isempty(idx_return)
        lambda_threshold_low(ib) = lambda_list(idx_return);
    end
    if ~isempty(idx_leave)
        lambda_threshold_high(ib) = lambda_list(idx_leave);
    end
end

%% ============================================================
% SUMMARY TABLE
% ============================================================

summary_table = table(b0_list(:), lambda_threshold_low, lambda_threshold_high, ...
    'VariableNames', {'b0','lambda_threshold_low','lambda_threshold_high'});

disp(' ');
disp(summary_table);

%% ============================================================
% PLOT
% ============================================================

figure;
plot(b0_list, lambda_threshold_low, '-o', 'LineWidth', 1.8, 'DisplayName', 'last return to A');
hold on;
plot(b0_list, lambda_threshold_high, '-s', 'LineWidth', 1.8, 'DisplayName', 'first non-return');
xlabel('b_0');
ylabel('Approximate threshold interval for \lambda^*');
title('Mapping the slice threshold \lambda^*(b_0)');
legend('Location','best');
grid on;
hold off;

%% ============================================================
% SAVE
% ============================================================

save('map_lambda_threshold_vs_b0_results.mat', ...
    'b0_list', 'lambda_list', 'return_matrix', ...
    'lambda_threshold_low', 'lambda_threshold_high', 'summary_table');