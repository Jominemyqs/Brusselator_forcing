clear; clc; close all;

%% ============================================================
% PARAMETERS
% ============================================================

sigma = 0.45;
b0    = 10.0;

% Long time horizon
Tfinal = 2000;
flag   = 0;

% Lambda test cases near threshold
lambda_list = [0.74, 0.75, 0.76];

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

% A = control final state
x  = C.ic_control(:,1);
uA = C.ic_control(:,2);
vA = C.ic_control(:,3);

% B = transitioned overshoot final state
uB = T.run_data.ic_forced(:,2);
vB = T.run_data.ic_forced(:,3);

N = length(x);

%% ============================================================
% BASELINE SOLVER PARAMETER
% ============================================================

par.sigma = sigma;
par.Bfun  = @(t) b0;

%% ============================================================
% CASE DEFINITIONS
% ============================================================

case_names = ["A", "B", "lambda_0.74", "lambda_0.75", "lambda_0.76"];
num_cases  = length(case_names);

ic_cases = cell(num_cases,1);

% Case 1: start from A
ic_cases{1} = [x, uA, vA];

% Case 2: start from B
ic_cases{2} = [x, uB, vB];

% Cases 3-5: interpolated states
for k = 1:length(lambda_list)
    lam = lambda_list(k);
    u0 = (1-lam)*uA + lam*uB;
    v0 = (1-lam)*vA + lam*vB;
    ic_cases{2+k} = [x, u0, v0];
end

%% ============================================================
% STORAGE
% ============================================================

t_store = cell(num_cases,1);
distA_store = cell(num_cases,1);
distB_store = cell(num_cases,1);
vfinal_store = cell(num_cases,1);
ufinal_store = cell(num_cases,1);

%% ============================================================
% MAIN LOOP
% ============================================================

for i = 1:num_cases
    fprintf('\n====================================\n');
    fprintf('Running case: %s\n', case_names(i));
    fprintf('====================================\n');

    ic0 = ic_cases{i};

    % Run solver
    [ic_final, S, V] = solve_brusselator_1d_forced(ic0, par, Tfinal, flag);

    % Recover u from initial/final only? We need full u trajectory too.
    % Since solve_brusselator_1d_forced only returns V and final ic,
    % we use V(t,:) for distance tracking in v-component only.
    %
    % Compare each time slice of V to vA and vB
    nt = length(S);
    distA = zeros(nt,1);
    distB = zeros(nt,1);

    for j = 1:nt
        vj = V(j,:)';
        distA(j) = norm(vj - vA, 2) / max(norm(vA,2),1e-12);
        distB(j) = norm(vj - vB, 2) / max(norm(vB,2),1e-12);
    end

    t_store{i}    = S;
    distA_store{i} = distA;
    distB_store{i} = distB;

    ufinal_store{i} = ic_final(:,2);
    vfinal_store{i} = ic_final(:,3);

    fprintf('Final distance to A: %.6e\n', distA(end));
    fprintf('Final distance to B: %.6e\n', distB(end));
end

%% ============================================================
% SUMMARY TABLE
% ============================================================

final_distA = zeros(num_cases,1);
final_distB = zeros(num_cases,1);

for i = 1:num_cases
    final_distA(i) = distA_store{i}(end);
    final_distB(i) = distB_store{i}(end);
end

summary_table = table(case_names(:), final_distA, final_distB, ...
    'VariableNames', {'case_name','final_dist_to_A','final_dist_to_B'});

disp(' ');
disp(summary_table);

%% ============================================================
% PLOT 1: DISTANCE TO A OVER TIME
% ============================================================

figure;
hold on;
for i = 1:num_cases
    plot(t_store{i}, distA_store{i}, 'LineWidth', 1.8, ...
        'DisplayName', char(case_names(i)));
end
xlabel('t');
ylabel('Relative distance to A');
title('Distance to A over time at baseline b_0');
legend('Location','best');
grid on;
hold off;

%% ============================================================
% PLOT 2: DISTANCE TO B OVER TIME
% ============================================================

figure;
hold on;
for i = 1:num_cases
    plot(t_store{i}, distB_store{i}, 'LineWidth', 1.8, ...
        'DisplayName', char(case_names(i)));
end
xlabel('t');
ylabel('Relative distance to B');
title('Distance to B over time at baseline b_0');
legend('Location','best');
grid on;
hold off;

%% ============================================================
% PLOT 3: FINAL PROFILES
% ============================================================

figure;
plot(x, vA, 'k-', 'LineWidth', 2.0, 'DisplayName', 'A');
hold on;
plot(x, vB, 'k--', 'LineWidth', 2.0, 'DisplayName', 'B');

for i = 1:num_cases
    plot(x, vfinal_store{i}, 'LineWidth', 1.4, ...
        'DisplayName', sprintf('%s final', case_names(i)));
end

xlabel('x');
ylabel('v(x,T)');
title(sprintf('Final profiles at T = %d', Tfinal));
legend('Location','eastoutside');
grid on;
hold off;

%% ============================================================
% SAVE OUTPUT
% ============================================================

save('long_run_baseline_AB_test_results.mat', ...
    'case_names', 't_store', 'distA_store', 'distB_store', ...
    'final_distA', 'final_distB', 'summary_table', ...
    'x', 'uA', 'vA', 'uB', 'vB', 'ufinal_store', 'vfinal_store', ...
    'Tfinal', 'lambda_list');