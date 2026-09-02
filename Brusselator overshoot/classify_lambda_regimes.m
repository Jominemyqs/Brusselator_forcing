clear; clc; close all;

%% ============================================================
% PARAMETERS
% ============================================================

sigma = 0.45;
b0    = 10.0;

% Long baseline evolution for classification
Tfinal_long = 3000;

% Perturbed restart length
Tfinal_pert = 1500;

flag = 0;

% Lambda grid to classify
lambda_list = 0.70:0.01:0.80;

% Perturbation size
pert_amp = 1e-2;

% Distance thresholds
tol_AorB      = 5e-2;   % "close to A/B"
tol_same      = 5e-2;   % restart returns to same candidate
tol_separate  = 1e-1;   % clearly not close

%% ============================================================
% FILE PATHS
% ============================================================

control_file    = fullfile('Thold_bmax_sweep_fast_outputs', 'control_results.mat');
transition_file = fullfile('Thold_bmax_sweep_fast_outputs', 'run_Thold_80_bmax_11.14.mat');

outdir = 'lambda_regime_classification_outputs';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

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
% STORAGE
% ============================================================

nL = length(lambda_list);

distA_long      = zeros(nL,1);
distB_long      = zeros(nL,1);
distSelf_restart = zeros(nL,1);
distA_restart   = zeros(nL,1);
distB_restart   = zeros(nL,1);

regime_label    = strings(nL,1);

candidate_v_store = cell(nL,1);
restart_v_store   = cell(nL,1);

%% ============================================================
% MAIN LOOP
% ============================================================

for i = 1:nL
    lambda = lambda_list(i);

    fprintf('\n====================================\n');
    fprintf('Running lambda = %.2f\n', lambda);
    fprintf('====================================\n');

    %% --------------------------------------------------------
    % Build interpolated initial condition
    %% --------------------------------------------------------
    u0 = (1-lambda)*uA + lambda*uB;
    v0 = (1-lambda)*vA + lambda*vB;
    ic0 = [x, u0, v0];

    %% --------------------------------------------------------
    % Run long baseline trajectory
    %% --------------------------------------------------------
    [ic_long, ~, ~] = solve_brusselator_1d_forced(ic0, par, Tfinal_long, flag);

    u_long = ic_long(:,2);
    v_long = ic_long(:,3);

    candidate_v_store{i} = v_long;

    distA_long(i) = norm(v_long - vA, 2) / max(norm(vA,2), 1e-12);
    distB_long(i) = norm(v_long - vB, 2) / max(norm(vB,2), 1e-12);

    %% --------------------------------------------------------
    % Perturb candidate final state and restart
    %% --------------------------------------------------------
    rng(1000 + i); % reproducible but different per lambda
    u_pert = u_long + pert_amp * randn(size(u_long));
    v_pert = v_long + pert_amp * randn(size(v_long));
    ic_pert = [x, u_pert, v_pert];

    [ic_restart, ~, ~] = solve_brusselator_1d_forced(ic_pert, par, Tfinal_pert, flag);

    u_restart = ic_restart(:,2);
    v_restart = ic_restart(:,3);

    restart_v_store{i} = v_restart;

    distSelf_restart(i) = norm(v_restart - v_long, 2) / max(norm(v_long,2), 1e-12);
    distA_restart(i)    = norm(v_restart - vA,    2) / max(norm(vA,2),    1e-12);
    distB_restart(i)    = norm(v_restart - vB,    2) / max(norm(vB,2),    1e-12);

    %% --------------------------------------------------------
    % Classification logic
    %% --------------------------------------------------------
    if distA_long(i) < tol_AorB && distA_restart(i) < tol_AorB
        regime_label(i) = "A";
    elseif distB_long(i) < tol_AorB && distB_restart(i) < tol_AorB
        regime_label(i) = "B";
    elseif distSelf_restart(i) < tol_same && ...
           distA_restart(i) > tol_separate && distB_restart(i) > tol_separate
        regime_label(i) = "RobustOther";
    else
        regime_label(i) = "Transitional";
    end

    fprintf('distA_long      = %.4e\n', distA_long(i));
    fprintf('distB_long      = %.4e\n', distB_long(i));
    fprintf('distSelf_restart= %.4e\n', distSelf_restart(i));
    fprintf('distA_restart   = %.4e\n', distA_restart(i));
    fprintf('distB_restart   = %.4e\n', distB_restart(i));
    fprintf('Regime          = %s\n', regime_label(i));
end

%% ============================================================
% SUMMARY TABLE
% ============================================================

summary_table = table(lambda_list(:), distA_long, distB_long, ...
    distSelf_restart, distA_restart, distB_restart, regime_label, ...
    'VariableNames', {'lambda','distA_long','distB_long', ...
    'distSelf_restart','distA_restart','distB_restart','regime_label'});

disp(' ');
disp(summary_table);

writetable(summary_table, fullfile(outdir, 'lambda_regime_summary.csv'));
save(fullfile(outdir, 'lambda_regime_summary.mat'), ...
    'lambda_list', 'distA_long', 'distB_long', ...
    'distSelf_restart', 'distA_restart', 'distB_restart', ...
    'regime_label', 'candidate_v_store', 'restart_v_store', ...
    'x', 'vA', 'vB');

%% ============================================================
% PLOT 1: DISTANCES AFTER LONG RUN
% ============================================================

fig = figure;
plot(lambda_list, distA_long, '-o', 'LineWidth', 1.8, 'DisplayName', 'distance to A');
hold on;
plot(lambda_list, distB_long, '-s', 'LineWidth', 1.8, 'DisplayName', 'distance to B');
yline(tol_AorB, '--', 'close-threshold', 'LineWidth', 1.2);
xlabel('\lambda');
ylabel('Distance after long run');
title('Long-run distances vs \lambda');
legend('Location','best');
grid on;
hold off;
exportgraphics(fig, fullfile(outdir, 'long_run_distances_vs_lambda.png'));
close(fig);

%% ============================================================
% PLOT 2: DISTANCE OF RESTART TO CANDIDATE
% ============================================================

fig = figure;
plot(lambda_list, distSelf_restart, '-o', 'LineWidth', 1.8);
yline(tol_same, '--', 'same-state threshold', 'LineWidth', 1.2);
xlabel('\lambda');
ylabel('Distance after perturbed restart to candidate state');
title('Restart stability vs \lambda');
grid on;
exportgraphics(fig, fullfile(outdir, 'restart_stability_vs_lambda.png'));
close(fig);

%% ============================================================
% PLOT 3: REGIME LABELS AS NUMBERS
% ============================================================

regime_num = zeros(nL,1);
for i = 1:nL
    switch regime_label(i)
        case "A"
            regime_num(i) = 1;
        case "B"
            regime_num(i) = 2;
        case "RobustOther"
            regime_num(i) = 3;
        case "Transitional"
            regime_num(i) = 4;
    end
end

fig = figure;
stairs(lambda_list, regime_num, 'LineWidth', 2);
xlabel('\lambda');
ylabel('Regime code');
title('Regime classification vs \lambda');
yticks([1 2 3 4]);
yticklabels({'A','B','RobustOther','Transitional'});
grid on;
exportgraphics(fig, fullfile(outdir, 'regime_classification_vs_lambda.png'));
close(fig);

disp(' ');
disp(['Done. Outputs saved in: ', outdir]);