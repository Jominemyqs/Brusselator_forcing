clear; clc; close all;

%% ============================================================
% PARAMETERS
% ============================================================

sigma  = 0.45;
b0     = 10.0;
Tfinal = 2000;
flag   = 0;

% Cases to inspect
lambda_list = [0.74, 0.75, 0.76];

% Save folder
outdir = 'post_threshold_regime_outputs';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

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
% BUILD CASES
% ============================================================

case_names = ["B", "lambda_0.74", "lambda_0.75", "lambda_0.76"];
num_cases  = length(case_names);

ic_cases = cell(num_cases,1);

% Case 1: B itself
ic_cases{1} = [x, uB, vB];

% Cases 2-4: interpolated states
for k = 1:length(lambda_list)
    lam = lambda_list(k);
    u0 = (1-lam)*uA + lam*uB;
    v0 = (1-lam)*vA + lam*vB;
    ic_cases{1+k} = [x, u0, v0];
end

%% ============================================================
% MAIN LOOP
% ============================================================

for i = 1:num_cases
    fprintf('\n====================================\n');
    fprintf('Running case: %s\n', case_names(i));
    fprintf('====================================\n');

    ic0 = ic_cases{i};

    [ic_final, S, V] = solve_brusselator_1d_forced(ic0, par, Tfinal, flag);

    v_final = ic_final(:,3);

    % Sample a few spatial points for time series
    idx1 = round(N/4);
    idx2 = round(N/2);
    idx3 = round(3*N/4);

    ts1 = V(:, idx1);
    ts2 = V(:, idx2);
    ts3 = V(:, idx3);

    % Late-time snapshot indices
    snap_times = [1000, 1500, 2000];
    snap_inds  = zeros(size(snap_times));
    for k = 1:length(snap_times)
        [~, snap_inds(k)] = min(abs(S - snap_times(k)));
    end

    case_dir = fullfile(outdir, char(case_names(i)));
    if ~exist(case_dir, 'dir')
        mkdir(case_dir);
    end

    %% --------------------------------------------------------
    % Space-time plot
    %% --------------------------------------------------------
    fig = figure('Visible','off');
    imagesc(x, S, V);
    axis xy;
    xlabel('x');
    ylabel('t');
    title(sprintf('Space-time plot of v(x,t): %s', case_names(i)));
    colorbar;
    box on;
    exportgraphics(fig, fullfile(case_dir, 'spacetime_v.png'));
    close(fig);

    %% --------------------------------------------------------
    % Late-time snapshots
    %% --------------------------------------------------------
    fig = figure('Visible','off');
    plot(x, vA, 'k-', 'LineWidth', 2.0, 'DisplayName', 'A');
    hold on;
    plot(x, vB, 'k--', 'LineWidth', 2.0, 'DisplayName', 'B');

    for k = 1:length(snap_inds)
        plot(x, V(snap_inds(k),:), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('t = %.0f', S(snap_inds(k))));
    end

    plot(x, v_final, 'LineWidth', 1.8, 'DisplayName', 'final');
    xlabel('x');
    ylabel('v(x,t)');
    title(sprintf('Late-time profiles: %s', case_names(i)));
    legend('Location','eastoutside');
    box on;
    exportgraphics(fig, fullfile(case_dir, 'late_time_profiles.png'));
    close(fig);

    %% --------------------------------------------------------
    % Time series at selected spatial points
    %% --------------------------------------------------------
    fig = figure('Visible','off');
    plot(S, ts1, 'LineWidth', 1.4, 'DisplayName', sprintf('x index %d', idx1));
    hold on;
    plot(S, ts2, 'LineWidth', 1.4, 'DisplayName', sprintf('x index %d', idx2));
    plot(S, ts3, 'LineWidth', 1.4, 'DisplayName', sprintf('x index %d', idx3));
    xlabel('t');
    ylabel('v(t,x_*)');
    title(sprintf('Time series at selected spatial points: %s', case_names(i)));
    legend('Location','best');
    box on;
    exportgraphics(fig, fullfile(case_dir, 'time_series_selected_points.png'));
    close(fig);

    %% --------------------------------------------------------
    % Save numeric data
    %% --------------------------------------------------------
    save(fullfile(case_dir, 'case_data.mat'), ...
        'x', 'S', 'V', 'ic0', 'ic_final', 'v_final', ...
        'snap_times', 'snap_inds', 'idx1', 'idx2', 'idx3', ...
        'ts1', 'ts2', 'ts3', 'uA', 'vA', 'uB', 'vB');
end

disp(' ');
disp(['Post-threshold characterization complete. Outputs saved in: ', outdir]);