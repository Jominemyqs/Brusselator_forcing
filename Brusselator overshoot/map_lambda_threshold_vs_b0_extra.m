clear; clc; close all;

%% ============================================================
% map_lambda_threshold_vs_b0_extra.m
%
% Targeted extra b0 mapping for Ritchie-inspired Brusselator project.
%
% Goal:
%   Add b0 = 9.9 and b0 = 10.1 to the existing threshold map.
%
% Existing thresholds:
%   b0 = 9.8  -> lambda* in (0.74, 0.75)
%   b0 = 10.0 -> lambda* in (0.72, 0.73)
%   b0 = 10.2 -> lambda* in (0.55, 0.60)
%
% This script computes:
%   b0 = 9.9
%   b0 = 10.1
%
% For each b0:
%   1. Compute baseline state A under frozen b = b0.
%   2. Compute forced post-overshoot state B.
%   3. Interpolate between A and B:
%        IC(lambda) = (1-lambda) A + lambda B.
%   4. Evolve under frozen b = b0.
%   5. Determine whether the trajectory returns to A.
%
% Requires:
%   solve_brusselator_1d_forced.m
%   overshoot_B.m
% ============================================================

%% -------------------- User settings --------------------

sigma = 0.45;

% Only run the missing extra b0 values
b0_list = [9.9, 10.1];

% Use same overshoot size as previous mapping:
% bmax = b0 + 1.2
delta_b = 1.2;

% Overshoot protocol
Tup   = 40;
Thold = 80;
Tdown = 40;

% Time horizons
Tfinal_forced = 300;     % for generating B
Tfinal_lambda = 2000;    % for testing return to A

% Lambda scan range
% This is broad enough to catch thresholds for both b0 = 9.9 and 10.1.
lambda_list = 0.50:0.02:0.80;

% Return threshold
tol_A = 1e-2;

% Solver plotting flag
flag = 0;

% Output folder
outdir = 'lambda_threshold_b0_extra_outputs';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

%% -------------------- Spatial grid and initial condition --------------------

N  = 400;
Lx = 40;
x  = linspace(0, Lx, N)';

a = 2.5;

% We will build the initial condition separately for each b0,
% because v* = b0/a depends on b0.
rng(1);

%% -------------------- Storage --------------------

nB = length(b0_list);
nL = length(lambda_list);

threshold_low  = NaN(nB,1);
threshold_high = NaN(nB,1);

all_distA = NaN(nB,nL);
all_return_flag = false(nB,nL);

%% ============================================================
% Main loop over b0 values
% ============================================================

for ib = 1:nB

    b0 = b0_list(ib);
    bmax = b0 + delta_b;

    fprintf('\n=================================================\n');
    fprintf('Processing b0 = %.3f, bmax = %.3f\n', b0, bmax);
    fprintf('=================================================\n');

    %% -------------------- Initial condition near homogeneous state --------------------

    u1_ss = a;
    u2_ss = b0/a;

    rng(1);  % keep initial perturbation comparable across b0
    u1 = u1_ss + 0.01*randn(N,1);
    u2 = u2_ss + 0.01*randn(N,1);

    ic0 = [x, u1, u2];

    %% -------------------- Compute baseline state A --------------------

    fprintf('\nComputing baseline state A...\n');

    par_control.sigma = sigma;
    par_control.Bfun  = @(t) b0;

    [icA, SA, VA] = solve_brusselator_1d_forced(ic0, par_control, Tfinal_lambda, flag);

    uA = icA(:,2);
    vA = icA(:,3);

    %% -------------------- Compute forced post-overshoot state B --------------------

    fprintf('Computing forced post-overshoot state B...\n');

    par_forced.sigma = sigma;
    par_forced.Bfun  = @(t) overshoot_B(t, b0, bmax, Tup, Thold, Tdown);

    [icB, SB, VB] = solve_brusselator_1d_forced(ic0, par_forced, Tfinal_forced, flag);

    uB = icB(:,2);
    vB = icB(:,3);

    %% -------------------- Save A and B for this b0 --------------------

    b0_tag = sprintf('b0_%04.1f', b0);
    b0_tag = strrep(b0_tag, '.', 'p');

    save(fullfile(outdir, [b0_tag '_AB_states.mat']), ...
        'b0','bmax','sigma','x','icA','icB','uA','vA','uB','vB', ...
        'SA','VA','SB','VB','Tup','Thold','Tdown','Tfinal_forced','Tfinal_lambda');

    %% -------------------- Lambda scan --------------------

    distA = NaN(nL,1);
    return_flag = false(nL,1);

    for il = 1:nL

        lam = lambda_list(il);

        fprintf('\n-------------------------------------\n');
        fprintf('b0 = %.3f, lambda = %.3f\n', b0, lam);
        fprintf('-------------------------------------\n');

        u0_lam = (1-lam)*uA + lam*uB;
        v0_lam = (1-lam)*vA + lam*vB;

        ic_lam = [x, u0_lam, v0_lam];

        [ic_final, ~, ~] = solve_brusselator_1d_forced(ic_lam, par_control, Tfinal_lambda, flag);

        v_final = ic_final(:,3);

        distA(il) = norm(v_final - vA, 2) / max(norm(vA,2), 1e-12);
        return_flag(il) = distA(il) < tol_A;

        fprintf('distA = %.6e, return = %d\n', distA(il), return_flag(il));

        % Save partial progress after each lambda in case the run is interrupted
        partial_table = table(lambda_list(:), distA, return_flag, ...
            'VariableNames', {'lambda','dist_to_A','returns_to_A'});

        save(fullfile(outdir, [b0_tag '_partial_results.mat']), ...
            'b0','bmax','sigma','lambda_list','distA','return_flag', ...
            'tol_A','Tfinal_lambda','Tfinal_forced','Tup','Thold','Tdown');

        writetable(partial_table, fullfile(outdir, [b0_tag '_partial_results.csv']));
    end

    %% -------------------- Determine threshold interval --------------------

    idx_return = find(return_flag == 1, 1, 'last');
    idx_leave  = find(return_flag == 0, 1, 'first');

    if ~isempty(idx_return)
        threshold_low(ib) = lambda_list(idx_return);
    end

    if ~isempty(idx_leave)
        threshold_high(ib) = lambda_list(idx_leave);
    end

    all_distA(ib,:) = distA(:)';
    all_return_flag(ib,:) = return_flag(:)';

    %% -------------------- Save b0-specific final result --------------------

    result_table = table(lambda_list(:), distA, return_flag, ...
        'VariableNames', {'lambda','dist_to_A','returns_to_A'});

    threshold_table_single = table(b0, bmax, threshold_low(ib), threshold_high(ib), ...
        'VariableNames', {'b0','bmax','lambda_threshold_low','lambda_threshold_high'});

    disp(' ');
    disp('Lambda scan result:');
    disp(result_table);

    disp(' ');
    disp('Threshold interval:');
    disp(threshold_table_single);

    save(fullfile(outdir, [b0_tag '_final_results.mat']), ...
        'b0','bmax','sigma','lambda_list','distA','return_flag', ...
        'threshold_table_single','tol_A','Tfinal_lambda','Tfinal_forced', ...
        'Tup','Thold','Tdown');

    writetable(result_table, fullfile(outdir, [b0_tag '_lambda_scan.csv']));
    writetable(threshold_table_single, fullfile(outdir, [b0_tag '_threshold.csv']));

    %% -------------------- Plot b0-specific scan --------------------

    fig = figure('Color','w','Position',[100 100 800 500]);
    plot(lambda_list, distA, '-o', 'LineWidth', 2, 'MarkerSize', 7);
    hold on;
    yline(tol_A, '--', 'return threshold', 'LineWidth', 1.2);

    if ~isnan(threshold_low(ib))
        xline(threshold_low(ib), '--', 'LineWidth', 1.2, 'HandleVisibility','off');
    end
    if ~isnan(threshold_high(ib))
        xline(threshold_high(ib), '--', 'LineWidth', 1.2, 'HandleVisibility','off');
    end

    xlabel('\lambda','Interpreter','latex');
    ylabel('Relative distance to $A$','Interpreter','latex');
    title(sprintf('Basin-slice scan for $b_0=%.1f$', b0), 'Interpreter','latex');
    grid on;
    box on;
    set(gca,'FontSize',13);

    exportgraphics(fig, fullfile(outdir, [b0_tag '_dist_to_A_vs_lambda.png']), 'Resolution', 300);
    close(fig);

end

%% ============================================================
% Combined summary table
% ============================================================

summary_table = table(b0_list(:), b0_list(:)+delta_b, threshold_low, threshold_high, ...
    'VariableNames', {'b0','bmax','lambda_threshold_low','lambda_threshold_high'});

disp(' ');
disp('=================================================');
disp('Extra b0 threshold mapping complete.');
disp('=================================================');
disp(summary_table);

writetable(summary_table, fullfile(outdir, 'extra_b0_threshold_summary.csv'));

save(fullfile(outdir, 'extra_b0_threshold_summary.mat'), ...
    'b0_list','delta_b','lambda_list','threshold_low','threshold_high', ...
    'all_distA','all_return_flag','summary_table','tol_A', ...
    'Tfinal_lambda','Tfinal_forced','Tup','Thold','Tdown');

%% ============================================================
% Combined plot: extra b0 values only
% ============================================================

fig = figure('Color','w','Position',[100 100 800 500]);
hold on;

for ib = 1:nB
    plot(lambda_list, all_distA(ib,:), '-o', 'LineWidth', 2, 'MarkerSize', 6, ...
        'DisplayName', sprintf('$b_0=%.1f$', b0_list(ib)));
end

yline(tol_A, '--', 'return threshold', 'LineWidth', 1.2, 'HandleVisibility','off');

xlabel('\lambda','Interpreter','latex');
ylabel('Relative distance to $A$','Interpreter','latex');
title('Extra basin-slice scans for intermediate $b_0$', 'Interpreter','latex');
legend('Interpreter','latex','Location','best');
grid on;
box on;
set(gca,'FontSize',13);

exportgraphics(fig, fullfile(outdir, 'extra_b0_dist_to_A_comparison.png'), 'Resolution', 300);
close(fig);

disp(['All outputs saved in folder: ', outdir]);