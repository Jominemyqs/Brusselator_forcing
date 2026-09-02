clear; clc; close all;

%% ============================================================
% 1. MASTER OUTPUT FOLDER
% ============================================================
master_outdir = 'bmax_threshold_refine_outputs';
if ~exist(master_outdir, 'dir')
    mkdir(master_outdir);
end

%% ============================================================
% 2. USER SETTINGS
% ============================================================

% Spatial grid
N  = 400;
Lx = 40;
x  = linspace(0, Lx, N)';

% Model parameter
sigma = 0.45;

% Baseline and forcing timing
b0     = 10.0;
Tup    = 40;
Thold  = 80;
Tdown  = 40;
Tfinal = 300;

% Refined sweep values near threshold
% bmax_list = [11.00, 11.10, 11.20, 11.30, 11.40, 11.45, 11.50];
bmax_list = [11.10, 11.12, 11.14, 11.16, 11.18, 11.20];

% Classification thresholds based on rel_diff_u2
recovery_thresh   = 1e-2;
transition_thresh = 1e-1;

% Solver plotting off
flag = 0;

%% ============================================================
% 3. INITIAL CONDITION
% ============================================================
a     = 2.5;
u1_ss = a;
u2_ss = b0 / a;

rng(1);
u1 = u1_ss + 0.01 * randn(N,1);
u2 = u2_ss + 0.01 * randn(N,1);
ic0 = [x, u1, u2];

%% ============================================================
% 4. CONTROL RUN (ONCE)
% ============================================================
disp('Running control simulation...');

par_control.sigma = sigma;
par_control.Bfun  = @(t) b0;

[ic_control, S_control, V_control] = solve_brusselator_1d_forced(ic0, par_control, Tfinal, flag);

u1_control_final = ic_control(:,2);
u2_control_final = ic_control(:,3);

control_outdir = fullfile(master_outdir, 'control');
if ~exist(control_outdir, 'dir')
    mkdir(control_outdir);
end

fig = figure;
imagesc(1:N, S_control, V_control);
axis xy;
xlabel('Spatial index');
ylabel('Time');
title('Control: space-time plot of u_2(x,t)');
colorbar;
box on;
saveas(fig, fullfile(control_outdir, 'control_spacetime_u2.png'));
close(fig);

fig = figure;
plot(x, u2_control_final, 'LineWidth', 1.8);
xlabel('x');
ylabel('u_2(x,T)');
title(sprintf('Control final profile at T = %.1f', Tfinal));
box on;
saveas(fig, fullfile(control_outdir, 'control_final_profile_u2.png'));
close(fig);

save(fullfile(control_outdir, 'control_results.mat'), ...
    'x', 'ic0', 'ic_control', 'S_control', 'V_control', ...
    'u1_control_final', 'u2_control_final', ...
    'sigma', 'b0', 'Tup', 'Thold', 'Tdown', 'Tfinal');

%% ============================================================
% 5. PREALLOCATE SUMMARY STORAGE
% ============================================================
nb = length(bmax_list);

rel_diff_u1 = zeros(nb,1);
rel_diff_u2 = zeros(nb,1);
max_diff_u1 = zeros(nb,1);
max_diff_u2 = zeros(nb,1);
class_label = strings(nb,1);

%% ============================================================
% 6. LOOP OVER BMAX
% ============================================================
for i = 1:nb
    bmax = bmax_list(i);

    fprintf('\n====================================\n');
    fprintf('Running refined threshold simulation for bmax = %.3f\n', bmax);
    fprintf('====================================\n');

    run_name = sprintf('bmax_%0.2f', bmax);
    run_outdir = fullfile(master_outdir, run_name);
    if ~exist(run_outdir, 'dir')
        mkdir(run_outdir);
    end

    par_forced.sigma = sigma;
    par_forced.Bfun  = @(t) overshoot_B(t, b0, bmax, Tup, Thold, Tdown);

    [ic_forced, S_forced, V_forced] = solve_brusselator_1d_forced(ic0, par_forced, Tfinal, flag);

    u1_forced_final = ic_forced(:,2);
    u2_forced_final = ic_forced(:,3);

    rel_diff_u1(i) = norm(u1_forced_final - u1_control_final, 2) / max(norm(u1_control_final, 2), 1e-12);
    rel_diff_u2(i) = norm(u2_forced_final - u2_control_final, 2) / max(norm(u2_control_final, 2), 1e-12);
    max_diff_u1(i) = max(abs(u1_forced_final - u1_control_final));
    max_diff_u2(i) = max(abs(u2_forced_final - u2_control_final));

    if rel_diff_u2(i) < recovery_thresh
        class_label(i) = "Recovery";
    elseif rel_diff_u2(i) > transition_thresh
        class_label(i) = "LastingTransition";
    else
        class_label(i) = "Intermediate";
    end

    fprintf('Relative L2 difference in u1 final state: %.6e\n', rel_diff_u1(i));
    fprintf('Relative L2 difference in u2 final state: %.6e\n', rel_diff_u2(i));
    fprintf('Max absolute difference in u1 final state: %.6e\n', max_diff_u1(i));
    fprintf('Max absolute difference in u2 final state: %.6e\n', max_diff_u2(i));
    fprintf('Classified as: %s\n', class_label(i));

    tt = linspace(0, Tfinal, 1500);
    bb = arrayfun(par_forced.Bfun, tt);

    snapshot_times = [0, 20, 40, 80, 120, 160, 220, 300];
    snapshot_inds_control = zeros(size(snapshot_times));
    snapshot_inds_forced  = zeros(size(snapshot_times));

    for k = 1:length(snapshot_times)
        [~, snapshot_inds_control(k)] = min(abs(S_control - snapshot_times(k)));
        [~, snapshot_inds_forced(k)]  = min(abs(S_forced  - snapshot_times(k)));
    end

    % Forcing plot
    fig = figure;
    plot(tt, b0*ones(size(tt)), 'LineWidth', 1.8, 'DisplayName', 'Control');
    hold on;
    plot(tt, bb, 'LineWidth', 1.8, 'DisplayName', 'Forced');
    xlabel('t');
    ylabel('b(t)');
    title(sprintf('Parameter profile comparison (b_{max}=%.2f)', bmax));
    legend('Location','best');
    box on;
    hold off;
    saveas(fig, fullfile(run_outdir, 'forcing_comparison.png'));
    close(fig);

    % Forced spacetime
    fig = figure;
    imagesc(1:N, S_forced, V_forced);
    axis xy;
    xlabel('Spatial index');
    ylabel('Time');
    title(sprintf('Forced: space-time plot of u_2(x,t), b_{max}=%.2f', bmax));
    colorbar;
    box on;
    saveas(fig, fullfile(run_outdir, 'forced_spacetime_u2.png'));
    close(fig);

    % Final profile comparison
    fig = figure;
    plot(x, u2_control_final, 'LineWidth', 1.8, 'DisplayName', 'Control final');
    hold on;
    plot(x, u2_forced_final, '--', 'LineWidth', 1.8, 'DisplayName', 'Forced final');
    xlabel('x');
    ylabel('u_2(x,T)');
    title(sprintf('Final profile comparison, b_{max}=%.2f', bmax));
    legend('Location','best');
    box on;
    hold off;
    saveas(fig, fullfile(run_outdir, 'final_profile_comparison_u2.png'));
    close(fig);

    % Final difference profile
    fig = figure;
    plot(x, u2_forced_final - u2_control_final, 'LineWidth', 1.8);
    xlabel('x');
    ylabel('\Delta u_2');
    title(sprintf('Final difference profile, b_{max}=%.2f', bmax));
    box on;
    saveas(fig, fullfile(run_outdir, 'final_difference_u2.png'));
    close(fig);

    % Snapshot comparison grid
    fig = figure;
    nrows = 2;
    ncols = ceil(length(snapshot_times)/2);

    for k = 1:length(snapshot_times)
        subplot(nrows, ncols, k);
        plot(x, V_control(snapshot_inds_control(k),:), 'LineWidth', 1.2, 'DisplayName', 'Control');
        hold on;
        plot(x, V_forced(snapshot_inds_forced(k),:), '--', 'LineWidth', 1.2, 'DisplayName', 'Forced');
        title(sprintf('t = %.1f', snapshot_times(k)));
        xlabel('x');
        ylabel('u_2');
        box on;
        hold off;
    end
    legend('Location','best');
    sgtitle(sprintf('Snapshot comparison, b_{max}=%.2f', bmax));
    saveas(fig, fullfile(run_outdir, 'snapshot_comparison_grid.png'));
    close(fig);

    % Save run data
    results = struct();
    results.x = x;
    results.ic0 = ic0;

    results.control.S = S_control;
    results.control.V = V_control;
    results.control.ic_final = ic_control;
    results.control.u1_final = u1_control_final;
    results.control.u2_final = u2_control_final;

    results.forced.S = S_forced;
    results.forced.V = V_forced;
    results.forced.ic_final = ic_forced;
    results.forced.u1_final = u1_forced_final;
    results.forced.u2_final = u2_forced_final;

    results.params.sigma  = sigma;
    results.params.b0     = b0;
    results.params.bmax   = bmax;
    results.params.Tup    = Tup;
    results.params.Thold  = Thold;
    results.params.Tdown  = Tdown;
    results.params.Tfinal = Tfinal;

    results.metrics.rel_diff_u1 = rel_diff_u1(i);
    results.metrics.rel_diff_u2 = rel_diff_u2(i);
    results.metrics.max_diff_u1 = max_diff_u1(i);
    results.metrics.max_diff_u2 = max_diff_u2(i);
    results.class_label = class_label(i);

    results.snapshot_times = snapshot_times;
    results.snapshot_inds_control = snapshot_inds_control;
    results.snapshot_inds_forced  = snapshot_inds_forced;

    save(fullfile(run_outdir, 'results.mat'), 'results');

    % Save text summary
    fid = fopen(fullfile(run_outdir, 'summary.txt'), 'w');
    fprintf(fid, 'Refined forced vs control Brusselator comparison\n\n');
    fprintf(fid, 'sigma  = %.6f\n', sigma);
    fprintf(fid, 'b0     = %.6f\n', b0);
    fprintf(fid, 'bmax   = %.6f\n', bmax);
    fprintf(fid, 'Tup    = %.6f\n', Tup);
    fprintf(fid, 'Thold  = %.6f\n', Thold);
    fprintf(fid, 'Tdown  = %.6f\n', Tdown);
    fprintf(fid, 'Tfinal = %.6f\n\n', Tfinal);

    fprintf(fid, 'Relative L2 difference in u1 final state: %.12e\n', rel_diff_u1(i));
    fprintf(fid, 'Relative L2 difference in u2 final state: %.12e\n', rel_diff_u2(i));
    fprintf(fid, 'Max absolute difference in u1 final state: %.12e\n', max_diff_u1(i));
    fprintf(fid, 'Max absolute difference in u2 final state: %.12e\n', max_diff_u2(i));
    fprintf(fid, 'Classification: %s\n', class_label(i));
    fclose(fid);
end

%% ============================================================
% 7. SUMMARY PLOTS
% ============================================================

% Relative differences vs bmax
fig = figure;
plot(bmax_list, rel_diff_u1, '-o', 'LineWidth', 1.8, 'MarkerSize', 7, 'DisplayName', 'u_1');
hold on;
plot(bmax_list, rel_diff_u2, '-s', 'LineWidth', 1.8, 'MarkerSize', 7, 'DisplayName', 'u_2');
yline(recovery_thresh, '--', 'Recovery threshold', 'LineWidth', 1.2);
yline(transition_thresh, '--', 'Transition threshold', 'LineWidth', 1.2);
xlabel('b_{max}');
ylabel('Relative final-state difference');
title('Refined threshold sweep: relative final-state difference vs b_{max}');
legend('Location','best');
box on;
hold off;
saveas(fig, fullfile(master_outdir, 'summary_rel_diff_vs_bmax.png'));
close(fig);

% Semilog plot for u2
fig = figure;
semilogy(bmax_list, rel_diff_u2, '-s', 'LineWidth', 1.8, 'MarkerSize', 7);
hold on;
yline(recovery_thresh, '--', 'Recovery threshold', 'LineWidth', 1.2);
yline(transition_thresh, '--', 'Transition threshold', 'LineWidth', 1.2);
xlabel('b_{max}');
ylabel('Relative final-state difference in u_2 (log scale)');
title('Refined threshold sweep: u_2 relative difference vs b_{max}');
box on;
hold off;
saveas(fig, fullfile(master_outdir, 'summary_rel_diff_u2_semilogy.png'));
close(fig);

% Max absolute differences vs bmax
fig = figure;
plot(bmax_list, max_diff_u1, '-o', 'LineWidth', 1.8, 'MarkerSize', 7, 'DisplayName', 'u_1');
hold on;
plot(bmax_list, max_diff_u2, '-s', 'LineWidth', 1.8, 'MarkerSize', 7, 'DisplayName', 'u_2');
xlabel('b_{max}');
ylabel('Max absolute final-state difference');
title('Refined threshold sweep: max absolute difference vs b_{max}');
legend('Location','best');
box on;
hold off;
saveas(fig, fullfile(master_outdir, 'summary_max_diff_vs_bmax.png'));
close(fig);

%% ============================================================
% 8. SAVE SUMMARY TABLE
% ============================================================

summary_table = table(bmax_list(:), rel_diff_u1, rel_diff_u2, max_diff_u1, max_diff_u2, class_label, ...
    'VariableNames', {'bmax','rel_diff_u1','rel_diff_u2','max_diff_u1','max_diff_u2','class_label'});

writetable(summary_table, fullfile(master_outdir, 'summary_metrics.csv'));
save(fullfile(master_outdir, 'summary_metrics.mat'), ...
    'bmax_list', 'rel_diff_u1', 'rel_diff_u2', 'max_diff_u1', 'max_diff_u2', 'class_label');

disp(' ');
disp('Refined threshold sweep complete.');
disp(['All outputs saved under folder: ', master_outdir]);
disp(summary_table);