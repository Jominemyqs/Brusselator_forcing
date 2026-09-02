clear; clc; close all;

%% ============================================================
% 1. CREATE OUTPUT FOLDER
% =============================================================
outdir = 'overshoot_compare_outputs';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

%% ============================================================
% 2. USER SETTINGS
% =============================================================

% Spatial grid
N  = 400;
Lx = 40;
x  = linspace(0, Lx, N)';

% Model parameter
par_forced.sigma  = 0.45;
par_control.sigma = 0.45;

% Baseline and overshoot parameters
b0    = 10.0;
bmax  = 13;
Tup   = 40;
Thold = 160;
Tdown = 40;

% Final time
Tfinal = 300;

% No internal plotting from solver
flag = 0;

%% ============================================================
% 3. DEFINE FORCING FUNCTIONS
% =============================================================

% Control: constant b(t) = b0
par_control.Bfun = @(t) b0;

% Forced: overshoot profile
par_forced.Bfun = @(t) overshoot_B(t, b0, bmax, Tup, Thold, Tdown);

%% ============================================================
% 4. INITIAL CONDITION
% =============================================================

% De Wit equilibrium at baseline b0
a     = 2.5;
u1_ss = a;
u2_ss = b0 / a;

rng(1);
u1 = u1_ss + 0.01 * randn(N,1);
u2 = u2_ss + 0.01 * randn(N,1);

ic0 = [x, u1, u2];

%% ============================================================
% 5. RUN CONTROL
% =============================================================

disp('Running control simulation...');
[ic_control, S_control, V_control] = solve_brusselator_1d_forced(ic0, par_control, Tfinal, flag);

u1_control_final = ic_control(:,2);
u2_control_final = ic_control(:,3);

%% ============================================================
% 6. RUN FORCED
% =============================================================

disp('Running forced simulation...');
[ic_forced, S_forced, V_forced] = solve_brusselator_1d_forced(ic0, par_forced, Tfinal, flag);

u1_forced_final = ic_forced(:,2);
u2_forced_final = ic_forced(:,3);

%% ============================================================
% 7. BUILD FORCING CURVES FOR PLOTTING
% =============================================================

tt = linspace(0, Tfinal, 1500);
b_control = arrayfun(par_control.Bfun, tt);
b_forced  = arrayfun(par_forced.Bfun, tt);

%% ============================================================
% 8. COMPUTE COMPARISON METRICS
% =============================================================

% Relative L2 differences in final states
rel_diff_u1 = norm(u1_forced_final - u1_control_final, 2) / max(norm(u1_control_final, 2), 1e-12);
rel_diff_u2 = norm(u2_forced_final - u2_control_final, 2) / max(norm(u2_control_final, 2), 1e-12);

% Max absolute differences
max_diff_u1 = max(abs(u1_forced_final - u1_control_final));
max_diff_u2 = max(abs(u2_forced_final - u2_control_final));

fprintf('\n=== Final-state comparison ===\n');
fprintf('Relative L2 difference in u1 final state: %.6e\n', rel_diff_u1);
fprintf('Relative L2 difference in u2 final state: %.6e\n', rel_diff_u2);
fprintf('Max absolute difference in u1 final state: %.6e\n', max_diff_u1);
fprintf('Max absolute difference in u2 final state: %.6e\n', max_diff_u2);

%% ============================================================
% 9. CHOOSE SNAPSHOT TIMES
% =============================================================

snapshot_times = [0, 20, 40, 80, 120, 160, 220, 300];

snapshot_inds_control = zeros(size(snapshot_times));
snapshot_inds_forced  = zeros(size(snapshot_times));

for k = 1:length(snapshot_times)
    [~, snapshot_inds_control(k)] = min(abs(S_control - snapshot_times(k)));
    [~, snapshot_inds_forced(k)]  = min(abs(S_forced  - snapshot_times(k)));
end

%% ============================================================
% 10. PLOT: FORCING COMPARISON
% =============================================================

fig = figure;
plot(tt, b_control, 'LineWidth', 1.8, 'DisplayName', 'Control');
hold on;
plot(tt, b_forced,  'LineWidth', 1.8, 'DisplayName', 'Forced');
xlabel('t');
ylabel('b(t)');
title('Control vs forced parameter profile');
legend('Location','best');
box on;
hold off;
saveas(fig, fullfile(outdir, 'forcing_comparison.png'));

%% ============================================================
% 11. PLOT: CONTROL SPACE-TIME
% =============================================================

fig = figure;
imagesc(1:N, S_control, V_control);
axis xy;
xlabel('Spatial index');
ylabel('Time');
title('Control: space-time plot of u_2(x,t)');
colorbar;
box on;
saveas(fig, fullfile(outdir, 'control_spacetime_u2.png'));

%% ============================================================
% 12. PLOT: FORCED SPACE-TIME
% =============================================================

fig = figure;
imagesc(1:N, S_forced, V_forced);
axis xy;
xlabel('Spatial index');
ylabel('Time');
title('Forced: space-time plot of u_2(x,t)');
colorbar;
box on;
saveas(fig, fullfile(outdir, 'forced_spacetime_u2.png'));

%% ============================================================
% 13. PLOT: FINAL PROFILE COMPARISON FOR u2
% =============================================================

fig = figure;
plot(x, u2_control_final, 'LineWidth', 1.8, 'DisplayName', 'Control final');
hold on;
plot(x, u2_forced_final,  '--', 'LineWidth', 1.8, 'DisplayName', 'Forced final');
xlabel('x');
ylabel('u_2(x,T)');
title(sprintf('Final profile comparison at T = %.1f', Tfinal));
legend('Location','best');
box on;
hold off;
saveas(fig, fullfile(outdir, 'final_profile_comparison_u2.png'));

%% ============================================================
% 14. PLOT: FINAL DIFFERENCE PROFILE FOR u2
% =============================================================

fig = figure;
plot(x, u2_forced_final - u2_control_final, 'LineWidth', 1.8);
xlabel('x');
ylabel('\Delta u_2 = u_{2,forced} - u_{2,control}');
title('Difference in final u_2 profiles');
box on;
saveas(fig, fullfile(outdir, 'final_difference_u2.png'));

%% ============================================================
% 15. PLOT: SNAPSHOT COMPARISON AT SELECTED TIMES
% =============================================================

fig = figure;
nrows = 2;
ncols = ceil(length(snapshot_times)/2);

for k = 1:length(snapshot_times)
    subplot(nrows, ncols, k);
    plot(x, V_control(snapshot_inds_control(k),:), 'LineWidth', 1.4, 'DisplayName', 'Control');
    hold on;
    plot(x, V_forced(snapshot_inds_forced(k),:), '--', 'LineWidth', 1.4, 'DisplayName', 'Forced');
    title(sprintf('t = %.1f', snapshot_times(k)));
    xlabel('x');
    ylabel('u_2');
    box on;
    hold off;
end

legend('Location','best');
sgtitle('Snapshot comparison: control vs forced');
saveas(fig, fullfile(outdir, 'snapshot_comparison_grid.png'));

%% ============================================================
% 16. SAVE NUMERICAL RESULTS
% =============================================================

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

results.params.sigma = par_forced.sigma;
results.params.b0    = b0;
results.params.bmax  = bmax;
results.params.Tup   = Tup;
results.params.Thold = Thold;
results.params.Tdown = Tdown;
results.params.Tfinal = Tfinal;

results.metrics.rel_diff_u1 = rel_diff_u1;
results.metrics.rel_diff_u2 = rel_diff_u2;
results.metrics.max_diff_u1 = max_diff_u1;
results.metrics.max_diff_u2 = max_diff_u2;

results.snapshot_times = snapshot_times;
results.snapshot_inds_control = snapshot_inds_control;
results.snapshot_inds_forced  = snapshot_inds_forced;

save(fullfile(outdir, 'control_forced_results.mat'), 'results');

%% ============================================================
% 17. SAVE A TEXT SUMMARY
% =============================================================

fid = fopen(fullfile(outdir, 'summary.txt'), 'w');
fprintf(fid, 'Control vs forced Brusselator comparison\n\n');
fprintf(fid, 'sigma  = %.6f\n', par_forced.sigma);
fprintf(fid, 'b0     = %.6f\n', b0);
fprintf(fid, 'bmax   = %.6f\n', bmax);
fprintf(fid, 'Tup    = %.6f\n', Tup);
fprintf(fid, 'Thold  = %.6f\n', Thold);
fprintf(fid, 'Tdown  = %.6f\n', Tdown);
fprintf(fid, 'Tfinal = %.6f\n\n', Tfinal);

fprintf(fid, 'Relative L2 difference in u1 final state: %.12e\n', rel_diff_u1);
fprintf(fid, 'Relative L2 difference in u2 final state: %.12e\n', rel_diff_u2);
fprintf(fid, 'Max absolute difference in u1 final state: %.12e\n', max_diff_u1);
fprintf(fid, 'Max absolute difference in u2 final state: %.12e\n', max_diff_u2);
fclose(fid);

disp('Done.');
disp(['All outputs saved in folder: ', outdir]);