clear; clc; close all;

%% ============================================================
% 1. MASTER OUTPUT FOLDER
% ============================================================
master_outdir = 'Thold_bmax_sweep_fast_outputs';
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
Tdown  = 40;
Tfinal = 300;

% Sweep values
Thold_list = [20, 80, 320];
bmax_list  = [11.08, 11.10, 11.12, 11.14, 11.16];

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

save(fullfile(master_outdir, 'control_results.mat'), ...
    'x', 'ic0', 'ic_control', 'S_control', 'V_control', ...
    'u1_control_final', 'u2_control_final', ...
    'sigma', 'b0', 'Tup', 'Tdown', 'Tfinal');

%% ============================================================
% 5. PREALLOCATE SUMMARY STORAGE
% ============================================================
nH = length(Thold_list);
nB = length(bmax_list);

rel_diff_u1 = zeros(nH, nB);
rel_diff_u2 = zeros(nH, nB);
max_diff_u1 = zeros(nH, nB);
max_diff_u2 = zeros(nH, nB);
class_label = strings(nH, nB);

% Optional: store final profiles only
u2_final_store = cell(nH, nB);

%% ============================================================
% 6. LOOP OVER THOLD AND BMAX
% ============================================================
for iH = 1:nH
    Thold = Thold_list(iH);

    for iB = 1:nB
        bmax = bmax_list(iB);

        fprintf('\n====================================\n');
        fprintf('Running simulation for Thold = %.1f, bmax = %.3f\n', Thold, bmax);
        fprintf('====================================\n');

        par_forced.sigma = sigma;
        par_forced.Bfun  = @(t) overshoot_B(t, b0, bmax, Tup, Thold, Tdown);

        [ic_forced, S_forced, V_forced] = solve_brusselator_1d_forced(ic0, par_forced, Tfinal, flag);

        u1_forced_final = ic_forced(:,2);
        u2_forced_final = ic_forced(:,3);

        rel_diff_u1(iH,iB) = norm(u1_forced_final - u1_control_final, 2) / max(norm(u1_control_final, 2), 1e-12);
        rel_diff_u2(iH,iB) = norm(u2_forced_final - u2_control_final, 2) / max(norm(u2_control_final, 2), 1e-12);
        max_diff_u1(iH,iB) = max(abs(u1_forced_final - u1_control_final));
        max_diff_u2(iH,iB) = max(abs(u2_forced_final - u2_control_final));

        if rel_diff_u2(iH,iB) < recovery_thresh
            class_label(iH,iB) = "Recovery";
        elseif rel_diff_u2(iH,iB) > transition_thresh
            class_label(iH,iB) = "LastingTransition";
        else
            class_label(iH,iB) = "Intermediate";
        end

        u2_final_store{iH,iB} = u2_forced_final;

        fprintf('Relative L2 difference in u2 final state: %.6e\n', rel_diff_u2(iH,iB));
        fprintf('Classified as: %s\n', class_label(iH,iB));

        % Save minimal per-run data
        run_data = struct();
        run_data.Thold = Thold;
        run_data.bmax  = bmax;
        run_data.S_forced = S_forced;
        run_data.ic_forced = ic_forced;
        run_data.rel_diff_u1 = rel_diff_u1(iH,iB);
        run_data.rel_diff_u2 = rel_diff_u2(iH,iB);
        run_data.max_diff_u1 = max_diff_u1(iH,iB);
        run_data.max_diff_u2 = max_diff_u2(iH,iB);
        run_data.class_label = class_label(iH,iB);

        save(fullfile(master_outdir, sprintf('run_Thold_%d_bmax_%0.2f.mat', Thold, bmax)), 'run_data');
    end
end

%% ============================================================
% 7. SUMMARY PLOTS ONLY
% ============================================================

% Relative difference vs bmax, one curve per Thold
fig = figure('Visible','off');
hold on;
for iH = 1:nH
    plot(bmax_list, rel_diff_u2(iH,:), '-o', 'LineWidth', 1.8, 'MarkerSize', 7, ...
        'DisplayName', sprintf('Thold = %d', Thold_list(iH)));
end
yline(recovery_thresh, '--', 'Recovery threshold', 'LineWidth', 1.2);
yline(transition_thresh, '--', 'Transition threshold', 'LineWidth', 1.2);
xlabel('b_{max}');
ylabel('Relative final-state difference in u_2');
title('Threshold shift with hold time');
legend('Location','best');
box on;
hold off;
exportgraphics(fig, fullfile(master_outdir, 'summary_rel_diff_u2_vs_bmax_by_Thold.png'));
close(fig);

% Semilog version
fig = figure('Visible','off');
hold on;
for iH = 1:nH
    semilogy(bmax_list, rel_diff_u2(iH,:), '-o', 'LineWidth', 1.8, 'MarkerSize', 7, ...
        'DisplayName', sprintf('Thold = %d', Thold_list(iH)));
end
yline(recovery_thresh, '--', 'Recovery threshold', 'LineWidth', 1.2);
yline(transition_thresh, '--', 'Transition threshold', 'LineWidth', 1.2);
xlabel('b_{max}');
ylabel('Relative final-state difference in u_2 (log scale)');
title('Threshold shift with hold time (semilog)');
legend('Location','best');
box on;
hold off;
exportgraphics(fig, fullfile(master_outdir, 'summary_rel_diff_u2_semilog_by_Thold.png'));
close(fig);

% Heatmap
fig = figure('Visible','off');
imagesc(bmax_list, Thold_list, rel_diff_u2);
axis xy;
xlabel('b_{max}');
ylabel('T_{hold}');
title('Heatmap of relative final-state difference in u_2');
colorbar;
exportgraphics(fig, fullfile(master_outdir, 'heatmap_rel_diff_u2.png'));
close(fig);

%% ============================================================
% 8. SAVE SUMMARY TABLE
% ============================================================
rows = [];
for iH = 1:nH
    for iB = 1:nB
        rows = [rows; ...
            Thold_list(iH), ...
            bmax_list(iB), ...
            rel_diff_u1(iH,iB), ...
            rel_diff_u2(iH,iB), ...
            max_diff_u1(iH,iB), ...
            max_diff_u2(iH,iB)];
    end
end

class_col = reshape(class_label.', [], 1);

summary_table = table( ...
    rows(:,1), rows(:,2), rows(:,3), rows(:,4), rows(:,5), rows(:,6), class_col, ...
    'VariableNames', {'Thold','bmax','rel_diff_u1','rel_diff_u2','max_diff_u1','max_diff_u2','class_label'});

writetable(summary_table, fullfile(master_outdir, 'summary_metrics.csv'));
save(fullfile(master_outdir, 'summary_metrics.mat'), ...
    'Thold_list', 'bmax_list', 'rel_diff_u1', 'rel_diff_u2', 'max_diff_u1', 'max_diff_u2', 'class_label', ...
    'u2_final_store', 'u2_control_final', 'x');

disp(' ');
disp('Fast Thold-bmax sweep complete.');
disp(['All outputs saved under folder: ', master_outdir]);
disp(summary_table);