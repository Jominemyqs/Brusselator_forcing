clear; clc; close all;

%% ============================================================
% USER SETTINGS
% ============================================================

sigma = 0.45;

% Target baseline
b0 = 10.2;

% Overshoot protocol
Tup    = 40;
Thold  = 80;
Tdown  = 40;
Tfinal_forced = 300;

% Baseline/frozen evolution time for lambda tests
Tfinal_lambda = 2000;

% Same overshoot amplitude used in previous map
delta_b = 1.2;
bmax = b0 + delta_b;

% Refined lambda range for b0 = 10.2
lambda_list = 0.00:0.05:0.65;

% Threshold for returning to A
tol_A = 1e-2;

flag = 0;

% Output folder
outdir = 'refine_b0_10p2_outputs';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

%% ============================================================
% GRID AND INITIAL CONDITION
% ============================================================

N  = 400;
Lx = 40;
x  = linspace(0, Lx, N)';

a = 2.5;

u1_ss = a;
u2_ss = b0 / a;

rng(1);
u1 = u1_ss + 0.01 * randn(N,1);
u2 = u2_ss + 0.01 * randn(N,1);
ic0 = [x, u1, u2];

%% ============================================================
% COMPUTE BASELINE STATE A
% ============================================================

fprintf('\n=================================================\n');
fprintf('Computing baseline state A for b0 = %.3f\n', b0);
fprintf('=================================================\n');

par_control.sigma = sigma;
par_control.Bfun  = @(t) b0;

[icA, SA, VA] = solve_brusselator_1d_forced(ic0, par_control, Tfinal_lambda, flag);

uA = icA(:,2);
vA = icA(:,3);

%% ============================================================
% COMPUTE SUPERCRITICAL OVERSHOOT STATE B
% ============================================================

fprintf('\n=================================================\n');
fprintf('Computing overshoot state B for b0 = %.3f, bmax = %.3f\n', b0, bmax);
fprintf('=================================================\n');

par_forced.sigma = sigma;
par_forced.Bfun  = @(t) overshoot_B(t, b0, bmax, Tup, Thold, Tdown);

[icB, SB, VB] = solve_brusselator_1d_forced(ic0, par_forced, Tfinal_forced, flag);

uB = icB(:,2);
vB = icB(:,3);

%% ============================================================
% SCAN LAMBDA
% ============================================================

nL = length(lambda_list);

distA = zeros(nL,1);
return_flag = false(nL,1);
v_final_store = cell(nL,1);

for il = 1:nL
    lam = lambda_list(il);

    fprintf('\n-------------------------------------\n');
    fprintf('Testing lambda = %.3f\n', lam);
    fprintf('-------------------------------------\n');

    u0_lam = (1-lam)*uA + lam*uB;
    v0_lam = (1-lam)*vA + lam*vB;
    ic_lam = [x, u0_lam, v0_lam];

    [ic_final, ~, ~] = solve_brusselator_1d_forced(ic_lam, par_control, Tfinal_lambda, flag);

    v_final = ic_final(:,3);
    v_final_store{il} = v_final;

    distA(il) = norm(v_final - vA, 2) / max(norm(vA,2), 1e-12);
    return_flag(il) = distA(il) < tol_A;

    fprintf('lambda = %.3f, distA = %.6e, return = %d\n', ...
        lam, distA(il), return_flag(il));
end

%% ============================================================
% EXTRACT THRESHOLD INTERVAL
% ============================================================

idx_return = find(return_flag == 1, 1, 'last');
idx_leave  = find(return_flag == 0, 1, 'first');

lambda_threshold_low  = NaN;
lambda_threshold_high = NaN;

if ~isempty(idx_return)
    lambda_threshold_low = lambda_list(idx_return);
end

if ~isempty(idx_leave)
    lambda_threshold_high = lambda_list(idx_leave);
end

summary_table = table(lambda_list(:), distA, return_flag, ...
    'VariableNames', {'lambda','dist_to_A','returns_to_A'});

threshold_table = table(b0, bmax, lambda_threshold_low, lambda_threshold_high, ...
    'VariableNames', {'b0','bmax','lambda_threshold_low','lambda_threshold_high'});

disp(' ');
disp('Lambda scan summary:');
disp(summary_table);

disp(' ');
disp('Threshold interval:');
disp(threshold_table);

%% ============================================================
% PLOTS
% ============================================================

fig = figure;
plot(lambda_list, distA, '-o', 'LineWidth', 1.8, 'MarkerSize', 7);
yline(tol_A, '--', 'return threshold', 'LineWidth', 1.2);
xlabel('\lambda');
ylabel('Relative distance to A');
title(sprintf('Refined basin-slice scan for b_0 = %.1f', b0));
grid on;
exportgraphics(fig, fullfile(outdir, 'dist_to_A_vs_lambda_b0_10p2.png'));

fig = figure;
stairs(lambda_list, double(return_flag), 'LineWidth', 2);
xlabel('\lambda');
ylabel('Return to A?');
yticks([0 1]);
yticklabels({'No','Yes'});
title(sprintf('Return classification for b_0 = %.1f', b0));
grid on;
exportgraphics(fig, fullfile(outdir, 'return_flag_vs_lambda_b0_10p2.png'));

%% ============================================================
% SAVE OUTPUT
% ============================================================

writetable(summary_table, fullfile(outdir, 'lambda_scan_summary_b0_10p2.csv'));
writetable(threshold_table, fullfile(outdir, 'threshold_interval_b0_10p2.csv'));

save(fullfile(outdir, 'refine_b0_10p2_results.mat'), ...
    'b0', 'bmax', 'sigma', 'lambda_list', 'distA', 'return_flag', ...
    'lambda_threshold_low', 'lambda_threshold_high', ...
    'x', 'uA', 'vA', 'uB', 'vB', 'v_final_store', ...
    'Tup', 'Thold', 'Tdown', 'Tfinal_forced', 'Tfinal_lambda');

disp(' ');
disp(['Done. Outputs saved in: ', outdir]);