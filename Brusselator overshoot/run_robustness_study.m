clear; clc; close all;

%RUN_ROBUSTNESS_STUDY Minimal sensitivity study around the reported bmax gap.
% New data are written below experiment_outputs/ and existing results are
% never overwritten. Each matrix point compares matched control and forced
% runs using the same initial condition.

study_id = 'robustness_reference_v2';
base_cfg = brusselator_default_config(struct( ...
    'experiment_name', study_id, ...
    'time', struct('post_forcing', 140)));

outdir = fullfile(base_cfg.output.root, study_id);
if isfolder(outdir)
    progress_file = fullfile(outdir, 'progress_summary.mat');
    if ~isfile(progress_file)
        error('run_robustness_study:OutputExists', ...
            'Existing directory has no resumable checkpoint: %s', outdir);
    end
    resume_data = load(progress_file, 'records', 'record_index');
    fprintf('Resuming %s after %d completed matched pairs.\n', ...
        study_id, resume_data.record_index);
else
    mkdir(outdir);
    mkdir(fullfile(outdir, 'raw_final_states'));
    brusselator_write_metadata(outdir, base_cfg, ...
        brusselator_run_metadata(base_cfg, mfilename, NaN));
end

% This screen tests spatial resolution, adaptive-solver tolerances, and
% seed sensitivity without changing the forcing protocol or recovery metric.
N_list = [200, 400, 800];
seed_list = [1, 2];
bmax_list = [11.12, 11.14];
solver_profiles = struct( ...
    'name', {'matlab_default', 'tight'}, ...
    'RelTol', {[], 1e-6}, ...
    'AbsTol', {[], 1e-9});

total_runs = numel(N_list) * numel(seed_list) * numel(bmax_list) * numel(solver_profiles);
empty_record = struct( ...
    'N', [], 'Lx', [], 'seed', [], 'solver_profile', '', ...
    'RelTol', [], 'AbsTol', [], 'bmax', [], 'forcing_end', [], ...
    'Tfinal', [], 'post_forcing_time', [], 'rel_diff_u', [], ...
    'rel_diff_v', [], 'max_diff_u', [], 'max_diff_v', [], ...
    'class_label', '', 'runtime_seconds', []);
if exist('resume_data', 'var')
    records = resume_data.records;
    record_index = resume_data.record_index;
    if numel(records) ~= total_runs || record_index < 0 || record_index > total_runs
        error('run_robustness_study:InvalidCheckpoint', ...
            'Checkpoint does not match the configured robustness matrix.');
    end
else
    records = repmat(empty_record, total_runs, 1);
    record_index = 0;
end
case_index = 0;

for iN = 1:numel(N_list)
    for iSeed = 1:numel(seed_list)
        for iSolver = 1:numel(solver_profiles)
            for iB = 1:numel(bmax_list)
                case_index = case_index + 1;
                if case_index <= record_index
                    fprintf('[%d/%d] already checkpointed; skipping.\n', ...
                        case_index, total_runs);
                    continue;
                end
                record_index = case_index;

                overrides = struct();
                overrides.grid = struct('N', N_list(iN));
                overrides.initial = struct('seed', seed_list(iSeed));
                overrides.forcing = struct('bmax', bmax_list(iB));
                overrides.solver = struct( ...
                    'RelTol', solver_profiles(iSolver).RelTol, ...
                    'AbsTol', solver_profiles(iSolver).AbsTol);
                cfg = brusselator_default_config(overrides);

                fprintf('\n[%d/%d] N=%d, seed=%d, solver=%s, bmax=%.2f\n', ...
                    record_index, total_runs, cfg.grid.N, cfg.initial.seed, ...
                    solver_profiles(iSolver).name, cfg.forcing.bmax);
                result = brusselator_run_pair(cfg);

                record = struct();
                record.N = cfg.grid.N;
                record.Lx = cfg.grid.Lx;
                record.seed = cfg.initial.seed;
                record.solver_profile = solver_profiles(iSolver).name;
                record.RelTol = cfg.solver.RelTol;
                record.AbsTol = cfg.solver.AbsTol;
                record.bmax = cfg.forcing.bmax;
                record.forcing_end = cfg.time.forcing_end;
                record.Tfinal = cfg.time.Tfinal;
                record.post_forcing_time = cfg.time.post_forcing;
                record.rel_diff_u = result.metrics.rel_diff_u;
                record.rel_diff_v = result.metrics.rel_diff_v;
                record.max_diff_u = result.metrics.max_diff_u;
                record.max_diff_v = result.metrics.max_diff_v;
                record.class_label = result.class_label;
                record.runtime_seconds = result.runtime_seconds;
                records(record_index) = record;

                % Retain raw final states and the exact IC, but not full
                % trajectories for every matrix point.
                raw = struct();
                raw.config = cfg;
                raw.ic0 = result.ic0;
                raw.control_final = result.control.ic_final;
                raw.forced_final = result.forced.ic_final;
                raw.metrics = result.metrics;
                raw.class_label = result.class_label;
                raw.runtime_seconds = result.runtime_seconds;
                raw_filename = sprintf('N_%d_seed_%d_%s_bmax_%0.2f.mat', ...
                    cfg.grid.N, cfg.initial.seed, solver_profiles(iSolver).name, cfg.forcing.bmax);
                save(fullfile(outdir, 'raw_final_states', raw_filename), 'raw');

                % A completed point remains inspectable if a later point
                % fails or MATLAB is interrupted.
                progress_table = struct2table(records(1:record_index), 'AsArray', true);
                writetable(progress_table, fullfile(outdir, 'progress_summary.csv'));
                save(fullfile(outdir, 'progress_summary.mat'), 'progress_table', ...
                    'records', 'record_index', 'N_list', 'seed_list', ...
                    'bmax_list', 'solver_profiles');
            end
        end
    end
end

summary_table = struct2table(records, 'AsArray', true);
writetable(summary_table, fullfile(outdir, 'robustness_summary.csv'));
save(fullfile(outdir, 'robustness_summary.mat'), 'summary_table', 'records', ...
    'N_list', 'seed_list', 'bmax_list', 'solver_profiles');

fig = figure('Color', 'w', 'Position', [100 100 900 550]);
hold on;
groups = unique(strcat("N=", string(summary_table.N), ", ", string(summary_table.solver_profile), ...
    ", seed=", string(summary_table.seed)));
for k = 1:numel(groups)
    rows = strcat("N=", string(summary_table.N), ", ", string(summary_table.solver_profile), ...
        ", seed=", string(summary_table.seed)) == groups(k);
    plot(summary_table.bmax(rows), summary_table.rel_diff_v(rows), '-o', ...
        'LineWidth', 1.3, 'DisplayName', groups(k));
end
yline(base_cfg.classification.recovery_threshold, '--', 'Recovery threshold');
yline(base_cfg.classification.transition_threshold, '--', 'Transition threshold');
set(gca, 'YScale', 'log');
xlabel('b_{max}');
ylabel('Relative final-state difference in v');
title('Robustness screen near the reported overshoot threshold');
legend('Location', 'eastoutside');
grid on; box on;
exportgraphics(fig, fullfile(outdir, 'robustness_rel_diff_v.png'), 'Resolution', 300);
close(fig);

metadata = brusselator_run_metadata(base_cfg, mfilename, sum(summary_table.runtime_seconds));
brusselator_write_metadata(outdir, base_cfg, metadata);

disp(summary_table);
disp(['Robustness study complete. Results saved in: ', outdir]);
