clear; clc; close all;

%RUN_EXTENDED_POSTFORCING_TRAJECTORIES Extend the two bracketing trajectories.
% The selected cases retain the historical forcing protocol, then integrate
% for 10,000 time units after b has returned to b0. Direct distances are to
% the matched control trajectory at the same time; no phase/translation
% alignment is applied in this first extension. A tight stiff solve is used
% because explicit ode45 is prohibitively slow at this time horizon.

study_id = 'extended_postforcing_reference_v2';
base_cfg = brusselator_default_config(struct( ...
    'experiment_name', study_id, ...
    'time', struct('post_forcing', 10000), ...
    'solver', struct('name', 'ode15s', 'output_dt', 1.0, ...
        'RelTol', 1e-6, 'AbsTol', 1e-9, 'MaxStep', 1.0)));

outdir = fullfile(base_cfg.output.root, study_id);
if isfolder(outdir)
    error('run_extended_postforcing_trajectories:OutputExists', ...
        'Refusing to overwrite existing experiment directory: %s', outdir);
end
mkdir(outdir);
brusselator_write_metadata(outdir, base_cfg, ...
    brusselator_run_metadata(base_cfg, mfilename, NaN));

[x, ic0] = brusselator_initial_condition(base_cfg);
control_par = brusselator_make_parameters(base_cfg, @(t) base_cfg.forcing.b0);

fprintf('Running matched control to t = %.0f ...\n', base_cfg.time.Tfinal);
control_timer = tic;
[ic_control, S_control, V_control] = solve_brusselator_1d_forced( ...
    ic0, control_par, base_cfg.time.Tfinal, 0);
control_runtime_seconds = toc(control_timer);
save(fullfile(outdir, 'progress_control.mat'), 'base_cfg', 'ic0', ...
    'S_control', 'V_control', 'ic_control', 'control_runtime_seconds', '-v7.3');

bmax_list = [11.12, 11.14];
empty_record = struct( ...
    'bmax', [], 'forcing_end', [], 'Tfinal', [], 'post_forcing_time', [], ...
    'final_rel_diff_v', [], 'min_postforcing_rel_diff_v', [], ...
    'median_last_1000_rel_diff_v', [], 'runtime_seconds', []);
case_records = repmat(empty_record, numel(bmax_list), 1);
forced_trajectories = cell(numel(bmax_list), 1);
distance_to_control = cell(numel(bmax_list), 1);
case_configs = cell(numel(bmax_list), 1);

for k = 1:numel(bmax_list)
    cfg = brusselator_default_config(struct( ...
        'experiment_name', study_id, ...
        'forcing', struct('bmax', bmax_list(k)), ...
        'time', struct('post_forcing', base_cfg.time.post_forcing), ...
        'solver', struct('output_dt', base_cfg.solver.output_dt)));
    forced_par = brusselator_make_parameters(cfg, @(t) overshoot_B(t, ...
        cfg.forcing.b0, cfg.forcing.bmax, cfg.forcing.Tup, ...
        cfg.forcing.Thold, cfg.forcing.Tdown));

    fprintf('Running bmax = %.2f to t = %.0f ...\n', ...
        cfg.forcing.bmax, cfg.time.Tfinal);
    case_timer = tic;
    [ic_forced, S_forced, V_forced] = solve_brusselator_1d_forced( ...
        ic0, forced_par, cfg.time.Tfinal, 0);
    runtime_seconds = toc(case_timer);

    if ~isequal(S_control, S_forced)
        error('run_extended_postforcing_trajectories:TimeMismatch', ...
            'Control and forced output times must match.');
    end
    distance = sqrt(sum((V_forced - V_control).^2, 2)) ./ ...
        max(sqrt(sum(V_control.^2, 2)), 1e-12);
    post_idx = S_forced >= cfg.time.forcing_end;

    record = struct();
    record.bmax = cfg.forcing.bmax;
    record.forcing_end = cfg.time.forcing_end;
    record.Tfinal = cfg.time.Tfinal;
    record.post_forcing_time = cfg.time.post_forcing;
    record.final_rel_diff_v = distance(end);
    record.min_postforcing_rel_diff_v = min(distance(post_idx));
    record.median_last_1000_rel_diff_v = median(distance(S_forced >= cfg.time.Tfinal - 1000));
    record.runtime_seconds = runtime_seconds;
    case_records(k) = record;

    forced_trajectories{k} = struct('bmax', cfg.forcing.bmax, ...
        'S', S_forced, 'V', V_forced, 'ic_final', ic_forced);
    distance_to_control{k} = distance;
    case_configs{k} = cfg;

    progress_table = struct2table(case_records(1:k), 'AsArray', true);
    writetable(progress_table, fullfile(outdir, 'progress_summary.csv'));
    save(fullfile(outdir, 'progress_forced_cases.mat'), 'case_configs', ...
        'forced_trajectories', 'distance_to_control', 'case_records', '-v7.3');
end

raw = struct();
raw.config = base_cfg;
raw.ic0 = ic0;
raw.control = struct('S', S_control, 'V', V_control, 'ic_final', ic_control, ...
    'runtime_seconds', control_runtime_seconds);
raw.case_configs = case_configs;
raw.forced_trajectories = forced_trajectories;
raw.distance_to_control = distance_to_control;
save(fullfile(outdir, 'raw_extended_trajectories.mat'), 'raw', '-v7.3');

summary_table = struct2table(case_records, 'AsArray', true);
writetable(summary_table, fullfile(outdir, 'extended_trajectory_summary.csv'));
plot_extended_postforcing_results(outdir);

metadata = brusselator_run_metadata(base_cfg, mfilename, ...
    control_runtime_seconds + sum(summary_table.runtime_seconds));
brusselator_write_metadata(outdir, base_cfg, metadata);

disp(summary_table);
disp(['Extended trajectories complete. Results saved in: ', outdir]);
