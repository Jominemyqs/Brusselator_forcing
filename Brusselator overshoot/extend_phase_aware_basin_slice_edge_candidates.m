function results = extend_phase_aware_basin_slice_edge_candidates()
%EXTEND_PHASE_AWARE_BASIN_SLICE_EDGE_CANDIDATES Continue unresolved slice cases.
%   Loads the versioned phase-aware basin-slice pilot and continues only its
%   unresolved coarse trajectories under the same frozen b=10 system.  This
%   distinguishes a finite-time transitional band from trajectories that
%   eventually approach one of the two already-defined endpoint outcomes.
%   The clear baseline and periodic-orbit cases are reused from saved output;
%   they are deliberately not recomputed or reclassified.

study_id = 'phase_aware_basin_slice_edge_extension_v1';
parent_file = fullfile('experiment_outputs', 'phase_aware_basin_slice_v1', ...
    'raw_phase_aware_basin_slice.mat');
outdir = fullfile('experiment_outputs', study_id);
if ~isfile(parent_file)
    error('extend_phase_aware_basin_slice_edge_candidates:MissingParentData', ...
        'Expected phase-aware slice pilot data at %s.', parent_file);
end
if isfolder(outdir)
    error('extend_phase_aware_basin_slice_edge_candidates:OutputExists', ...
        'Refusing to overwrite existing experiment directory: %s', outdir);
end

timer = tic;
loaded = load(parent_file, 'results');
parent = loaded.results;
cfg = parent.configuration;
cfg.experiment_name = study_id;
cfg.edge_extension = struct( ...
    'parent_file', parent_file, ...
    'selection', 'coarse cases unresolved within the initial 80-unit window', ...
    'original_duration', cfg.basin_slice.coarse_duration, ...
    'continuation_duration', cfg.basin_slice.edge_extension_duration - ...
        cfg.basin_slice.coarse_duration, ...
    'total_duration', cfg.basin_slice.edge_extension_duration, ...
    'classification_metric', cfg.basin_slice.metric, ...
    'late_window', cfg.basin_slice.late_window);
if cfg.edge_extension.continuation_duration <= 0
    error('extend_phase_aware_basin_slice_edge_candidates:InvalidDuration', ...
        'The configured edge-extension duration must exceed the coarse duration.');
end

mkdir(outdir);
brusselator_write_metadata(outdir, cfg, ...
    brusselator_run_metadata(cfg, mfilename, NaN));

candidate_indices = find(strcmp({parent.records.stage}, 'coarse') & ...
    strcmp({parent.records.classification}, 'unresolved_within_test_window'));
if isempty(candidate_indices)
    error('extend_phase_aware_basin_slice_edge_candidates:NoCandidates', ...
        'The parent pilot has no unresolved coarse cases to extend.');
end

frozen_par = brusselator_make_parameters(cfg, @(t) cfg.basin_slice.frozen_b);
record_template = empty_record_template();
records = repmat(record_template, numel(candidate_indices), 1);
trajectory_data = cell(numel(candidate_indices), 1);
for k = 1:numel(candidate_indices)
    parent_index = candidate_indices(k);
    parent_record = parent.records(parent_index);
    old_data = parent.trajectory_data{parent_index};
    fprintf(['Extending parent lambda = %.8f from t = %.0f to t = %.0f ', ...
        'under frozen b = %.2f ...\n'], parent_record.lambda, ...
        old_data.S(end), cfg.edge_extension.total_duration, cfg.basin_slice.frozen_b);
    case_timer = tic;
    [final_state, S_new, V_new, U_new] = solve_brusselator_1d_forced( ...
        old_data.final_state, frozen_par, cfg.edge_extension.continuation_duration, 0);
    runtime = toc(case_timer);

    [new_sample_times, new_dist_A, new_dist_orbit, new_A_reflection, ...
        new_orbit_lag, new_orbit_reflection] = sample_distances( ...
        S_new, U_new, V_new, parent.endpoint_A, parent.template, cfg);
    absolute_new_times = old_data.S(end) + new_sample_times;
    late = absolute_new_times >= cfg.edge_extension.total_duration - ...
        cfg.basin_slice.late_window;
    if ~any(late)
        error('extend_phase_aware_basin_slice_edge_candidates:LateWindowEmpty', ...
            'No samples lie in the requested final comparison window.');
    end

    record = record_template;
    record.parent_case_index = parent_index;
    record.lambda = parent_record.lambda;
    record.original_duration = old_data.S(end);
    record.continuation_duration = cfg.edge_extension.continuation_duration;
    record.total_duration = cfg.edge_extension.total_duration;
    record.start_distance_to_A = new_dist_A(1);
    record.start_distance_to_orbit = new_dist_orbit(1);
    record.final_distance_to_A = new_dist_A(end);
    record.final_distance_to_orbit = new_dist_orbit(end);
    record.late_median_distance_to_A = median(new_dist_A(late));
    record.late_maximum_distance_to_A = max(new_dist_A(late));
    record.late_median_distance_to_orbit = median(new_dist_orbit(late));
    record.late_maximum_distance_to_orbit = max(new_dist_orbit(late));
    record.final_A_uses_reflection = new_A_reflection(end);
    record.final_orbit_uses_reflection = new_orbit_reflection(end);
    record.final_orbit_phase_lag = new_orbit_lag(end);
    record.runtime_seconds = runtime;
    record.classification = classify_record(record, cfg);
    records(k) = record;

    trajectory_data{k} = struct( ...
        'parent_case_index', parent_index, 'old_S', old_data.S, ...
        'old_U', old_data.U, 'old_V', old_data.V, ...
        'old_sample_times', old_data.sample_times, ...
        'old_distance_to_A', old_data.distance_to_A, ...
        'old_distance_to_orbit', old_data.distance_to_orbit, ...
        'continuation_S', S_new, 'continuation_U', U_new, ...
        'continuation_V', V_new, 'continuation_final_state', final_state, ...
        'continuation_absolute_sample_times', absolute_new_times, ...
        'continuation_distance_to_A', new_dist_A, ...
        'continuation_distance_to_orbit', new_dist_orbit, ...
        'continuation_A_uses_reflection', new_A_reflection, ...
        'continuation_orbit_phase_lag', new_orbit_lag, ...
        'continuation_orbit_uses_reflection', new_orbit_reflection);
    write_checkpoint(outdir, cfg, parent, records, trajectory_data);
    fprintf('  %s; final-window medians d_A = %.3e, d_orbit = %.3e (%.1f s).\n', ...
        record.classification, record.late_median_distance_to_A, ...
        record.late_median_distance_to_orbit, runtime);
end

summary_table = struct2table(records, 'AsArray', true);
writetable(summary_table, fullfile(outdir, 'edge_extension_summary.csv'));
make_extension_figure(outdir, records, trajectory_data, cfg);

results = struct();
results.configuration = cfg;
results.parent_file = parent_file;
results.parent_summary = parent.summary;
results.records = records;
results.trajectory_data = trajectory_data;
save(fullfile(outdir, 'raw_edge_extension_trajectories.mat'), 'results', '-v7.3');

metadata = brusselator_run_metadata(cfg, mfilename, toc(timer));
brusselator_write_metadata(outdir, cfg, metadata);
disp(summary_table);
fprintf('Edge-candidate extension complete. Results saved in: %s\n', outdir);
end

function record = empty_record_template()
record = struct( ...
    'parent_case_index', NaN, ...
    'lambda', NaN, ...
    'original_duration', NaN, ...
    'continuation_duration', NaN, ...
    'total_duration', NaN, ...
    'start_distance_to_A', NaN, ...
    'start_distance_to_orbit', NaN, ...
    'final_distance_to_A', NaN, ...
    'final_distance_to_orbit', NaN, ...
    'late_median_distance_to_A', NaN, ...
    'late_maximum_distance_to_A', NaN, ...
    'late_median_distance_to_orbit', NaN, ...
    'late_maximum_distance_to_orbit', NaN, ...
    'final_A_uses_reflection', false, ...
    'final_orbit_uses_reflection', false, ...
    'final_orbit_phase_lag', NaN, ...
    'runtime_seconds', NaN, ...
    'classification', 'unclassified');
end

function [sample_times, dist_A, dist_orbit, A_reflection, orbit_lag, ...
        orbit_reflection] = sample_distances(S, U, V, A, template, cfg)
stride = max(1, round(cfg.basin_slice.distance_sample_interval / ...
    cfg.solver.output_dt));
sample_indices = unique([1:stride:numel(S), numel(S)]);
sample_times = S(sample_indices);
dist_A = zeros(numel(sample_indices), 1);
dist_orbit = zeros(numel(sample_indices), 1);
A_reflection = false(numel(sample_indices), 1);
orbit_lag = zeros(numel(sample_indices), 1);
orbit_reflection = false(numel(sample_indices), 1);
for j = 1:numel(sample_indices)
    [dist_A(j), A_reflection(j)] = distance_to_baseline( ...
        U(sample_indices(j),:), V(sample_indices(j),:), A);
    [dist_orbit(j), orbit_lag(j), orbit_reflection(j)] = ...
        distance_to_orbit(U(sample_indices(j),:), V(sample_indices(j),:), template);
end
end

function [distance, uses_reflection] = distance_to_baseline(u, v, A)
direct = sqrt(sum((u - A.U).^2) + sum((v - A.V).^2));
reflected = sqrt(sum((u - fliplr(A.U)).^2) + ...
    sum((v - fliplr(A.V)).^2));
normalizer = max(sqrt(sum(A.U.^2) + sum(A.V.^2)), eps);
if reflected < direct
    distance = reflected / normalizer;
    uses_reflection = true;
else
    distance = direct / normalizer;
    uses_reflection = false;
end
end

function [distance, phase_lag, uses_reflection] = distance_to_orbit(u, v, template)
direct = sqrt(sum((template.U - u).^2 + (template.V - v).^2, 2)) ./ ...
    max(sqrt(sum(template.U.^2 + template.V.^2, 2)), eps);
U_reflected = fliplr(template.U);
V_reflected = fliplr(template.V);
reflected = sqrt(sum((U_reflected - u).^2 + (V_reflected - v).^2, 2)) ./ ...
    max(sqrt(sum(U_reflected.^2 + V_reflected.^2, 2)), eps);
[direct_minimum, direct_index] = min(direct);
[reflection_minimum, reflection_index] = min(reflected);
if reflection_minimum < direct_minimum
    distance = reflection_minimum;
    phase_lag = template.phase_times(reflection_index);
    uses_reflection = true;
else
    distance = direct_minimum;
    phase_lag = template.phase_times(direct_index);
    uses_reflection = false;
end
end

function label = classify_record(record, cfg)
near_A = record.late_median_distance_to_A < cfg.classification.recovery_threshold && ...
    record.late_maximum_distance_to_A < cfg.basin_slice.late_max_threshold;
near_orbit = record.late_median_distance_to_orbit < ...
    cfg.classification.recovery_threshold && ...
    record.late_maximum_distance_to_orbit < cfg.basin_slice.late_max_threshold;
if near_A && ~near_orbit
    label = 'returns_to_baseline_neighborhood';
elseif near_orbit && ~near_A
    label = 'returns_to_periodic_orbit_neighborhood';
elseif near_A && near_orbit
    label = 'ambiguous_both_endpoint_neighborhoods';
else
    label = 'persistent_unresolved_edge_candidate';
end
end

function write_checkpoint(outdir, cfg, parent, records, trajectory_data)
writetable(struct2table(records, 'AsArray', true), ...
    fullfile(outdir, 'progress_summary.csv'));
save(fullfile(outdir, 'progress_edge_extension.mat'), 'cfg', 'parent', ...
    'records', 'trajectory_data', '-v7.3');
end

function make_extension_figure(outdir, records, trajectories, cfg)
fig = figure('Color', 'w', 'Position', [100 100 1000 650]);
tiledlayout(2, 1, 'TileSpacing', 'compact');
nexttile; hold on;
for k = 1:numel(records)
    data = trajectories{k};
    times = [data.old_sample_times; data.continuation_absolute_sample_times(2:end)];
    values = [data.old_distance_to_A; data.continuation_distance_to_A(2:end)];
    semilogy(times, max(values, eps), 'LineWidth', 1.15, ...
        'DisplayName', sprintf('\\lambda=%.5f (%s)', records(k).lambda, ...
        records(k).classification));
end
yline(cfg.classification.recovery_threshold, '--', 'HandleVisibility', 'off');
ylabel('d_A'); title('Extended unresolved slice trajectories: baseline distance');
legend('Location', 'eastoutside', 'FontSize', 7); grid on; box on;
nexttile; hold on;
for k = 1:numel(records)
    data = trajectories{k};
    times = [data.old_sample_times; data.continuation_absolute_sample_times(2:end)];
    values = [data.old_distance_to_orbit; data.continuation_distance_to_orbit(2:end)];
    semilogy(times, max(values, eps), 'LineWidth', 1.15, 'HandleVisibility', 'off');
end
yline(cfg.classification.recovery_threshold, '--', 'HandleVisibility', 'off');
xlabel('Frozen b=10 time'); ylabel('d_{orbit}');
title('Distance to candidate periodic orbit modulo phase/reflection');
grid on; box on;
exportgraphics(fig, fullfile(outdir, 'extended_edge_candidate_distances.png'), ...
    'Resolution', 300);
close(fig);
end
