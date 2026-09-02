function results = analyze_local_stability_orbit_distance(input_dir, output_dir)
%ANALYZE_LOCAL_STABILITY_ORBIT_DISTANCE Refine orbit distances after a run.
%   Re-evaluates saved perturbation trajectories against a densely interpolated
%   version of the reference periodic orbit. This removes the O(0.05) phase
%   sampling floor of the initial local-stability runner without rerunning any
%   PDE trajectories.

if nargin < 1 || isempty(input_dir)
    input_dir = fullfile('experiment_outputs', 'periodic_orbit_local_stability_v1');
end
if nargin < 2 || isempty(output_dir)
    output_dir = fullfile('experiment_outputs', 'periodic_orbit_local_stability_reanalysis_v1');
end
input_file = fullfile(input_dir, 'raw_local_stability_trajectories.mat');
if ~isfile(input_file)
    error('analyze_local_stability_orbit_distance:MissingInput', ...
        'Expected saved local-stability trajectories at %s.', input_file);
end
if isfolder(output_dir)
    error('analyze_local_stability_orbit_distance:OutputExists', ...
        'Refusing to overwrite existing analysis directory: %s', output_dir);
end
mkdir(output_dir);

timer = tic;
loaded = load(input_file, 'results');
parent = loaded.results;
cfg = parent.configuration;
phase_step = 0.0025;
template = densely_interpolate_template(parent.template, ...
    cfg.local_stability.observed_direct_period, phase_step);

record_template = struct( ...
    'phase_fraction', NaN, 'phase_start_time', NaN, 'direction', '', ...
    'amplitude_relative_full_state', NaN, 'initial_orbit_distance', NaN, ...
    'final_orbit_distance', NaN, 'median_last_10_orbit_distance', NaN, ...
    'maximum_last_10_orbit_distance', NaN, 'best_final_lag', NaN, ...
    'best_final_reflection', false, 'classification', 'unclassified');
records = repmat(record_template, numel(parent.records), 1);
distance_histories = cell(numel(records), 1);

for k = 1:numel(records)
    data = parent.trajectory_data{k};
    sample_indices = nearest_time_indices(data.S, data.sample_times);
    distances = zeros(numel(sample_indices), 1);
    lags = zeros(numel(sample_indices), 1);
    reflections = false(numel(sample_indices), 1);
    for j = 1:numel(sample_indices)
        [distances(j), lags(j), reflections(j)] = dense_template_distance( ...
            data.U(sample_indices(j),:), data.V(sample_indices(j),:), template);
    end

    original = parent.records(k);
    record = record_template;
    record.phase_fraction = original.phase_fraction;
    record.phase_start_time = original.phase_start_time;
    record.direction = original.direction;
    record.amplitude_relative_full_state = original.amplitude_relative_full_state;
    record.initial_orbit_distance = distances(1);
    record.final_orbit_distance = distances(end);
    late = data.sample_times >= data.S(end) - 10;
    record.median_last_10_orbit_distance = median(distances(late));
    record.maximum_last_10_orbit_distance = max(distances(late));
    record.best_final_lag = lags(end);
    record.best_final_reflection = reflections(end);
    record.classification = classify_record(record, cfg.classification);
    records(k) = record;
    distance_histories{k} = struct('sample_times', data.sample_times, ...
        'orbit_distance', distances, 'best_lag', lags, ...
        'best_reflection', reflections);
end

summary_table = struct2table(records, 'AsArray', true);
writetable(summary_table, fullfile(output_dir, 'refined_local_stability_summary.csv'));
make_refined_figure(output_dir, records, distance_histories, cfg);

results = struct();
results.input_file = input_file;
results.phase_step = phase_step;
results.template = template;
results.records = records;
results.distance_histories = distance_histories;
save(fullfile(output_dir, 'refined_local_stability.mat'), 'results', '-v7.3');

analysis_cfg = cfg;
analysis_cfg.experiment_name = 'periodic_orbit_local_stability_reanalysis_v1';
analysis_cfg.analysis = struct('input_file', input_file, ...
    'template_phase_step', phase_step, ...
    'symmetries', 'identity, reflection');
metadata = brusselator_run_metadata(analysis_cfg, mfilename, toc(timer));
brusselator_write_metadata(output_dir, analysis_cfg, metadata);
disp(summary_table);
end

function dense = densely_interpolate_template(template, period, phase_step)
query = (0:phase_step:period)';
base_phase = template.phase_times(:);
dense = struct();
dense.phase_times = query;
dense.U = interp1(base_phase, template.U, query, 'linear');
dense.V = interp1(base_phase, template.V, query, 'linear');
end

function indices = nearest_time_indices(S, sample_times)
indices = zeros(numel(sample_times), 1);
for j = 1:numel(sample_times)
    [~, indices(j)] = min(abs(S - sample_times(j)));
end
end

function [distance, lag, reflection] = dense_template_distance(u, v, template)
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
    lag = template.phase_times(reflection_index);
    reflection = true;
else
    distance = direct_minimum;
    lag = template.phase_times(direct_index);
    reflection = false;
end
end

function label = classify_record(record, classification)
if record.median_last_10_orbit_distance < classification.recovery_threshold && ...
        record.maximum_last_10_orbit_distance < classification.transition_threshold
    label = 'returns_to_candidate_orbit_neighborhood';
elseif record.final_orbit_distance > classification.transition_threshold
    label = 'does_not_return_within_test_window';
else
    label = 'intermediate_or_inconclusive';
end
end

function make_refined_figure(outdir, records, histories, cfg)
fig = figure('Color', 'w', 'Position', [100 100 1000 600]);
hold on;
for k = 1:numel(records)
    semilogy(histories{k}.sample_times, max(histories{k}.orbit_distance, eps), ...
        'LineWidth', 1.0, 'DisplayName', sprintf('p=%0.2f, %s, %0.0e', ...
        records(k).phase_fraction, records(k).direction, ...
        records(k).amplitude_relative_full_state));
end
yline(cfg.classification.recovery_threshold, '--', 'recovery-scale threshold', ...
    'HandleVisibility', 'off');
xlabel('Frozen b=10 integration time');
ylabel('Interpolated orbit distance modulo phase/reflection');
title('Refined local stability metric for the candidate periodic orbit');
legend('Location', 'eastoutside', 'FontSize', 7); grid on; box on;
exportgraphics(fig, fullfile(outdir, 'refined_orbital_distance_by_case.png'), ...
    'Resolution', 300);
close(fig);
end
