function results = run_extended_even_sector_stability()
%RUN_EXTENDED_EVEN_SECTOR_STABILITY Resolve the slow reflection-even response.
%   Extends the potentially slow reflection-even local-stability direction for
%   roughly 40 direct periods. A zero-perturbation restart establishes the
%   numerical orbit-distance floor, while three smooth perturbation amplitudes
%   test whether the observed short-run growth later decays or persists.

study_id = 'extended_even_sector_stability_v1';
parent_dir = fullfile('experiment_outputs', 'periodic_orbit_local_stability_v1');
parent_file = fullfile(parent_dir, 'raw_local_stability_trajectories.mat');
outdir = fullfile('experiment_outputs', study_id);
if ~isfile(parent_file)
    error('run_extended_even_sector_stability:MissingParentData', ...
        'Expected local-stability data at %s.', parent_file);
end
if isfolder(outdir)
    error('run_extended_even_sector_stability:OutputExists', ...
        'Refusing to overwrite existing experiment directory: %s', outdir);
end

loaded = load(parent_file, 'results');
parent = loaded.results;
cfg = parent.configuration;
duration = 160;
amplitudes = [0, 1e-5, 1e-4, 1e-3];
direction_index = find(strcmp({parent.directions.name}, 'reflection_even_smooth'), 1);
if isempty(direction_index)
    error('run_extended_even_sector_stability:MissingDirection', ...
        'The parent test has no reflection-even direction.');
end
direction = parent.directions(direction_index);
cfg.experiment_name = study_id;
cfg.local_stability_extension = struct( ...
    'parent_file', parent_file, ...
    'duration', duration, ...
    'direct_period', cfg.local_stability.observed_direct_period, ...
    'amplitudes_relative_full_state', amplitudes, ...
    'direction', direction.name, ...
    'orbit_template_phase_step', 0.0025);

mkdir(outdir);
brusselator_write_metadata(outdir, cfg, ...
    brusselator_run_metadata(cfg, mfilename, NaN));

template = densely_interpolate_template(parent.template, ...
    cfg.local_stability.observed_direct_period, ...
    cfg.local_stability_extension.orbit_template_phase_step);
x = linspace(0, cfg.grid.Lx, cfg.grid.N)';
base_u = parent.template.U(1,:);
base_v = parent.template.V(1,:);
frozen_par = brusselator_make_parameters(cfg, @(t) cfg.local_stability.frozen_b);
sample_stride = round(2 / cfg.solver.output_dt);

record_template = struct( ...
    'amplitude_relative_full_state', NaN, 'initial_orbit_distance', NaN, ...
    'final_orbit_distance', NaN, 'median_last_20_orbit_distance', NaN, ...
    'maximum_last_20_orbit_distance', NaN, 'late_to_initial_ratio', NaN, ...
    'runtime_seconds', NaN, 'classification', 'unclassified');
records = repmat(record_template, numel(amplitudes), 1);
trajectory_data = cell(numel(amplitudes), 1);

for k = 1:numel(amplitudes)
    [u0, v0] = apply_relative_perturbation(base_u, base_v, direction, amplitudes(k));
    fprintf('Even-sector amplitude %.1e for %.0f frozen time units ...\n', ...
        amplitudes(k), duration);
    case_timer = tic;
    [final_state, S, V, U] = solve_brusselator_1d_forced( ...
        [x, u0', v0'], frozen_par, duration, 0);
    runtime = toc(case_timer);
    sample_indices = unique([1:sample_stride:numel(S), numel(S)]);
    distances = zeros(numel(sample_indices), 1);
    lags = zeros(numel(sample_indices), 1);
    reflections = false(numel(sample_indices), 1);
    for j = 1:numel(sample_indices)
        [distances(j), lags(j), reflections(j)] = dense_template_distance( ...
            U(sample_indices(j),:), V(sample_indices(j),:), template);
    end
    late = S(sample_indices) >= duration - 20;

    record = record_template;
    record.amplitude_relative_full_state = amplitudes(k);
    record.initial_orbit_distance = distances(1);
    record.final_orbit_distance = distances(end);
    record.median_last_20_orbit_distance = median(distances(late));
    record.maximum_last_20_orbit_distance = max(distances(late));
    if amplitudes(k) > 0
        record.late_to_initial_ratio = record.median_last_20_orbit_distance / amplitudes(k);
        if record.median_last_20_orbit_distance < cfg.classification.recovery_threshold
            record.classification = 'remains_in_candidate_orbit_neighborhood';
        else
            record.classification = 'leaves_candidate_orbit_neighborhood';
        end
    else
        record.classification = 'unperturbed_restart_reference';
    end
    record.runtime_seconds = runtime;
    records(k) = record;
    trajectory_data{k} = struct('S', S, 'U', U, 'V', V, ...
        'final_state', final_state, 'sample_times', S(sample_indices), ...
        'orbit_distance', distances, 'best_lag', lags, ...
        'best_reflection', reflections);

    writetable(struct2table(records(1:k), 'AsArray', true), ...
        fullfile(outdir, 'progress_summary.csv'));
    save(fullfile(outdir, 'progress_even_sector.mat'), 'cfg', 'records', ...
        'trajectory_data', 'template', 'direction', '-v7.3');
end

summary_table = struct2table(records, 'AsArray', true);
writetable(summary_table, fullfile(outdir, 'extended_even_sector_summary.csv'));
make_figure(outdir, records, trajectory_data, cfg);
results = struct('configuration', cfg, 'records', records, 'template', template, ...
    'direction', direction, 'trajectory_data', {trajectory_data});
save(fullfile(outdir, 'raw_extended_even_sector_trajectories.mat'), 'results', '-v7.3');
metadata = brusselator_run_metadata(cfg, mfilename, sum([records.runtime_seconds]));
brusselator_write_metadata(outdir, cfg, metadata);
disp(summary_table);
fprintf('Extended even-sector stability test complete. Results saved in: %s\n', outdir);
end

function template = densely_interpolate_template(raw_template, period, phase_step)
query = (0:phase_step:period)';
template = struct('phase_times', query, ...
    'U', interp1(raw_template.phase_times(:), raw_template.U, query, 'linear'), ...
    'V', interp1(raw_template.phase_times(:), raw_template.V, query, 'linear'));
end

function [u, v] = apply_relative_perturbation(base_u, base_v, direction, amplitude)
base_norm = sqrt(sum(base_u.^2) + sum(base_v.^2));
direction_norm = sqrt(sum(direction.u.^2) + sum(direction.v.^2));
scale = amplitude * base_norm / max(direction_norm, eps);
u = base_u + scale * direction.u;
v = base_v + scale * direction.v;
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

function make_figure(outdir, records, data, cfg)
fig = figure('Color', 'w', 'Position', [100 100 900 500]);
hold on;
for k = 1:numel(records)
    semilogy(data{k}.sample_times, max(data{k}.orbit_distance, eps), ...
        'LineWidth', 1.1, 'DisplayName', sprintf('relative amplitude %.0e', ...
        records(k).amplitude_relative_full_state));
end
yline(cfg.classification.recovery_threshold, '--', 'recovery-scale threshold', ...
    'HandleVisibility', 'off');
xlabel('Frozen b=10 integration time');
ylabel('Interpolated orbit distance modulo phase/reflection');
title('Extended reflection-even local stability test');
legend('Location', 'best'); grid on; box on;
exportgraphics(fig, fullfile(outdir, 'extended_even_sector_orbit_distance.png'), ...
    'Resolution', 300);
close(fig);
end
