function results = run_periodic_orbit_local_stability()
%RUN_PERIODIC_ORBIT_LOCAL_STABILITY Probe attraction near the periodic state.
%   The frozen-b=10 departure trajectory is treated as a candidate periodic
%   orbit with the observed spatiotemporal reflection symmetry. Smooth,
%   Neumann-compatible perturbations are applied at three standardized orbit
%   phases. Each perturbed trajectory is compared against the whole reference
%   orbit modulo temporal phase and x -> L-x reflection.
%
%   This is a local numerical stability test, not a Floquet calculation. A
%   return to the orbit neighborhood supports local attraction; failure to
%   return is reported rather than reclassified silently.

study_id = 'periodic_orbit_local_stability_v1';
parent_dir = fullfile('experiment_outputs', 'frozen_full_state_confirmation_v1');
parent_file = fullfile(parent_dir, 'frozen_full_state_trajectories.mat');
outdir = fullfile('experiment_outputs', study_id);
if ~isfile(parent_file)
    error('run_periodic_orbit_local_stability:MissingParentData', ...
        'Expected full-state frozen output at %s.', parent_file);
end
if isfolder(outdir)
    error('run_periodic_orbit_local_stability:OutputExists', ...
        'Refusing to overwrite existing experiment directory: %s', outdir);
end

loaded = load(parent_file, 'results');
parent = loaded.results;
if ~isfinite(parent.forced_section.period) || parent.forced_section.period <= 0
    error('run_periodic_orbit_local_stability:MissingPeriod', ...
        'The parent result has no valid Poincare crossing period.');
end

duration = 40;
phase_fractions = [0, 0.25, 0.5];
amplitudes = [1e-4, 1e-3];
solver_settings = struct('name', 'ode15s', 'output_dt', 0.05, ...
    'RelTol', 1e-8, 'AbsTol', 1e-10, 'MaxStep', 0.05);
cfg = brusselator_default_config(struct( ...
    'experiment_name', study_id, ...
    'model', parent.configuration.model, ...
    'grid', parent.configuration.grid, ...
    'forcing', parent.configuration.forcing, ...
    'initial', parent.configuration.initial, ...
    'solver', solver_settings));
cfg.local_stability = struct( ...
    'parent_file', parent_file, ...
    'frozen_b', parent.configuration.confirmation.frozen_b, ...
    'source_bmax', parent.configuration.confirmation.source_bmax, ...
    'observed_section_period', parent.forced_section.period, ...
    'observed_direct_period', 2 * parent.forced_section.period, ...
    'duration', duration, ...
    'phase_fractions', phase_fractions, ...
    'amplitudes_relative_full_state', amplitudes, ...
    'directions', {{'reflection_odd_smooth', 'reflection_even_smooth'}}, ...
    'orbit_metric', 'minimum_relative_l2_over_phase_and_reflection');

mkdir(outdir);
brusselator_write_metadata(outdir, cfg, ...
    brusselator_run_metadata(cfg, mfilename, NaN));

S_parent = parent.trajectories.S;
U_parent = parent.trajectories.U_forced;
V_parent = parent.trajectories.V_forced;
full_period = cfg.local_stability.observed_direct_period;
anchor_time = select_anchor_time(parent.forced_section.times, S_parent, duration, full_period);
phase_start_indices = zeros(numel(phase_fractions), 1);
for p = 1:numel(phase_fractions)
    phase_time = anchor_time + phase_fractions(p) * full_period;
    [~, phase_start_indices(p)] = min(abs(S_parent - phase_time));
end

template = build_orbit_template(S_parent, U_parent, V_parent, ...
    phase_start_indices(1), full_period);
directions = smooth_neumann_directions(cfg.grid.Lx, cfg.grid.N);
time_sample_stride = round(2 / cfg.solver.output_dt);
frozen_par = brusselator_make_parameters(cfg, @(t) cfg.local_stability.frozen_b);

record_template = struct( ...
    'phase_fraction', NaN, 'phase_start_time', NaN, 'direction', '', ...
    'amplitude_relative_full_state', NaN, 'initial_orbit_distance', NaN, ...
    'final_orbit_distance', NaN, 'median_last_10_orbit_distance', NaN, ...
    'maximum_last_10_orbit_distance', NaN, 'best_final_lag', NaN, ...
    'best_final_reflection', false, 'runtime_seconds', NaN, ...
    'classification', 'unclassified');
records = repmat(record_template, ...
    numel(phase_fractions) * numel(directions) * numel(amplitudes), 1);
trajectory_data = cell(numel(records), 1);

case_number = 0;
for p = 1:numel(phase_fractions)
    start_index = phase_start_indices(p);
    base_u = U_parent(start_index,:);
    base_v = V_parent(start_index,:);
    x = linspace(0, cfg.grid.Lx, cfg.grid.N)';
    for d = 1:numel(directions)
        for a = 1:numel(amplitudes)
            case_number = case_number + 1;
            [u0, v0] = apply_relative_perturbation(base_u, base_v, ...
                directions(d), amplitudes(a));
            initial_state = [x, u0', v0'];
            fprintf(['Case %d/%d: phase %.2f, %s, relative amplitude %.1e ...\n'], ...
                case_number, numel(records), phase_fractions(p), ...
                directions(d).name, amplitudes(a));
            case_timer = tic;
            [final_state, S, V, U] = solve_brusselator_1d_forced( ...
                initial_state, frozen_par, duration, 0);
            runtime = toc(case_timer);

            sample_indices = unique([1:time_sample_stride:numel(S), numel(S)]);
            orbit_distance = zeros(numel(sample_indices), 1);
            best_lag = zeros(numel(sample_indices), 1);
            best_reflection = false(numel(sample_indices), 1);
            for j = 1:numel(sample_indices)
                [orbit_distance(j), best_lag(j), best_reflection(j)] = ...
                    distance_to_template(U(sample_indices(j),:), ...
                    V(sample_indices(j),:), template);
            end
            last_ten = S(sample_indices) >= duration - 10;

            record = record_template;
            record.phase_fraction = phase_fractions(p);
            record.phase_start_time = S_parent(start_index);
            record.direction = directions(d).name;
            record.amplitude_relative_full_state = amplitudes(a);
            record.initial_orbit_distance = orbit_distance(1);
            record.final_orbit_distance = orbit_distance(end);
            record.median_last_10_orbit_distance = median(orbit_distance(last_ten));
            record.maximum_last_10_orbit_distance = max(orbit_distance(last_ten));
            record.best_final_lag = best_lag(end);
            record.best_final_reflection = best_reflection(end);
            record.runtime_seconds = runtime;
            record.classification = classify_stability_record(record, cfg.classification);
            records(case_number) = record;

            trajectory_data{case_number} = struct( ...
                'S', S, 'U', U, 'V', V, 'final_state', final_state, ...
                'sample_times', S(sample_indices), ...
                'orbit_distance', orbit_distance, 'best_lag', best_lag, ...
                'best_reflection', best_reflection);

            progress_table = struct2table(records(1:case_number), 'AsArray', true);
            writetable(progress_table, fullfile(outdir, 'progress_summary.csv'));
            save(fullfile(outdir, 'progress_local_stability.mat'), 'cfg', ...
                'records', 'trajectory_data', 'template', 'directions', ...
                'phase_start_indices', '-v7.3');
        end
    end
end

summary_table = struct2table(records, 'AsArray', true);
writetable(summary_table, fullfile(outdir, 'local_stability_summary.csv'));
make_local_stability_figures(outdir, records, trajectory_data, cfg);

results = struct();
results.configuration = cfg;
results.records = records;
results.template = template;
results.directions = directions;
results.phase_start_indices = phase_start_indices;
results.trajectory_data = trajectory_data;
save(fullfile(outdir, 'raw_local_stability_trajectories.mat'), 'results', '-v7.3');

metadata = brusselator_run_metadata(cfg, mfilename, sum([records.runtime_seconds]));
brusselator_write_metadata(outdir, cfg, metadata);
disp(summary_table);
fprintf('Local stability test complete. Results saved in: %s\n', outdir);
end

function anchor_time = select_anchor_time(crossing_times, S, duration, full_period)
target = 0.55 * S(end);
candidates = crossing_times(crossing_times >= target & ...
    crossing_times + full_period / 2 + duration <= S(end));
if isempty(candidates)
    candidates = crossing_times(crossing_times + full_period / 2 + duration <= S(end));
end
if isempty(candidates)
    error('run_periodic_orbit_local_stability:InsufficientReferenceWindow', ...
        'The parent trajectory is too short for the requested duration.');
end
anchor_time = candidates(1);
end

function template = build_orbit_template(S, U, V, start_index, full_period)
dt = median(diff(S));
count = round(full_period / dt) + 1;
indices = start_index:(start_index + count - 1);
if indices(end) > numel(S)
    error('run_periodic_orbit_local_stability:TemplateOutOfRange', ...
        'The requested periodic-orbit template extends beyond the parent data.');
end
template = struct('U', U(indices,:), 'V', V(indices,:), ...
    'phase_times', S(indices) - S(start_index));
end

function directions = smooth_neumann_directions(Lx, N)
x = linspace(0, Lx, N);
directions(1).name = 'reflection_odd_smooth';
directions(1).u = cos(pi * x / Lx) + 0.35 * cos(3 * pi * x / Lx);
directions(1).v = -0.7 * cos(pi * x / Lx) + 0.25 * cos(3 * pi * x / Lx);
directions(2).name = 'reflection_even_smooth';
directions(2).u = cos(2 * pi * x / Lx) - 0.4 * cos(4 * pi * x / Lx);
directions(2).v = 0.6 * cos(2 * pi * x / Lx) + 0.2 * cos(4 * pi * x / Lx);
end

function [u, v] = apply_relative_perturbation(base_u, base_v, direction, amplitude)
base_norm = sqrt(sum(base_u.^2) + sum(base_v.^2));
direction_norm = sqrt(sum(direction.u.^2) + sum(direction.v.^2));
scale = amplitude * base_norm / max(direction_norm, eps);
u = base_u + scale * direction.u;
v = base_v + scale * direction.v;
end

function [distance, phase_lag, uses_reflection] = distance_to_template(u, v, template)
U_reflected = fliplr(template.U);
V_reflected = fliplr(template.V);
direct = sqrt(sum((template.U - u).^2 + (template.V - v).^2, 2)) ./ ...
    max(sqrt(sum(template.U.^2 + template.V.^2, 2)), eps);
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

function label = classify_stability_record(record, classification)
if record.median_last_10_orbit_distance < classification.recovery_threshold && ...
        record.maximum_last_10_orbit_distance < classification.transition_threshold
    label = 'returns_to_candidate_orbit_neighborhood';
elseif record.final_orbit_distance > classification.transition_threshold
    label = 'does_not_return_within_test_window';
else
    label = 'intermediate_or_inconclusive';
end
end

function make_local_stability_figures(outdir, records, trajectories, cfg)
fig = figure('Color', 'w', 'Position', [100 100 1000 600]);
hold on;
for k = 1:numel(records)
    data = trajectories{k};
    semilogy(data.sample_times, max(data.orbit_distance, eps), 'LineWidth', 1.0, ...
        'DisplayName', sprintf('p=%0.2f, %s, %0.0e', ...
        records(k).phase_fraction, records(k).direction, ...
        records(k).amplitude_relative_full_state));
end
yline(cfg.classification.recovery_threshold, '--', 'recovery-scale threshold', ...
    'HandleVisibility', 'off');
xlabel('Frozen b=10 integration time');
ylabel('Distance to periodic orbit modulo phase/reflection');
title('Local perturbations of the candidate periodic orbit');
legend('Location', 'eastoutside', 'FontSize', 7); grid on; box on;
exportgraphics(fig, fullfile(outdir, 'orbital_distance_by_case.png'), 'Resolution', 300);
close(fig);

phases = unique([records.phase_fraction]);
amplitudes = unique([records.amplitude_relative_full_state]);
directions = unique({records.direction});
for d = 1:numel(directions)
    values = nan(numel(phases), numel(amplitudes));
    for p = 1:numel(phases)
        for a = 1:numel(amplitudes)
            idx = find([records.phase_fraction] == phases(p) & ...
                strcmp({records.direction}, directions{d}) & ...
                [records.amplitude_relative_full_state] == amplitudes(a), 1);
            values(p,a) = records(idx).median_last_10_orbit_distance;
        end
    end
    fig = figure('Color', 'w', 'Position', [100 100 600 450]);
    imagesc(log10(max(values, eps))); axis xy; colorbar;
    xticks(1:numel(amplitudes)); xticklabels(arrayfun(@(x) sprintf('%.0e', x), ...
        amplitudes, 'UniformOutput', false));
    yticks(1:numel(phases)); yticklabels(arrayfun(@(x) sprintf('%.2f', x), ...
        phases, 'UniformOutput', false));
    xlabel('Relative perturbation amplitude'); ylabel('Orbit phase fraction');
    title(sprintf('log_{10} median late orbital distance: %s', directions{d}));
    exportgraphics(fig, fullfile(outdir, ...
        ['late_distance_', directions{d}, '.png']), 'Resolution', 300);
    close(fig);
end
end
