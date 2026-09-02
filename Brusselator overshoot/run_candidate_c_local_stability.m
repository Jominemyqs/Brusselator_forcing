function results = run_candidate_c_local_stability()
%RUN_CANDIDATE_C_LOCAL_STABILITY Probe attraction near candidate periodic C.
%   Takes the freshly confirmed lambda=0.8 periodic-pattern candidate and
%   applies smooth Neumann-compatible perturbations at standardized phases.
%   Every trajectory is compared with the full candidate-C orbit modulo phase
%   and spatial reflection. This is a finite-direction, finite-time stability
%   screen; it is not a Floquet calculation or proof of an attracting basin.

study_id = 'candidate_c_local_stability_v1';
parent_file = fullfile('experiment_outputs', ...
    'edge_candidate_frozen_confirmation_v1', ...
    'candidate_confirmation_trajectory.mat');
outdir = fullfile('experiment_outputs', study_id);
if ~isfile(parent_file)
    error('run_candidate_c_local_stability:MissingParentData', ...
        'Expected candidate-C confirmation data at %s.', parent_file);
end
if isfolder(outdir)
    error('run_candidate_c_local_stability:OutputExists', ...
        'Refusing to overwrite existing experiment directory: %s', outdir);
end

loaded = load(parent_file, 'results');
parent = loaded.results;
if ~isfinite(parent.section.period) || parent.section.period <= 0
    error('run_candidate_c_local_stability:MissingPeriod', ...
        'The candidate-C confirmation has no valid Poincare period.');
end

cfg = parent.configuration;
cfg.experiment_name = study_id;
cfg.candidate_c_local_stability = struct( ...
    'parent_file', parent_file, ...
    'frozen_b', parent.configuration.candidate_confirmation.frozen_b, ...
    'lambda_source', parent.configuration.candidate_confirmation.lambda, ...
    'observed_direct_period', parent.section.period, ...
    'duration', 40, ...
    'late_window', 10, ...
    'template_phase_step', 0.005, ...
    'phase_fractions', [0, 0.25, 0.5], ...
    'amplitudes_relative_full_state', [1e-4, 1e-3], ...
    'directions', {{'reflection_odd_smooth', 'reflection_even_smooth'}}, ...
    'orbit_metric', 'minimum_relative_l2_over_phase_and_reflection', ...
    'interpretation', ['finite-time local-return screen around candidate C; ', ...
        'not a Floquet or basin calculation']);

mkdir(outdir);
brusselator_write_metadata(outdir, cfg, ...
    brusselator_run_metadata(cfg, mfilename, NaN));

S_parent = parent.trajectory.S;
U_parent = parent.trajectory.U;
V_parent = parent.trajectory.V;
[anchor_index, anchor_time] = select_anchor(parent.section.times, S_parent, ...
    cfg.candidate_c_local_stability.duration, ...
    cfg.candidate_c_local_stability.observed_direct_period);
template = build_dense_template(S_parent, U_parent, V_parent, ...
    parent.section.states(anchor_index,:), anchor_time, ...
    cfg.candidate_c_local_stability.observed_direct_period, ...
    cfg.candidate_c_local_stability.template_phase_step);
directions = smooth_neumann_directions(cfg.grid.Lx, cfg.grid.N);
x = linspace(0, cfg.grid.Lx, cfg.grid.N)';
frozen_par = brusselator_make_parameters(cfg, ...
    @(t) cfg.candidate_c_local_stability.frozen_b);

phase_fractions = cfg.candidate_c_local_stability.phase_fractions;
amplitudes = cfg.candidate_c_local_stability.amplitudes_relative_full_state;
record_template = empty_record_template();
case_count = numel(phase_fractions) * numel(directions) * numel(amplitudes);
records = repmat(record_template, case_count, 1);
trajectory_data = cell(case_count, 1);
case_number = 0;
for p = 1:numel(phase_fractions)
    [base_u, base_v] = template_state_at_phase(template, ...
        phase_fractions(p) * cfg.candidate_c_local_stability.observed_direct_period);
    for d = 1:numel(directions)
        for a = 1:numel(amplitudes)
            case_number = case_number + 1;
            [u0, v0] = apply_relative_perturbation(base_u, base_v, ...
                directions(d), amplitudes(a));
            initial_state = [x, u0', v0'];
            fprintf('Case %d/%d: phase %.2f, %s, relative amplitude %.1e ...\n', ...
                case_number, case_count, phase_fractions(p), ...
                directions(d).name, amplitudes(a));
            case_timer = tic;
            [final_state, S, V, U] = solve_brusselator_1d_forced( ...
                initial_state, frozen_par, cfg.candidate_c_local_stability.duration, 0);
            runtime = toc(case_timer);
            [sample_times, distances, lags, reflections] = sample_orbit_distances( ...
                S, U, V, template, cfg);
            late = sample_times >= cfg.candidate_c_local_stability.duration - ...
                cfg.candidate_c_local_stability.late_window;

            record = record_template;
            record.phase_fraction = phase_fractions(p);
            record.phase_start_time = anchor_time + ...
                phase_fractions(p) * cfg.candidate_c_local_stability.observed_direct_period;
            record.direction = directions(d).name;
            record.amplitude_relative_full_state = amplitudes(a);
            record.initial_orbit_distance = distances(1);
            record.final_orbit_distance = distances(end);
            record.median_last_10_orbit_distance = median(distances(late));
            record.maximum_last_10_orbit_distance = max(distances(late));
            record.best_final_lag = lags(end);
            record.best_final_reflection = reflections(end);
            record.runtime_seconds = runtime;
            record.classification = classify_record(record, cfg);
            records(case_number) = record;
            trajectory_data{case_number} = struct('S', S, 'U', U, 'V', V, ...
                'final_state', final_state, 'sample_times', sample_times, ...
                'orbit_distance', distances, 'best_lag', lags, ...
                'best_reflection', reflections);

            write_checkpoint(outdir, cfg, records(1:case_number), ...
                trajectory_data(1:case_number), template, directions, anchor_time);
            fprintf('  %s; late median orbit distance %.3e (%.1f s).\n', ...
                record.classification, record.median_last_10_orbit_distance, runtime);
        end
    end
end

summary_table = struct2table(records, 'AsArray', true);
writetable(summary_table, fullfile(outdir, 'candidate_c_local_stability_summary.csv'));
make_figures(outdir, records, trajectory_data, cfg);

results = struct();
results.configuration = cfg;
results.template = template;
results.directions = directions;
results.anchor_time = anchor_time;
results.records = records;
results.trajectory_data = trajectory_data;
save(fullfile(outdir, 'raw_candidate_c_local_stability.mat'), 'results', '-v7.3');

metadata = brusselator_run_metadata(cfg, mfilename, sum([records.runtime_seconds]));
brusselator_write_metadata(outdir, cfg, metadata);
disp(summary_table);
fprintf('Candidate-C local stability screen complete. Results saved in: %s\n', outdir);
end

function [index, time] = select_anchor(crossing_times, S, duration, period)
candidates = find(crossing_times >= 0.55 * S(end) & ...
    crossing_times + period + duration <= S(end));
if isempty(candidates)
    candidates = find(crossing_times + period + duration <= S(end));
end
if isempty(candidates)
    error('run_candidate_c_local_stability:InsufficientReferenceWindow', ...
        'The parent confirmation is too short for the requested stability duration.');
end
index = candidates(1);
time = crossing_times(index);
end

function template = build_dense_template(S, U, V, section_state, start_time, period, phase_step)
phase_times = (0:phase_step:period)';
if phase_times(end) < period
    phase_times(end+1, 1) = period;
end
query_times = start_time + phase_times;
if query_times(end) > S(end)
    error('run_candidate_c_local_stability:TemplateOutOfRange', ...
        'The requested template extends beyond the saved confirmation data.');
end
template = struct('phase_times', phase_times, ...
    'U', interp1(S, U, query_times, 'linear'), ...
    'V', interp1(S, V, query_times, 'linear'));
N = size(U,2);
template.U(1,:) = section_state(1:N);
template.V(1,:) = section_state(N+1:end);
end

function [u, v] = template_state_at_phase(template, phase)
u = interp1(template.phase_times, template.U, phase, 'linear');
v = interp1(template.phase_times, template.V, phase, 'linear');
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

function [sample_times, distances, lags, reflections] = sample_orbit_distances(S, U, V, template, cfg)
stride = max(1, round(1 / cfg.solver.output_dt));
indices = unique([1:stride:numel(S), numel(S)]);
sample_times = S(indices);
distances = zeros(numel(indices), 1);
lags = zeros(numel(indices), 1);
reflections = false(numel(indices), 1);
for j = 1:numel(indices)
    [distances(j), lags(j), reflections(j)] = distance_to_template( ...
        U(indices(j),:), V(indices(j),:), template);
end
end

function [distance, phase_lag, uses_reflection] = distance_to_template(u, v, template)
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

function record = empty_record_template()
record = struct('phase_fraction', NaN, 'phase_start_time', NaN, ...
    'direction', '', 'amplitude_relative_full_state', NaN, ...
    'initial_orbit_distance', NaN, 'final_orbit_distance', NaN, ...
    'median_last_10_orbit_distance', NaN, ...
    'maximum_last_10_orbit_distance', NaN, 'best_final_lag', NaN, ...
    'best_final_reflection', false, 'runtime_seconds', NaN, ...
    'classification', 'unclassified');
end

function label = classify_record(record, cfg)
if record.median_last_10_orbit_distance < cfg.classification.recovery_threshold && ...
        record.maximum_last_10_orbit_distance < cfg.basin_slice.late_max_threshold
    label = 'returns_to_candidate_C_orbit_neighborhood';
elseif record.final_orbit_distance > cfg.classification.transition_threshold
    label = 'does_not_return_within_test_window';
else
    label = 'intermediate_or_inconclusive';
end
end

function write_checkpoint(outdir, cfg, records, trajectory_data, template, directions, anchor_time)
writetable(struct2table(records, 'AsArray', true), ...
    fullfile(outdir, 'progress_summary.csv'));
save(fullfile(outdir, 'progress_candidate_c_local_stability.mat'), 'cfg', ...
    'records', 'trajectory_data', 'template', 'directions', 'anchor_time', '-v7.3');
end

function make_figures(outdir, records, trajectories, cfg)
fig = figure('Color', 'w', 'Position', [100 100 1050 650]);
hold on;
for k = 1:numel(records)
    data = trajectories{k};
    semilogy(data.sample_times, max(data.orbit_distance, eps), 'LineWidth', 1.0, ...
        'DisplayName', sprintf('p=%.2f, %s, %.0e', ...
        records(k).phase_fraction, records(k).direction, ...
        records(k).amplitude_relative_full_state));
end
yline(cfg.classification.recovery_threshold, '--', 'recovery threshold', ...
    'HandleVisibility', 'off');
xlabel('Frozen b=10 integration time');
ylabel('Distance to candidate-C orbit modulo phase/reflection');
title('Finite-time local perturbation screen for candidate C');
legend('Location', 'eastoutside', 'FontSize', 7); grid on; box on;
exportgraphics(fig, fullfile(outdir, 'candidate_c_orbit_distance_by_case.png'), ...
    'Resolution', 300);
close(fig);
end
