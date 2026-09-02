function results = run_phase_aware_basin_slice()
%RUN_PHASE_AWARE_BASIN_SLICE Test a full-state slice between identified outcomes.
%   This experiment replaces the earlier final-v-component slice diagnostic.
%   Its endpoints are (A) the stationary frozen b=10 baseline state and (B) a
%   Poincare-standardized state on the candidate reflection-symmetric periodic
%   orbit. Each full-state interpolation X_lambda = (1-lambda)A + lambda B is
%   evolved with frozen b=10 dynamics and classified against *both* endpoints:
%   the stationary state modulo reflection, and the periodic-orbit template
%   modulo temporal phase and reflection.
%
%   A sharp finite-time transition in this one slice would be evidence
%   consistent with a separatrix between the two sampled outcomes. It is not a
%   computation of a basin boundary. Unresolved trajectories are retained and
%   receive a longer edge-candidate extension instead of being forced into an
%   outcome class.

study_id = 'phase_aware_basin_slice_v1';
parent_file = fullfile('experiment_outputs', ...
    'frozen_full_state_confirmation_v1', 'frozen_full_state_trajectories.mat');
outdir = fullfile('experiment_outputs', study_id);
if ~isfile(parent_file)
    error('run_phase_aware_basin_slice:MissingParentData', ...
        'Expected frozen full-state data at %s.', parent_file);
end
if isfolder(outdir)
    error('run_phase_aware_basin_slice:OutputExists', ...
        'Refusing to overwrite existing experiment directory: %s', outdir);
end

timer = tic;
loaded = load(parent_file, 'results');
parent = loaded.results;
if ~isfinite(parent.forced_section.period) || parent.forced_section.period <= 0
    error('run_phase_aware_basin_slice:MissingPeriod', ...
        'The parent result has no valid forced Poincare crossing period.');
end

solver_settings = struct('name', 'ode15s', 'output_dt', 0.1, ...
    'RelTol', 1e-8, 'AbsTol', 1e-10, 'MaxStep', 0.1);
cfg = brusselator_default_config(struct( ...
    'experiment_name', study_id, ...
    'model', parent.configuration.model, ...
    'grid', parent.configuration.grid, ...
    'forcing', parent.configuration.forcing, ...
    'initial', parent.configuration.initial, ...
    'solver', solver_settings));
cfg.basin_slice = struct( ...
    'parent_file', parent_file, ...
    'frozen_b', parent.configuration.confirmation.frozen_b, ...
    'source_bmax', parent.configuration.confirmation.source_bmax, ...
    'endpoint_A', 'late stationary control state from frozen confirmation', ...
    'endpoint_B', 'Poincare-standardized candidate periodic-orbit state', ...
    'observed_section_period', parent.forced_section.period, ...
    'observed_direct_period', 2 * parent.forced_section.period, ...
    'template_phase_step', 0.005, ...
    'initial_lambdas', [0, 0.65, 0.70, 0.725, 0.75, 0.80, 1], ...
    'coarse_duration', 80, ...
    'bisection_duration', 80, ...
    'edge_extension_duration', 240, ...
    'maximum_bisection_steps', 4, ...
    'distance_sample_interval', 1, ...
    'late_window', 20, ...
    'late_max_threshold', 2e-2, ...
    'metric', ['relative_l2_full_state; A minimized over spatial reflection; ', ...
        'B minimized over phase and spatial reflection']);

mkdir(outdir);
brusselator_write_metadata(outdir, cfg, ...
    brusselator_run_metadata(cfg, mfilename, NaN));

[A, B, template, endpoint_info] = construct_endpoints_and_template(parent, cfg);
cfg.basin_slice.endpoint_B_section_index = endpoint_info.section_index;
cfg.basin_slice.endpoint_B_section_time = endpoint_info.section_time;
cfg.basin_slice.endpoint_B_template_end_time = endpoint_info.template_end_time;
brusselator_write_metadata(outdir, cfg, ...
    brusselator_run_metadata(cfg, mfilename, NaN));

x = linspace(0, cfg.grid.Lx, cfg.grid.N)';
frozen_par = brusselator_make_parameters(cfg, @(t) cfg.basin_slice.frozen_b);
record_template = empty_record_template();
records = record_template([]);
trajectory_data = cell(0, 1);

fprintf('Phase-aware basin slice: %d coarse states at frozen b = %.2f.\n', ...
    numel(cfg.basin_slice.initial_lambdas), cfg.basin_slice.frozen_b);
for k = 1:numel(cfg.basin_slice.initial_lambdas)
    [record, data] = execute_slice_case(cfg.basin_slice.initial_lambdas(k), ...
        cfg.basin_slice.coarse_duration, 'coarse', A, B, template, cfg, ...
        frozen_par, x, record_template);
    records(end+1, 1) = record; %#ok<AGROW>
    trajectory_data{end+1, 1} = data; %#ok<AGROW>
    write_checkpoint(outdir, cfg, A, B, template, endpoint_info, records, trajectory_data);
end

[lower_lambda, upper_lambda, bracket_status] = locate_outcome_bracket(records);
bisection_status = bracket_status;
for step = 1:cfg.basin_slice.maximum_bisection_steps
    if ~strcmp(bisection_status, 'bracket_found')
        break;
    end
    midpoint = (lower_lambda + upper_lambda) / 2;
    fprintf('Bisection %d/%d at lambda = %.8f in [%.8f, %.8f].\n', ...
        step, cfg.basin_slice.maximum_bisection_steps, midpoint, ...
        lower_lambda, upper_lambda);
    [record, data] = execute_slice_case(midpoint, ...
        cfg.basin_slice.bisection_duration, 'bisection', A, B, template, ...
        cfg, frozen_par, x, record_template);
    records(end+1, 1) = record; %#ok<AGROW>
    trajectory_data{end+1, 1} = data; %#ok<AGROW>
    write_checkpoint(outdir, cfg, A, B, template, endpoint_info, records, trajectory_data);

    if strcmp(record.classification, 'returns_to_baseline_neighborhood')
        lower_lambda = midpoint;
    elseif strcmp(record.classification, 'returns_to_periodic_orbit_neighborhood')
        upper_lambda = midpoint;
    else
        fprintf(['Midpoint is unresolved at %.0f time units; extending to %.0f ', ...
            'time units as an edge-candidate check.\n'], ...
            cfg.basin_slice.bisection_duration, cfg.basin_slice.edge_extension_duration);
        [extension_record, extension_data] = execute_slice_case(midpoint, ...
            cfg.basin_slice.edge_extension_duration, 'edge_extension', A, B, ...
            template, cfg, frozen_par, x, record_template);
        records(end+1, 1) = extension_record; %#ok<AGROW>
        trajectory_data{end+1, 1} = extension_data; %#ok<AGROW>
        write_checkpoint(outdir, cfg, A, B, template, endpoint_info, records, trajectory_data);
        if strcmp(extension_record.classification, 'returns_to_baseline_neighborhood')
            lower_lambda = midpoint;
        elseif strcmp(extension_record.classification, ...
                'returns_to_periodic_orbit_neighborhood')
            upper_lambda = midpoint;
        else
            bisection_status = 'stopped_at_persistent_unresolved_edge_candidate';
        end
    end
end

summary = summarize_slice(records, lower_lambda, upper_lambda, bisection_status, cfg);
summary_table = struct2table(records, 'AsArray', true);
writetable(summary_table, fullfile(outdir, 'phase_aware_basin_slice_summary.csv'));
writetable(struct2table(summary, 'AsArray', true), ...
    fullfile(outdir, 'phase_aware_basin_slice_bracket.csv'));
make_slice_figures(outdir, records, trajectory_data, cfg);

results = struct();
results.configuration = cfg;
results.endpoint_A = A;
results.endpoint_B = B;
results.endpoint_info = endpoint_info;
results.template = template;
results.records = records;
results.trajectory_data = trajectory_data;
results.summary = summary;
save(fullfile(outdir, 'raw_phase_aware_basin_slice.mat'), 'results', '-v7.3');

metadata = brusselator_run_metadata(cfg, mfilename, toc(timer));
brusselator_write_metadata(outdir, cfg, metadata);
disp(summary_table);
disp(struct2table(summary, 'AsArray', true));
fprintf('Phase-aware basin slice complete. Results saved in: %s\n', outdir);
end

function record = empty_record_template()
record = struct( ...
    'lambda', NaN, ...
    'stage', '', ...
    'duration', NaN, ...
    'initial_distance_to_A', NaN, ...
    'initial_distance_to_orbit', NaN, ...
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

function [A, B, template, info] = construct_endpoints_and_template(parent, cfg)
S = parent.trajectories.S;
U_control = parent.trajectories.U_control;
V_control = parent.trajectories.V_control;
U_forced = parent.trajectories.U_forced;
V_forced = parent.trajectories.V_forced;

A = struct('U', U_control(end,:), 'V', V_control(end,:));
direct_period = cfg.basin_slice.observed_direct_period;
section_times = parent.forced_section.times(:);
valid = find(section_times >= 0.55 * S(end) & ...
    section_times + direct_period <= S(end));
if isempty(valid)
    valid = find(section_times + direct_period <= S(end));
end
if isempty(valid)
    error('run_phase_aware_basin_slice:InsufficientTemplateWindow', ...
        'The frozen forced trajectory is too short for one direct period.');
end
section_index = valid(1);
section_time = section_times(section_index);
state = parent.forced_section.states(section_index,:);
N = cfg.grid.N;
if numel(state) ~= 2 * N
    error('run_phase_aware_basin_slice:StateDimensionMismatch', ...
        'Poincare state dimension does not match the configured grid.');
end
B = struct('U', state(1:N), 'V', state(N+1:end));

phase_times = (0:cfg.basin_slice.template_phase_step:direct_period)';
if phase_times(end) < direct_period
    phase_times(end+1, 1) = direct_period;
end
query_times = section_time + phase_times;
template = struct();
template.phase_times = phase_times;
template.U = interp1(S, U_forced, query_times, 'linear');
template.V = interp1(S, V_forced, query_times, 'linear');
% The Poincare-interpolated endpoint is the exact template phase zero.
template.U(1,:) = B.U;
template.V(1,:) = B.V;

info = struct('section_index', section_index, 'section_time', section_time, ...
    'template_end_time', query_times(end), ...
    'observed_direct_period', direct_period);
end

function [record, data] = execute_slice_case(lambda, duration, stage, A, B, ...
        template, cfg, frozen_par, x, record_template)
u0 = (1 - lambda) * A.U + lambda * B.U;
v0 = (1 - lambda) * A.V + lambda * B.V;
initial_state = [x, u0', v0'];
fprintf('  %s: lambda = %.8f, duration = %.0f ...\n', stage, lambda, duration);
case_timer = tic;
[final_state, S, V, U] = solve_brusselator_1d_forced( ...
    initial_state, frozen_par, duration, 0);
runtime = toc(case_timer);

stride = max(1, round(cfg.basin_slice.distance_sample_interval / ...
    cfg.solver.output_dt));
sample_indices = unique([1:stride:numel(S), numel(S)]);
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
sample_times = S(sample_indices);
late = sample_times >= duration - cfg.basin_slice.late_window;
if ~any(late)
    error('run_phase_aware_basin_slice:LateWindowEmpty', ...
        'The requested duration is shorter than the late comparison window.');
end

record = record_template;
record.lambda = lambda;
record.stage = stage;
record.duration = duration;
record.initial_distance_to_A = dist_A(1);
record.initial_distance_to_orbit = dist_orbit(1);
record.final_distance_to_A = dist_A(end);
record.final_distance_to_orbit = dist_orbit(end);
record.late_median_distance_to_A = median(dist_A(late));
record.late_maximum_distance_to_A = max(dist_A(late));
record.late_median_distance_to_orbit = median(dist_orbit(late));
record.late_maximum_distance_to_orbit = max(dist_orbit(late));
record.final_A_uses_reflection = A_reflection(end);
record.final_orbit_uses_reflection = orbit_reflection(end);
record.final_orbit_phase_lag = orbit_lag(end);
record.runtime_seconds = runtime;
record.classification = classify_slice_record(record, cfg);

fprintf('    %s; late medians: d_A = %.3e, d_orbit = %.3e (%.1f s).\n', ...
    record.classification, record.late_median_distance_to_A, ...
    record.late_median_distance_to_orbit, runtime);
data = struct('S', S, 'U', U, 'V', V, 'final_state', final_state, ...
    'sample_times', sample_times, 'distance_to_A', dist_A, ...
    'distance_to_orbit', dist_orbit, 'A_uses_reflection', A_reflection, ...
    'orbit_phase_lag', orbit_lag, 'orbit_uses_reflection', orbit_reflection);
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

function label = classify_slice_record(record, cfg)
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
    label = 'unresolved_within_test_window';
end
end

function [lower_lambda, upper_lambda, status] = locate_outcome_bracket(records)
lower_lambda = NaN;
upper_lambda = NaN;
status = 'no_ordered_baseline_to_periodic_bracket_in_coarse_scan';
coarse = records(strcmp({records.stage}, 'coarse'));
if isempty(coarse)
    return;
end
[~, ordering] = sort([coarse.lambda]);
coarse = coarse(ordering);
for k = 1:(numel(coarse) - 1)
    if strcmp(coarse(k).classification, 'returns_to_baseline_neighborhood') && ...
            strcmp(coarse(k+1).classification, ...
            'returns_to_periodic_orbit_neighborhood')
        lower_lambda = coarse(k).lambda;
        upper_lambda = coarse(k+1).lambda;
        status = 'bracket_found';
        return;
    end
end
end

function summary = summarize_slice(records, lower_lambda, upper_lambda, status, cfg)
summary = struct();
summary.number_of_trajectories = numel(records);
summary.number_baseline = sum(strcmp({records.classification}, ...
    'returns_to_baseline_neighborhood'));
summary.number_periodic_orbit = sum(strcmp({records.classification}, ...
    'returns_to_periodic_orbit_neighborhood'));
summary.number_unresolved = sum(strcmp({records.classification}, ...
    'unresolved_within_test_window'));
summary.bisection_status = status;
summary.lower_baseline_lambda = lower_lambda;
summary.upper_periodic_lambda = upper_lambda;
if strcmp(status, 'bracket_found')
    summary.finite_time_transition_interval_width = upper_lambda - lower_lambda;
    summary.finite_time_transition_midpoint = (lower_lambda + upper_lambda) / 2;
else
    summary.finite_time_transition_interval_width = NaN;
    summary.finite_time_transition_midpoint = NaN;
end
summary.frozen_b = cfg.basin_slice.frozen_b;
summary.coarse_duration = cfg.basin_slice.coarse_duration;
summary.bisection_duration = cfg.basin_slice.bisection_duration;
summary.classification_threshold = cfg.classification.recovery_threshold;
summary.interpretation = ['finite-time full-state slice classification; ', ...
    'not a full basin-boundary computation'];
end

function write_checkpoint(outdir, cfg, A, B, template, endpoint_info, records, trajectory_data)
progress_table = struct2table(records, 'AsArray', true);
writetable(progress_table, fullfile(outdir, 'progress_summary.csv'));
save(fullfile(outdir, 'progress_phase_aware_basin_slice.mat'), 'cfg', 'A', ...
    'B', 'template', 'endpoint_info', 'records', 'trajectory_data', '-v7.3');
end

function make_slice_figures(outdir, records, trajectories, cfg)
fig = figure('Color', 'w', 'Position', [100 100 950 500]);
lambdas = [records.lambda];
semilogy(lambdas, max([records.late_median_distance_to_A], eps), 'o-', ...
    'LineWidth', 1.2, 'DisplayName', 'late median distance to A'); hold on;
semilogy(lambdas, max([records.late_median_distance_to_orbit], eps), 's-', ...
    'LineWidth', 1.2, 'DisplayName', 'late median distance to periodic orbit');
yline(cfg.classification.recovery_threshold, '--', 'classification threshold', ...
    'HandleVisibility', 'off');
xlabel('\lambda in X_\lambda = (1-\lambda)A + \lambda B');
ylabel('Relative full-state distance');
title('Phase- and reflection-aware frozen b=10 slice outcomes');
legend('Location', 'best'); grid on; box on;
exportgraphics(fig, fullfile(outdir, 'late_distance_by_lambda.png'), ...
    'Resolution', 300);
close(fig);

fig = figure('Color', 'w', 'Position', [100 100 1050 650]);
tiledlayout(2, 1, 'TileSpacing', 'compact');
nexttile; hold on;
for k = 1:numel(records)
    data = trajectories{k};
    semilogy(data.sample_times, max(data.distance_to_A, eps), 'LineWidth', 0.9, ...
        'DisplayName', sprintf('\\lambda=%.5f (%s)', records(k).lambda, ...
        records(k).classification));
end
yline(cfg.classification.recovery_threshold, '--', 'HandleVisibility', 'off');
ylabel('d_A'); grid on; box on;
title('Frozen-slice trajectories: distance to stationary baseline');
legend('Location', 'eastoutside', 'FontSize', 6);
nexttile; hold on;
for k = 1:numel(records)
    data = trajectories{k};
    semilogy(data.sample_times, max(data.distance_to_orbit, eps), 'LineWidth', 0.9, ...
        'HandleVisibility', 'off');
end
yline(cfg.classification.recovery_threshold, '--', 'HandleVisibility', 'off');
xlabel('Frozen b=10 integration time'); ylabel('d_{orbit}'); grid on; box on;
title('Distance to candidate periodic orbit modulo phase/reflection');
exportgraphics(fig, fullfile(outdir, 'slice_distance_histories.png'), ...
    'Resolution', 300);
close(fig);
end
