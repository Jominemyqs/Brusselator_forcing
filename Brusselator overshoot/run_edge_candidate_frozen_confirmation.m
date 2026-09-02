function results = run_edge_candidate_frozen_confirmation()
%RUN_EDGE_CANDIDATE_FROZEN_CONFIRMATION Restart and test candidate state C.
%   Restarts the saved lambda=0.8 late state after the initial 240-unit slice
%   experiment, then evolves it for a fresh 160-unit frozen-b=10 segment. The
%   late half of that segment is tested for full-state Poincare recurrence
%   modulo spatial reflection. This checks persistence beyond the original
%   trajectory, but is not a stability or basin-boundary computation.

study_id = 'edge_candidate_frozen_confirmation_v1';
parent_file = fullfile('experiment_outputs', ...
    'phase_aware_basin_slice_edge_extension_v1', ...
    'raw_edge_extension_trajectories.mat');
outdir = fullfile('experiment_outputs', study_id);
if ~isfile(parent_file)
    error('run_edge_candidate_frozen_confirmation:MissingParentData', ...
        'Expected edge-extension data at %s.', parent_file);
end
if isfolder(outdir)
    error('run_edge_candidate_frozen_confirmation:OutputExists', ...
        'Refusing to overwrite existing experiment directory: %s', outdir);
end

loaded = load(parent_file, 'results');
parent = loaded.results;
candidate_index = find(abs([parent.records.lambda] - 0.8) < 1e-12, 1);
if isempty(candidate_index)
    error('run_edge_candidate_frozen_confirmation:MissingCandidate', ...
        'No lambda = 0.8 state is available in the edge-extension output.');
end
if ~strcmp(parent.records(candidate_index).classification, ...
        'persistent_unresolved_edge_candidate')
    error('run_edge_candidate_frozen_confirmation:UnexpectedClass', ...
        'The lambda = 0.8 state is not the documented unresolved candidate.');
end

cfg = parent.configuration;
cfg.experiment_name = study_id;
cfg.candidate_confirmation = struct( ...
    'parent_file', parent_file, ...
    'lambda', parent.records(candidate_index).lambda, ...
    'source_total_time', parent.records(candidate_index).total_duration, ...
    'duration', 160, ...
    'late_analysis_start', 80, ...
    'frozen_b', cfg.basin_slice.frozen_b, ...
    'metric', 'relative_l2_full_state modulo identity/reflection', ...
    'purpose', ['fresh frozen-system continuation of candidate C; ', ...
        'finite-time recurrence confirmation']);

mkdir(outdir);
brusselator_write_metadata(outdir, cfg, ...
    brusselator_run_metadata(cfg, mfilename, NaN));

initial_state = parent.trajectory_data{candidate_index}.continuation_final_state;
frozen_par = brusselator_make_parameters(cfg, @(t) cfg.candidate_confirmation.frozen_b);
fprintf(['Restarting candidate C from lambda = %.3f, source time %.0f, ', ...
    'for %.0f frozen time units ...\n'], cfg.candidate_confirmation.lambda, ...
    cfg.candidate_confirmation.source_total_time, cfg.candidate_confirmation.duration);
run_timer = tic;
[final_state, S, V, U] = solve_brusselator_1d_forced( ...
    initial_state, frozen_par, cfg.candidate_confirmation.duration, 0);
runtime = toc(run_timer);

late = S >= cfg.candidate_confirmation.late_analysis_start;
section = full_state_poincare_section(S(late), U(late,:), V(late,:));
recurrence = recurrence_scan(section, cfg.classification.recovery_threshold);
classification = classify_confirmation(section, recurrence);

summary = struct();
summary.lambda = cfg.candidate_confirmation.lambda;
summary.source_total_time = cfg.candidate_confirmation.source_total_time;
summary.restart_duration = cfg.candidate_confirmation.duration;
summary.late_analysis_start = cfg.candidate_confirmation.late_analysis_start;
summary.output_dt = median(diff(S));
summary.number_section_crossings = numel(section.times);
summary.section_crossing_period = section.period;
summary.section_crossing_period_cv = section.period_cv;
summary.first_close_return_crossings = recurrence.first_close.return_crossings;
summary.first_close_return_time = recurrence.first_close.return_time;
summary.first_close_return_rel_diff_full_state = ...
    recurrence.first_close.best_median_rel_diff_full_state;
summary.first_close_return_uses_reflection = ...
    recurrence.first_close.best_uses_reflection;
summary.best_return_crossings = recurrence.best.return_crossings;
summary.best_return_rel_diff_full_state = ...
    recurrence.best.best_median_rel_diff_full_state;
summary.best_return_uses_reflection = recurrence.best.best_uses_reflection;
summary.interpretation = classification;

writetable(struct2table(summary, 'AsArray', true), ...
    fullfile(outdir, 'candidate_confirmation_summary.csv'));
writetable(struct2table(recurrence.records, 'AsArray', true), ...
    fullfile(outdir, 'multi_return_recurrence.csv'));
make_figures(outdir, S, V, recurrence, cfg);

results = struct();
results.configuration = cfg;
results.parent_record = parent.records(candidate_index);
results.summary = summary;
results.section = section;
results.recurrence = recurrence;
results.trajectory = struct('S', S, 'U', U, 'V', V, ...
    'initial_state', initial_state, 'final_state', final_state);
save(fullfile(outdir, 'candidate_confirmation_trajectory.mat'), 'results', '-v7.3');

metadata = brusselator_run_metadata(cfg, mfilename, runtime);
brusselator_write_metadata(outdir, cfg, metadata);
disp(struct2table(summary, 'AsArray', true));
fprintf('Candidate-C frozen confirmation complete. Results saved in: %s\n', outdir);
end

function section = full_state_poincare_section(S, U, V)
observable = mean(V, 2);
observable = observable - mean(observable);
cross_index = find(observable(1:end-1) <= 0 & observable(2:end) > 0);
times = zeros(numel(cross_index), 1);
states = zeros(numel(cross_index), 2 * size(U,2));
for k = 1:numel(cross_index)
    j = cross_index(k);
    alpha = -observable(j) / (observable(j+1) - observable(j));
    times(k) = S(j) + alpha * (S(j+1) - S(j));
    u = (1 - alpha) * U(j,:) + alpha * U(j+1,:);
    v = (1 - alpha) * V(j,:) + alpha * V(j+1,:);
    states(k,:) = [u, v];
end
if numel(times) < 4
    error('run_edge_candidate_frozen_confirmation:InsufficientCrossings', ...
        'Need at least four late Poincare crossings for recurrence testing.');
end
periods = diff(times);
section = struct('times', times, 'states', states, 'periods', periods, ...
    'period', median(periods), ...
    'period_cv', std(periods) / max(mean(periods), eps));
end

function recurrence = recurrence_scan(section, threshold)
count = size(section.states,1);
max_lag = min(30, count - 2);
record_template = struct('return_crossings', NaN, 'return_time', NaN, ...
    'direct_median_rel_diff_full_state', NaN, ...
    'reflection_median_rel_diff_full_state', NaN, ...
    'best_median_rel_diff_full_state', NaN, 'best_uses_reflection', false, ...
    'direct_max_rel_diff_full_state', NaN, ...
    'reflection_max_rel_diff_full_state', NaN);
records = repmat(record_template, max_lag, 1);
for lag = 1:max_lag
    first = section.states(1:end-lag,:);
    second = section.states(1+lag:end,:);
    direct = rowwise_relative_distance(second, first);
    reflected = rowwise_relative_distance(reflect_full_states(second), first);
    record = record_template;
    record.return_crossings = lag;
    record.return_time = median(section.times(1+lag:end) - section.times(1:end-lag));
    record.direct_median_rel_diff_full_state = median(direct);
    record.reflection_median_rel_diff_full_state = median(reflected);
    if record.reflection_median_rel_diff_full_state < ...
            record.direct_median_rel_diff_full_state
        record.best_median_rel_diff_full_state = ...
            record.reflection_median_rel_diff_full_state;
        record.best_uses_reflection = true;
    else
        record.best_median_rel_diff_full_state = ...
            record.direct_median_rel_diff_full_state;
    end
    record.direct_max_rel_diff_full_state = max(direct);
    record.reflection_max_rel_diff_full_state = max(reflected);
    records(lag) = record;
end
[~, best_index] = min([records.best_median_rel_diff_full_state]);
best = records(best_index);
close_index = find([records.best_median_rel_diff_full_state] < threshold, 1);
if isempty(close_index)
    first_close = struct('return_crossings', NaN, 'return_time', NaN, ...
        'best_median_rel_diff_full_state', NaN, 'best_uses_reflection', false);
else
    first_close = records(close_index);
end
recurrence = struct('records', records, 'best', best, 'first_close', first_close);
end

function label = classify_confirmation(section, recurrence)
if recurrence.first_close.return_crossings > 0 && section.period_cv < 0.05
    if recurrence.first_close.best_uses_reflection
        label = ['fresh_restart_supports_reflection_related_candidate_periodic_state; ', ...
            'not a stability proof'];
    else
        label = ['fresh_restart_supports_direct_candidate_periodic_state; ', ...
            'not a stability proof'];
    end
elseif isfinite(section.period_cv) && section.period_cv < 0.02
    label = ['fresh_restart_has_regular_section_timing_without_close_recurrence; ', ...
        'inconclusive'];
else
    label = ['fresh_restart_does_not_confirm_close_recurrence; ', ...
        'inconclusive'];
end
end

function distances = rowwise_relative_distance(A, reference)
numerator = sqrt(sum((A - reference).^2, 2));
denominator = sqrt(sum(reference.^2, 2));
distances = numerator ./ max(denominator, eps);
end

function reflected = reflect_full_states(states)
N = size(states,2) / 2;
if N ~= floor(N)
    error('run_edge_candidate_frozen_confirmation:InvalidStateDimension', ...
        'Full states must concatenate equally sized u and v fields.');
end
reflected = [fliplr(states(:,1:N)), fliplr(states(:,N+1:end))];
end

function make_figures(outdir, S, V, recurrence, cfg)
fig = figure('Color', 'w', 'Position', [100 100 900 500]);
records = recurrence.records;
semilogy([records.return_crossings], ...
    max([records.direct_median_rel_diff_full_state], eps), '-o', ...
    'LineWidth', 1.15, 'DisplayName', 'direct'); hold on;
semilogy([records.return_crossings], ...
    max([records.reflection_median_rel_diff_full_state], eps), '-o', ...
    'LineWidth', 1.15, 'DisplayName', 'reflection');
yline(cfg.classification.recovery_threshold, '--', 'recurrence threshold', ...
    'HandleVisibility', 'off');
xlabel('Poincare return crossings');
ylabel('Median relative full-state return distance');
title('Fresh frozen restart: candidate-C late-window recurrence');
legend('Location', 'best'); grid on; box on;
exportgraphics(fig, fullfile(outdir, 'multi_return_recurrence.png'), ...
    'Resolution', 300);
close(fig);

fig = figure('Color', 'w', 'Position', [100 100 950 500]);
imagesc(linspace(0, cfg.grid.Lx, size(V,2)), S, V); axis xy;
xlabel('x'); ylabel('fresh frozen b=10 time');
title('Candidate-C confirmation: v(x,t)'); colorbar; box on;
exportgraphics(fig, fullfile(outdir, 'candidate_confirmation_v_spacetime.png'), ...
    'Resolution', 300);
close(fig);
end
