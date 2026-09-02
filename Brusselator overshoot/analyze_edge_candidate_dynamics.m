function results = analyze_edge_candidate_dynamics()
%ANALYZE_EDGE_CANDIDATE_DYNAMICS Characterize the persistent slice outcome.
%   Analyzes the saved lambda = 0.8 frozen-b=10 trajectory that remained far
%   from both established endpoint neighborhoods through time 240. The
%   analysis uses the late window only and asks whether it has a close
%   full-state Poincare recurrence, including the spatial reflection symmetry.
%   It also saves temporal spectral features and a two-dimensional Poincare
%   projection. This is an outcome-characterization diagnostic, not evidence
%   of another attractor unless recurrence and longer independent integrations
%   subsequently support that conclusion.

study_id = 'edge_candidate_dynamics_v1';
parent_file = fullfile('experiment_outputs', ...
    'phase_aware_basin_slice_edge_extension_v1', ...
    'raw_edge_extension_trajectories.mat');
outdir = fullfile('experiment_outputs', study_id);
if ~isfile(parent_file)
    error('analyze_edge_candidate_dynamics:MissingParentData', ...
        'Expected edge-extension data at %s.', parent_file);
end
if isfolder(outdir)
    error('analyze_edge_candidate_dynamics:OutputExists', ...
        'Refusing to overwrite existing experiment directory: %s', outdir);
end

timer = tic;
loaded = load(parent_file, 'results');
parent = loaded.results;
candidate_index = find(abs([parent.records.lambda] - 0.8) < 1e-12, 1);
if isempty(candidate_index)
    error('analyze_edge_candidate_dynamics:MissingLambdaPoint', ...
        'The edge-extension output contains no lambda = 0.8 trajectory.');
end
if ~strcmp(parent.records(candidate_index).classification, ...
        'persistent_unresolved_edge_candidate')
    error('analyze_edge_candidate_dynamics:UnexpectedCandidateClass', ...
        'The lambda = 0.8 trajectory is not marked as a persistent unresolved case.');
end

cfg = parent.configuration;
cfg.experiment_name = study_id;
cfg.edge_dynamics = struct( ...
    'parent_file', parent_file, ...
    'lambda', parent.records(candidate_index).lambda, ...
    'analysis_start_time', 160, ...
    'analysis_end_time', parent.records(candidate_index).total_duration, ...
    'observable_for_section', 'spatial_mean_v minus late-window mean', ...
    'section_direction', 'upward_crossing', ...
    'symmetries', 'identity, spatial reflection x -> L-x', ...
    'recurrence_threshold', cfg.classification.recovery_threshold, ...
    'purpose', ['finite-time recurrence and spectral characterization of ', ...
        'the unresolved slice trajectory']);

mkdir(outdir);
brusselator_write_metadata(outdir, cfg, ...
    brusselator_run_metadata(cfg, mfilename, NaN));

data = parent.trajectory_data{candidate_index};
[S, U, V] = join_trajectory_segments(data);
late = S >= cfg.edge_dynamics.analysis_start_time & ...
    S <= cfg.edge_dynamics.analysis_end_time;
if sum(late) < 20
    error('analyze_edge_candidate_dynamics:InsufficientLateData', ...
        'The requested late window has too few saved output states.');
end
S_late = S(late);
U_late = U(late,:);
V_late = V(late,:);
dt = median(diff(S_late));
if max(abs(diff(S_late) - dt)) > 100 * eps(max(S_late))
    error('analyze_edge_candidate_dynamics:NonuniformOutput', ...
        'Late trajectory output must be uniformly spaced for the spectral analysis.');
end

section = full_state_poincare_section(S_late, U_late, V_late);
recurrence = recurrence_scan(section, cfg.edge_dynamics.recurrence_threshold);
spectra = spectral_features(S_late, U_late, V_late, cfg.grid.Lx);
projection = section_projection(section.states);
classification = classify_candidate(section, recurrence);

summary = struct();
summary.lambda = cfg.edge_dynamics.lambda;
summary.analysis_start_time = cfg.edge_dynamics.analysis_start_time;
summary.analysis_end_time = cfg.edge_dynamics.analysis_end_time;
summary.analysis_duration = S_late(end) - S_late(1);
summary.output_dt = dt;
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
summary.best_return_time = recurrence.best.return_time;
summary.best_return_rel_diff_full_state = ...
    recurrence.best.best_median_rel_diff_full_state;
summary.best_return_uses_reflection = recurrence.best.best_uses_reflection;
summary.dominant_mean_v_frequency = spectra(2).dominant_frequency;
summary.dominant_mean_v_period = spectra(2).dominant_period;
summary.interpretation = classification;

writetable(struct2table(summary, 'AsArray', true), ...
    fullfile(outdir, 'edge_candidate_dynamics_summary.csv'));
writetable(struct2table(recurrence.records, 'AsArray', true), ...
    fullfile(outdir, 'multi_return_recurrence.csv'));
writetable(struct2table(spectra, 'AsArray', true), ...
    fullfile(outdir, 'temporal_spectral_features.csv'));
write_section_table(outdir, section, projection);
make_figures(outdir, S_late, U_late, V_late, section, recurrence, ...
    spectra, projection, cfg);

results = struct();
results.configuration = cfg;
results.parent_record = parent.records(candidate_index);
results.summary = summary;
results.section = section;
results.recurrence = recurrence;
results.spectra = spectra;
results.projection = projection;
save(fullfile(outdir, 'edge_candidate_dynamics.mat'), 'results', '-v7.3');

metadata = brusselator_run_metadata(cfg, mfilename, toc(timer));
brusselator_write_metadata(outdir, cfg, metadata);
disp(struct2table(summary, 'AsArray', true));
fprintf('Edge-candidate dynamics analysis complete. Results saved in: %s\n', outdir);
end

function [S, U, V] = join_trajectory_segments(data)
S = [data.old_S; data.old_S(end) + data.continuation_S(2:end)];
U = [data.old_U; data.continuation_U(2:end,:)];
V = [data.old_V; data.continuation_V(2:end,:)];
if size(U,1) ~= numel(S) || size(V,1) ~= numel(S)
    error('analyze_edge_candidate_dynamics:SegmentMismatch', ...
        'Saved trajectory segment dimensions are inconsistent.');
end
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
if numel(times) >= 2
    periods = diff(times);
    period = median(periods);
    period_cv = std(periods) / max(mean(periods), eps);
else
    periods = zeros(0,1);
    period = NaN;
    period_cv = NaN;
end
section = struct('observable', observable, 'times', times, 'states', states, ...
    'periods', periods, 'period', period, 'period_cv', period_cv);
end

function recurrence = recurrence_scan(section, threshold)
count = size(section.states, 1);
if count < 4
    error('analyze_edge_candidate_dynamics:InsufficientCrossings', ...
        'Need at least four Poincare crossings for a recurrence scan.');
end
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

function spectra = spectral_features(S, U, V, Lx)
dt = median(diff(S));
x = linspace(0, Lx, size(U,2));
features = struct();
features(1).name = 'spatial_mean_u';
features(1).signal = mean(U, 2);
features(2).name = 'spatial_mean_v';
features(2).signal = mean(V, 2);
features(3).name = 'u_neumann_cosine_mode_1';
features(3).signal = mean(U .* cos(pi * x / Lx), 2);
features(4).name = 'v_neumann_cosine_mode_1';
features(4).signal = mean(V .* cos(pi * x / Lx), 2);
features(5).name = 'v_neumann_cosine_mode_2';
features(5).signal = mean(V .* cos(2 * pi * x / Lx), 2);

record_template = struct('feature', '', 'dominant_frequency', NaN, ...
    'dominant_period', NaN, 'dominant_relative_power', NaN, ...
    'relative_temporal_variation', NaN);
spectra = repmat(record_template, numel(features), 1);
for k = 1:numel(features)
    signal = features(k).signal(:);
    centered = signal - mean(signal);
    count = numel(centered);
    Fourier = fft(centered);
    power = abs(Fourier).^2;
    positive = 2:(floor(count/2) + 1);
    [peak_power, local_index] = max(power(positive));
    peak_index = positive(local_index);
    frequency = (peak_index - 1) / (count * dt);
    record = record_template;
    record.feature = features(k).name;
    record.dominant_frequency = frequency;
    record.dominant_period = 1 / max(frequency, eps);
    record.dominant_relative_power = peak_power / max(sum(power(positive)), eps);
    record.relative_temporal_variation = std(signal) / max(rms(signal), eps);
    spectra(k) = record;
end
end

function projection = section_projection(states)
if isempty(states)
    projection = struct('scores', zeros(0,2), 'explained_variance', zeros(0,1));
    return;
end
centered = states - mean(states, 1);
[left, singular_values, ~] = svd(centered, 'econ');
scores = left * singular_values;
variances = diag(singular_values).^2;
explained = variances / max(sum(variances), eps);
projection = struct('scores', scores(:,1:min(2,size(scores,2))), ...
    'explained_variance', explained);
end

function label = classify_candidate(section, recurrence)
if recurrence.first_close.return_crossings > 0
    if recurrence.first_close.best_uses_reflection
        label = ['finite_time_candidate_reflection_related_recurrence; ', ...
            'requires independent confirmation'];
    else
        label = ['finite_time_candidate_direct_recurrence; ', ...
            'requires independent confirmation'];
    end
elseif isfinite(section.period_cv) && section.period_cv < 0.02
    label = ['regular_section_timing_without_close_full_state_recurrence; ', ...
        'not a demonstrated periodic orbit'];
else
    label = ['no_close_full_state_recurrence_in_late_window; ', ...
        'persistent finite_time_unresolved_dynamics'];
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
    error('analyze_edge_candidate_dynamics:InvalidStateDimension', ...
        'Full states must concatenate equally sized u and v fields.');
end
reflected = [fliplr(states(:,1:N)), fliplr(states(:,N+1:end))];
end

function write_section_table(outdir, section, projection)
table_data = table(section.times, 'VariableNames', {'crossing_time'});
if ~isempty(projection.scores)
    table_data.pc1 = projection.scores(:,1);
    if size(projection.scores,2) >= 2
        table_data.pc2 = projection.scores(:,2);
    else
        table_data.pc2 = zeros(height(table_data), 1);
    end
end
writetable(table_data, fullfile(outdir, 'poincare_intersections.csv'));
end

function make_figures(outdir, S, U, V, section, recurrence, spectra, projection, cfg)
fig = figure('Color', 'w', 'Position', [100 100 900 500]);
records = recurrence.records;
semilogy([records.return_crossings], ...
    max([records.direct_median_rel_diff_full_state], eps), '-o', ...
    'LineWidth', 1.15, 'DisplayName', 'direct'); hold on;
semilogy([records.return_crossings], ...
    max([records.reflection_median_rel_diff_full_state], eps), '-o', ...
    'LineWidth', 1.15, 'DisplayName', 'reflection');
yline(cfg.edge_dynamics.recurrence_threshold, '--', 'recurrence threshold', ...
    'HandleVisibility', 'off');
xlabel('Poincare return crossings');
ylabel('Median relative full-state return distance');
title(sprintf('Late-window recurrence: \\lambda = %.3f', cfg.edge_dynamics.lambda));
legend('Location', 'best'); grid on; box on;
exportgraphics(fig, fullfile(outdir, 'multi_return_recurrence.png'), ...
    'Resolution', 300);
close(fig);

fig = figure('Color', 'w', 'Position', [100 100 1000 650]);
tiledlayout(3, 1, 'TileSpacing', 'compact');
nexttile;
plot(S, mean(U,2), 'LineWidth', 1.0); hold on;
plot(S, mean(V,2), 'LineWidth', 1.0);
ylabel('spatial mean'); title('Late-window mean-field time series');
legend('u', 'v', 'Location', 'best'); grid on; box on;
nexttile;
x = linspace(0, cfg.grid.Lx, size(V,2));
plot(S, mean(V .* cos(pi * x / cfg.grid.Lx), 2), 'LineWidth', 1.0); hold on;
plot(S, mean(V .* cos(2 * pi * x / cfg.grid.Lx), 2), 'LineWidth', 1.0);
ylabel('v-mode coefficient'); legend('mode 1', 'mode 2', 'Location', 'best');
grid on; box on;
nexttile;
feature_names = {spectra.feature};
bar([spectra.dominant_frequency]);
xticks(1:numel(feature_names)); xticklabels(feature_names); xtickangle(25);
ylabel('dominant frequency'); grid on; box on;
exportgraphics(fig, fullfile(outdir, 'temporal_features.png'), 'Resolution', 300);
close(fig);

fig = figure('Color', 'w', 'Position', [100 100 800 550]);
if size(projection.scores,1) >= 2 && size(projection.scores,2) >= 2
    scatter(projection.scores(:,1), projection.scores(:,2), 36, section.times, ...
        'filled'); colorbar;
    xlabel(sprintf('PC1 (%.1f%%)', 100 * projection.explained_variance(1)));
    ylabel(sprintf('PC2 (%.1f%%)', 100 * projection.explained_variance(2)));
else
    plot(section.times, projection.scores(:,1), 'o-');
    xlabel('Poincare crossing time'); ylabel('PC1');
end
title('Late-window full-state Poincare intersections'); grid on; box on;
exportgraphics(fig, fullfile(outdir, 'poincare_projection.png'), 'Resolution', 300);
close(fig);

fig = figure('Color', 'w', 'Position', [100 100 950 500]);
imagesc(linspace(0, cfg.grid.Lx, size(V,2)), S, V); axis xy;
xlabel('x'); ylabel('frozen b=10 time');
title(sprintf('Late-window v(x,t): \\lambda = %.3f', cfg.edge_dynamics.lambda));
colorbar; box on;
exportgraphics(fig, fullfile(outdir, 'late_v_spacetime.png'), 'Resolution', 300);
close(fig);
end
