function results = analyze_phase_aware_attractors(input_dir, output_dir)
%ANALYZE_PHASE_AWARE_ATTRACTORS Test whether late trajectories share an orbit.
%   RESULTS = ANALYZE_PHASE_AWARE_ATTRACTORS() analyzes the saved trajectories
%   from EXTENDED_POSTFORCING_REFERENCE_V2 without rerunning a simulation. It
%   minimizes the late-time v-field distance over temporal lag and spatial
%   reflection, then compares temporal periods, time-averaged profiles,
%   cosine-mode power spectra, and Poincare recurrences.
%
%   The available reference output stores v(x,t), not u(x,t). Therefore this
%   is a deliberately conservative first discriminator: an apparent match is
%   evidence consistent with a common orbit but must later be confirmed using
%   the complete (u,v) state. A persistent v-only separation is nevertheless
%   sufficient to rule out exact equality of the full states.
%
%   Optional INPUT_DIR is the directory containing raw_extended_trajectories.mat.
%   OUTPUT_DIR defaults to experiment_outputs/phase_aware_attractor_analysis_v3
%   and must not already exist. This avoids overwriting previous analyses.

if nargin < 1 || isempty(input_dir)
    input_dir = fullfile('experiment_outputs', 'extended_postforcing_reference_v2');
end
if nargin < 2 || isempty(output_dir)
    output_dir = fullfile('experiment_outputs', 'phase_aware_attractor_analysis_v3');
end

raw_file = fullfile(input_dir, 'raw_extended_trajectories.mat');
if ~isfile(raw_file)
    error('analyze_phase_aware_attractors:MissingRawData', ...
        'Could not find %s.', raw_file);
end
if isfolder(output_dir)
    error('analyze_phase_aware_attractors:OutputExists', ...
        'Refusing to overwrite existing analysis directory: %s', output_dir);
end
mkdir(output_dir);

timer = tic;
loaded = load(raw_file, 'raw');
raw = loaded.raw;
cfg = raw.config;

analysis = struct();
analysis.input_dir = input_dir;
analysis.input_file = raw_file;
analysis.component = 'v';
analysis.late_window = 1000;
analysis.minimum_period = 2.01;
analysis.maximum_period = 250;
analysis.fine_lag_spacing = 0.05;
analysis.window_lengths = [100, 250, 500, 1000];
analysis.same_orbit_threshold = cfg.classification.recovery_threshold;
analysis.distinct_threshold = cfg.classification.transition_threshold;
analysis.period_relative_tolerance = 0.05;
analysis.stationary_relative_variation_threshold = 1e-4;
analysis.minimum_samples_per_period = 10;

if cfg.time.Tfinal - cfg.time.forcing_end < analysis.late_window
    error('analyze_phase_aware_attractors:InsufficientPostForcingData', ...
        'The saved trajectory has less than %.0f post-forcing time units.', ...
        analysis.late_window);
end

S = raw.control.S;
late_start = cfg.time.Tfinal - analysis.late_window;
late_idx = S >= late_start;
time = S(late_idx);
V_control = raw.control.V(late_idx,:);

control_dynamics = characterize_temporal_dynamics(time, V_control, analysis);
control_period = control_dynamics.period;
control_section = control_dynamics.section;

case_count = numel(raw.forced_trajectories);
record_template = struct( ...
    'bmax', NaN, ...
    'window_start', time(1), ...
    'window_end', time(end), ...
    'control_period_fft', control_period, ...
    'forced_period_fft', NaN, ...
    'control_relative_temporal_variation', control_dynamics.relative_variation, ...
    'forced_relative_temporal_variation', NaN, ...
    'control_temporal_class', control_dynamics.temporal_class, ...
    'forced_temporal_class', 'unclassified', ...
    'forced_samples_per_estimated_period', NaN, ...
    'control_period_poincare', control_dynamics.section_period, ...
    'forced_period_poincare', NaN, ...
    'period_relative_difference', NaN, ...
    'same_time_rel_diff_v', NaN, ...
    'phase_aware_rel_diff_v', NaN, ...
    'best_lag', NaN, ...
    'best_reflection', false, ...
    'mean_profile_rel_diff_v', NaN, ...
    'profile_uses_reflection', false, ...
    'spatial_spectrum_rel_diff_v', NaN, ...
    'control_poincare_recurrence_rel_diff_v', control_section.recurrence, ...
    'forced_poincare_recurrence_rel_diff_v', NaN, ...
    'poincare_cross_distance_rel_diff_v', NaN, ...
    'interpretation', 'unresolved');
records = repmat(record_template, case_count, 1);
window_records = struct('bmax', {}, 'window_length', {}, ...
    'phase_aware_rel_diff_v', {}, 'best_lag', {}, 'best_reflection', {});
section_records = struct('trajectory', {}, 'bmax', {}, 'crossing_time', {});
section_records = append_section_records(section_records, ...
    'control', NaN, control_section.times);
case_results = cell(case_count, 1);

for k = 1:case_count
    forced = raw.forced_trajectories{k};
    if ~isequal(forced.S, S)
        error('analyze_phase_aware_attractors:TimeMismatch', ...
            'Forced and control output times must agree.');
    end
    V_forced = forced.V(late_idx,:);
    forced_dynamics = characterize_temporal_dynamics(time, V_forced, analysis);
    forced_period = forced_dynamics.period;
    forced_section = forced_dynamics.section;

    max_lag = choose_max_lag(control_period, forced_period, time, analysis);
    alignment = phase_reflection_alignment(time, V_forced, V_control, ...
        max_lag, analysis.fine_lag_spacing);
    profile = compare_spatial_statistics(V_forced, V_control, cfg.grid.Lx);
    section_distance = poincare_cross_distance(control_section.states, ...
        forced_section.states);
    period_difference = relative_period_difference(control_period, forced_period);

    record = record_template;
    record.bmax = forced.bmax;
    record.forced_period_fft = forced_period;
    record.forced_relative_temporal_variation = forced_dynamics.relative_variation;
    record.forced_temporal_class = forced_dynamics.temporal_class;
    record.forced_samples_per_estimated_period = forced_dynamics.samples_per_period;
    record.forced_period_poincare = forced_dynamics.section_period;
    record.period_relative_difference = period_difference;
    record.same_time_rel_diff_v = alignment.same_time_distance;
    record.phase_aware_rel_diff_v = alignment.best_distance;
    record.best_lag = alignment.best_lag;
    record.best_reflection = alignment.best_reflection;
    record.mean_profile_rel_diff_v = profile.mean_profile_distance;
    record.profile_uses_reflection = profile.uses_reflection;
    record.spatial_spectrum_rel_diff_v = profile.spectrum_distance;
    record.forced_poincare_recurrence_rel_diff_v = forced_section.recurrence;
    record.poincare_cross_distance_rel_diff_v = section_distance;
    record.interpretation = interpret_case(record, analysis);
    records(k) = record;

    for w = analysis.window_lengths
        if w > time(end) - time(1)
            continue
        end
        window_idx = time >= time(end) - w;
        window_alignment = phase_reflection_alignment(time(window_idx), ...
            V_forced(window_idx,:), V_control(window_idx,:), max_lag, ...
            analysis.fine_lag_spacing);
        window_records(end+1) = struct( ... %#ok<AGROW>
            'bmax', forced.bmax, 'window_length', w, ...
            'phase_aware_rel_diff_v', window_alignment.best_distance, ...
            'best_lag', window_alignment.best_lag, ...
            'best_reflection', window_alignment.best_reflection);
    end

    section_records = append_section_records(section_records, ...
        'forced', forced.bmax, forced_section.times);

    case_results{k} = struct( ...
        'bmax', forced.bmax, ...
        'alignment', alignment, ...
        'profile', profile, ...
        'control_dynamics', control_dynamics, ...
        'forced_dynamics', forced_dynamics, ...
        'forced_section', forced_section);
    make_case_figures(output_dir, time, V_control, V_forced, control_dynamics, ...
        forced_dynamics, control_section, ...
        forced_section, alignment, profile, forced.bmax, late_start, cfg.grid.Lx);
end

summary_table = struct2table(records, 'AsArray', true);
writetable(summary_table, fullfile(output_dir, 'phase_aware_summary.csv'));
if ~isempty(window_records)
    writetable(struct2table(window_records, 'AsArray', true), ...
        fullfile(output_dir, 'windowed_phase_alignment.csv'));
end
if ~isempty(section_records)
    writetable(struct2table(section_records, 'AsArray', true), ...
        fullfile(output_dir, 'poincare_crossings.csv'));
end

results = struct();
results.analysis = analysis;
results.control_dynamics = control_dynamics;
results.control_section = control_section;
results.summary = records;
results.windowed_alignment = window_records;
results.case_results = case_results;
save(fullfile(output_dir, 'phase_aware_analysis.mat'), 'results', '-v7.3');

analysis_cfg = cfg;
analysis_cfg.experiment_name = 'phase_aware_attractor_analysis_v3';
analysis_cfg.analysis = analysis;
metadata = brusselator_run_metadata(analysis_cfg, mfilename, toc(timer));
brusselator_write_metadata(output_dir, analysis_cfg, metadata);

disp(summary_table);
fprintf(['Phase-aware v-field analysis complete. This is not a full-state ', ...
    'attractor identification; results saved in: %s\n'], output_dir);
end

function max_lag = choose_max_lag(control_period, forced_period, time, analysis)
candidate_periods = [control_period, forced_period];
candidate_periods = candidate_periods(isfinite(candidate_periods));
if isempty(candidate_periods)
    max_lag = analysis.maximum_period;
else
    max_lag = ceil(1.25 * max(candidate_periods));
end
max_lag = min(max_lag, floor((time(end) - time(1)) / 4));
max_lag = max(max_lag, analysis.minimum_period);
end

function alignment = phase_reflection_alignment(time, V_forced, V_control, max_lag, fine_step)
dt = median(diff(time));
coarse_lags = (-max_lag:dt:max_lag)';
coarse_no_reflection = zeros(size(coarse_lags));
coarse_reflection = zeros(size(coarse_lags));
for j = 1:numel(coarse_lags)
    coarse_no_reflection(j) = shifted_distance(time, V_forced, V_control, ...
        coarse_lags(j), false);
    coarse_reflection(j) = shifted_distance(time, V_forced, V_control, ...
        coarse_lags(j), true);
end

[minimum_no_reflection, idx_no_reflection] = min(coarse_no_reflection);
[minimum_reflection, idx_reflection] = min(coarse_reflection);
if minimum_reflection < minimum_no_reflection
    coarse_best_lag = coarse_lags(idx_reflection);
    coarse_reflection_flag = true;
else
    coarse_best_lag = coarse_lags(idx_no_reflection);
    coarse_reflection_flag = false;
end

fine_lags = unique((coarse_best_lag - dt):fine_step:(coarse_best_lag + dt));
fine_lags = fine_lags(fine_lags >= -max_lag & fine_lags <= max_lag);
fine_distances = zeros(size(fine_lags));
for j = 1:numel(fine_lags)
    fine_distances(j) = shifted_distance(time, V_forced, V_control, ...
        fine_lags(j), coarse_reflection_flag);
end
[best_distance, best_idx] = min(fine_distances);

alignment = struct();
alignment.coarse_lags = coarse_lags;
alignment.coarse_no_reflection = coarse_no_reflection;
alignment.coarse_reflection = coarse_reflection;
alignment.fine_lags = fine_lags;
alignment.fine_distances = fine_distances;
alignment.same_time_distance = shifted_distance(time, V_forced, V_control, 0, false);
alignment.best_distance = best_distance;
alignment.best_lag = fine_lags(best_idx);
alignment.best_reflection = coarse_reflection_flag;
end

function distance = shifted_distance(time, V_forced, V_control, lag, reflect_control)
valid = (time + lag >= time(1)) & (time + lag <= time(end));
forced_values = V_forced(valid,:);
control_values = interp1(time, V_control, time(valid) + lag, 'linear');
if reflect_control
    control_values = fliplr(control_values);
end
numerator = sum((forced_values - control_values).^2, 'all');
denominator = sum(control_values.^2, 'all');
distance = sqrt(numerator / max(denominator, eps));
end

function dynamics = characterize_temporal_dynamics(time, V, analysis)
% Use a physical spatial average when it is time-dependent; otherwise use
% the leading spatial-PCA coordinate only to visualize residual motion.
% A nearly stationary field has no meaningful temporal phase or section.
dt = median(diff(time));
W = V - mean(V, 1);
relative_variation = sqrt(mean(sum(W.^2, 2))) / ...
    max(sqrt(mean(sum(V.^2, 2))), eps);
q = mean(V, 2);
q = q - mean(q);
if std(q) <= eps
    [~,~,spatial_vectors] = svd(W, 'econ');
    q = W * spatial_vectors(:,1);
    q = q - mean(q);
end

[period, spectrum] = dominant_period_from_signal(time, q, ...
    analysis.minimum_period, analysis.maximum_period);

section = struct('times', zeros(0,1), 'states', zeros(0, size(V,2)), ...
    'period', NaN, 'recurrence', NaN);
section_period = NaN;
samples_per_period = NaN;
if relative_variation < analysis.stationary_relative_variation_threshold
    temporal_class = 'numerically_stationary_v';
    period = NaN;
else
    candidate_section = poincare_recurrence(time, V, q);
    section_period = candidate_section.period;
    if isfinite(section_period)
        samples_per_period = section_period / dt;
    elseif isfinite(period)
        samples_per_period = period / dt;
    end
    if samples_per_period < analysis.minimum_samples_per_period
        temporal_class = 'time_dependent_v_undersampled';
    else
        temporal_class = 'time_dependent_v_resolved';
        section = candidate_section;
    end
end

dynamics = struct('q', q, 'spectrum', spectrum, 'period', period, ...
    'section_period', section_period, 'samples_per_period', samples_per_period, ...
    'relative_variation', relative_variation, ...
    'temporal_class', temporal_class, 'section', section);
end

function [period, spectrum] = dominant_period_from_signal(time, q, minimum_period, maximum_period)
dt = median(diff(time));
n = numel(q);
transform = fft(q);
positive_idx = 1:(floor(n/2) + 1);
frequencies = (positive_idx - 1)' / (n * dt);
power = abs(transform(positive_idx)).^2;
valid = frequencies > 0 & frequencies >= 1 / maximum_period & ...
    frequencies <= 1 / minimum_period;
period = NaN;
if any(valid) && any(power(valid) > 0)
    valid_idx = find(valid);
    [~, local_idx] = max(power(valid));
    frequency_idx = valid_idx(local_idx);
    period = 1 / frequencies(frequency_idx);
end

if sum(power(2:end)) > 0
    normalized_power = power / sum(power(2:end));
else
    normalized_power = power;
end
spectrum = struct('frequency', frequencies, 'power', normalized_power);
end

function section = poincare_recurrence(time, V, q)
cross_idx = find(q(1:end-1) <= 0 & q(2:end) > 0);
crossing_times = zeros(numel(cross_idx), 1);
crossing_states = zeros(numel(cross_idx), size(V,2));
for j = 1:numel(cross_idx)
    i = cross_idx(j);
    alpha = -q(i) / (q(i+1) - q(i));
    crossing_times(j) = time(i) + alpha * (time(i+1) - time(i));
    crossing_states(j,:) = (1 - alpha) * V(i,:) + alpha * V(i+1,:);
end
if numel(crossing_times) >= 2
    period = median(diff(crossing_times));
else
    period = NaN;
end
if size(crossing_states,1) >= 2
    differences = sqrt(sum(diff(crossing_states).^2, 2));
    norms = sqrt(sum(crossing_states(2:end,:).^2, 2));
    recurrence = median(differences ./ max(norms, eps));
else
    recurrence = NaN;
end
section = struct('times', crossing_times, 'states', crossing_states, ...
    'period', period, 'recurrence', recurrence);
end

function profile = compare_spatial_statistics(V_forced, V_control, Lx)
mean_forced = mean(V_forced, 1);
mean_control = mean(V_control, 1);
profile_direct = relative_vector_distance(mean_forced, mean_control);
profile_reflected = relative_vector_distance(mean_forced, fliplr(mean_control));
if profile_reflected < profile_direct
    profile_distance = profile_reflected;
    uses_reflection = true;
else
    profile_distance = profile_direct;
    uses_reflection = false;
end

mode_count = min(30, floor((size(V_control,2) - 1) / 2));
grid = linspace(0, Lx, size(V_control,2));
mode_numbers = 0:mode_count;
basis = cos((mode_numbers' * pi / Lx) * grid);
weights = ones(1, numel(grid));
weights([1,end]) = 0.5;
weighted_basis = basis .* weights;
power_control = mean((V_control * weighted_basis').^2, 1);
power_forced = mean((V_forced * weighted_basis').^2, 1);
power_control = power_control / max(sum(power_control), eps);
power_forced = power_forced / max(sum(power_forced), eps);

profile = struct();
profile.mean_forced = mean_forced;
profile.mean_control = mean_control;
profile.mean_profile_distance = profile_distance;
profile.uses_reflection = uses_reflection;
profile.mode_numbers = mode_numbers;
profile.power_control = power_control;
profile.power_forced = power_forced;
profile.spectrum_distance = relative_vector_distance(power_forced, power_control);
end

function distance = relative_vector_distance(A, B)
distance = norm(A - B) / max(norm(B), eps);
end

function distance = poincare_cross_distance(states_A, states_B)
if isempty(states_A) || isempty(states_B)
    distance = NaN;
    return
end
states_A = states_A(max(1, end-9):end,:);
states_B = states_B(max(1, end-9):end,:);
distance = inf;
for j = 1:size(states_A,1)
    for k = 1:size(states_B,1)
        direct = relative_vector_distance(states_A(j,:), states_B(k,:));
        reflected = relative_vector_distance(states_A(j,:), fliplr(states_B(k,:)));
        distance = min([distance, direct, reflected]);
    end
end
end

function difference = relative_period_difference(period_A, period_B)
if ~(isfinite(period_A) && isfinite(period_B))
    difference = NaN;
    return
end
difference = abs(period_A - period_B) / max([period_A, period_B, eps]);
end

function label = interpret_case(record, analysis)
if strcmp(record.control_temporal_class, 'numerically_stationary_v')
    if record.phase_aware_rel_diff_v < analysis.same_orbit_threshold
        label = 'consistent_with_baseline_stationary_v_state';
    elseif ~strcmp(record.forced_temporal_class, 'numerically_stationary_v') && ...
            record.phase_aware_rel_diff_v > analysis.distinct_threshold
        label = 'candidate_distinct_time_dependent_v_outcome_requires_fine_sampling';
    else
        label = 'unresolved_v_only_dynamics';
    end
    return
end

is_recurrent = record.control_poincare_recurrence_rel_diff_v < ...
        analysis.same_orbit_threshold && ...
    record.forced_poincare_recurrence_rel_diff_v < analysis.same_orbit_threshold;
period_matches = isfinite(record.period_relative_difference) && ...
    record.period_relative_difference <= analysis.period_relative_tolerance;
if is_recurrent && period_matches && ...
        record.phase_aware_rel_diff_v < analysis.same_orbit_threshold
    label = 'consistent_with_same_v_orbit_up_to_phase_or_reflection';
elseif is_recurrent && record.phase_aware_rel_diff_v > analysis.distinct_threshold
    label = 'candidate_distinct_recurrent_v_outcome';
else
    label = 'unresolved_v_only_dynamics';
end
end

function records = append_section_records(records, trajectory, bmax, times)
for j = 1:numel(times)
    records(end+1) = struct( ... %#ok<AGROW>
        'trajectory', trajectory, 'bmax', bmax, 'crossing_time', times(j));
end
end

function make_case_figures(output_dir, time, V_control, V_forced, control_dynamics, ...
        forced_dynamics, control_section, ...
        forced_section, alignment, profile, bmax, late_start, Lx)
tag = sprintf('bmax_%0.2f', bmax);
control_q = control_dynamics.q;
forced_q = forced_dynamics.q;
control_spectrum = control_dynamics.spectrum;
forced_spectrum = forced_dynamics.spectrum;

fig = figure('Color', 'w', 'Position', [100 100 1050 750]);
subplot(2,2,1);
plot(alignment.coarse_lags, alignment.coarse_no_reflection, 'LineWidth', 1.25, ...
    'DisplayName', 'identity'); hold on;
plot(alignment.coarse_lags, alignment.coarse_reflection, 'LineWidth', 1.25, ...
    'DisplayName', 'reflection');
xline(alignment.best_lag, '--', 'best lag', 'HandleVisibility', 'off');
xlabel('Control time lag'); ylabel('Relative v-field distance');
title('Phase/reflection alignment'); legend('Location', 'best'); grid on;

subplot(2,2,2);
x = linspace(0, Lx, size(V_control,2));
plot(x, profile.mean_control, 'LineWidth', 1.3, 'DisplayName', 'control'); hold on;
plot(x, profile.mean_forced, 'LineWidth', 1.3, 'DisplayName', 'forced');
if profile.uses_reflection
    plot(x, fliplr(profile.mean_control), '--', 'LineWidth', 1.1, ...
        'DisplayName', 'reflected control');
end
xlabel('x'); ylabel('Time-averaged v');
title('Late-time mean spatial profile'); legend('Location', 'best'); grid on;

subplot(2,2,3);
plot(control_spectrum.frequency, control_spectrum.power, 'LineWidth', 1.25, ...
    'DisplayName', 'control'); hold on;
plot(forced_spectrum.frequency, forced_spectrum.power, 'LineWidth', 1.25, ...
    'DisplayName', 'forced');
xlim([0, 1/5]); xlabel('Temporal frequency'); ylabel('Normalized power');
title('Dominant temporal spectra'); legend('Location', 'best'); grid on;

subplot(2,2,4);
semilogy(profile.mode_numbers, max(profile.power_control, eps), '-o', ...
    'LineWidth', 1.1, 'DisplayName', 'control'); hold on;
semilogy(profile.mode_numbers, max(profile.power_forced, eps), '-o', ...
    'LineWidth', 1.1, 'DisplayName', 'forced');
xlabel('Neumann cosine mode'); ylabel('Normalized mean power');
title('Spatial cosine-mode power'); legend('Location', 'best'); grid on;
sgtitle(sprintf('Late-window phase-aware comparison, b_{max}=%.2f (t >= %.0f)', ...
    bmax, late_start));
exportgraphics(fig, fullfile(output_dir, ['phase_aware_comparison_', tag, '.png']), ...
    'Resolution', 300);
close(fig);

fig = figure('Color', 'w', 'Position', [100 100 1000 500]);
subplot(1,2,1);
plot(time, control_q, 'LineWidth', 1.0, 'DisplayName', 'control'); hold on;
plot(time, forced_q, 'LineWidth', 1.0, 'DisplayName', 'forced');
xlabel('t'); ylabel('Centered spatial-average v (or fallback coordinate)');
title('Physical temporal observable');
legend('Location', 'best'); grid on;

subplot(1,2,2);
if isempty(control_section.times) && isempty(forced_section.times)
    axis off;
    text(0.5, 0.55, 'No Poincare section computed', ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    text(0.5, 0.40, sprintf('control: %s\nforced: %s', ...
        control_dynamics.temporal_class, forced_dynamics.temporal_class), ...
        'HorizontalAlignment', 'center');
else
    plot(control_section.times, zeros(size(control_section.times)), 'o', ...
        'DisplayName', 'control crossings'); hold on;
    plot(forced_section.times, ones(size(forced_section.times)), 'o', ...
        'DisplayName', 'forced crossings');
    ylim([-0.5, 1.5]); yticks([0 1]); yticklabels({'control', 'forced'});
    xlabel('t'); title('Positive Poincare-section crossings');
    legend('Location', 'best'); grid on;
end
sgtitle(sprintf('Poincare recurrence diagnostic, b_{max}=%.2f', bmax));
exportgraphics(fig, fullfile(output_dir, ['poincare_recurrence_', tag, '.png']), ...
    'Resolution', 300);
close(fig);
end
