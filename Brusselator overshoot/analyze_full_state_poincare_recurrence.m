function results = analyze_full_state_poincare_recurrence(input_dir, output_dir)
%ANALYZE_FULL_STATE_POINCARE_RECURRENCE Scan multi-return Poincare distances.
%   Tests whether the fine frozen-b=10 departure trajectory repeats after a
%   higher number of spatial-mean section crossings, either directly or after
%   the Neumann-compatible reflection x -> L-x. This distinguishes a simple
%   period-one section return from candidate higher-period/symmetry-related
%   recurrence, while retaining the finite-time caveat for all conclusions.

if nargin < 1 || isempty(input_dir)
    input_dir = fullfile('experiment_outputs', 'frozen_full_state_confirmation_v1');
end
if nargin < 2 || isempty(output_dir)
    output_dir = fullfile('experiment_outputs', 'full_state_poincare_recurrence_v3');
end
input_file = fullfile(input_dir, 'frozen_full_state_trajectories.mat');
if ~isfile(input_file)
    error('analyze_full_state_poincare_recurrence:MissingInput', ...
        'Expected full-state confirmation data at %s.', input_file);
end
if isfolder(output_dir)
    error('analyze_full_state_poincare_recurrence:OutputExists', ...
        'Refusing to overwrite existing analysis directory: %s', output_dir);
end
mkdir(output_dir);

timer = tic;
loaded = load(input_file, 'results');
parent = loaded.results;
section = parent.forced_section;
states = section.states;
if size(states,1) < 12 || ~isfinite(section.period)
    error('analyze_full_state_poincare_recurrence:InsufficientCrossings', ...
        'Need at least 12 forced Poincare crossings for a multi-return scan.');
end

max_lag = min(80, size(states,1) - 10);
record_template = struct('return_crossings', NaN, 'return_time', NaN, ...
    'direct_median_rel_diff_full_state', NaN, ...
    'reflection_median_rel_diff_full_state', NaN, ...
    'best_median_rel_diff_full_state', NaN, 'best_uses_reflection', false, ...
    'direct_max_rel_diff_full_state', NaN, ...
    'reflection_max_rel_diff_full_state', NaN);
records = repmat(record_template, max_lag, 1);
for lag = 1:max_lag
    first = states(1:end-lag,:);
    second = states(1+lag:end,:);
    direct = rowwise_relative_distance(second, first);
    reflected = rowwise_relative_distance(reflect_full_states(second), first);
    record = record_template;
    record.return_crossings = lag;
    record.return_time = lag * section.period;
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

summary_table = struct2table(records, 'AsArray', true);
writetable(summary_table, fullfile(output_dir, 'multi_return_recurrence.csv'));
[best_distance, best_index] = min([records.best_median_rel_diff_full_state]);
best = records(best_index);
close_index = find([records.best_median_rel_diff_full_state] < ...
    parent.configuration.classification.recovery_threshold, 1);
if isempty(close_index)
    first_close = struct('return_crossings', NaN, 'return_time', NaN, ...
        'best_median_rel_diff_full_state', NaN, 'best_uses_reflection', false, ...
        'direct_repeat_crossings', NaN, 'direct_repeat_time', NaN, ...
        'direct_repeat_rel_diff_full_state', NaN);
else
    first_close = records(close_index);
    if first_close.best_uses_reflection && 2 * close_index <= numel(records)
        direct_repeat = records(2 * close_index);
        first_close.direct_repeat_crossings = direct_repeat.return_crossings;
        first_close.direct_repeat_time = direct_repeat.return_time;
        first_close.direct_repeat_rel_diff_full_state = ...
            direct_repeat.direct_median_rel_diff_full_state;
    else
        first_close.direct_repeat_crossings = first_close.return_crossings;
        first_close.direct_repeat_time = first_close.return_time;
        first_close.direct_repeat_rel_diff_full_state = ...
            first_close.direct_median_rel_diff_full_state;
    end
end
classification = classify_recurrence(first_close, parent.configuration.classification);

summary = struct();
summary.section_crossing_period = section.period;
summary.section_crossing_period_cv = section.period_cv;
summary.first_close_return_crossings = first_close.return_crossings;
summary.first_close_return_time = first_close.return_time;
summary.first_close_return_rel_diff_full_state = ...
    first_close.best_median_rel_diff_full_state;
summary.first_close_return_uses_reflection = first_close.best_uses_reflection;
summary.direct_repeat_crossings = first_close.direct_repeat_crossings;
summary.direct_repeat_time = first_close.direct_repeat_time;
summary.direct_repeat_rel_diff_full_state = ...
    first_close.direct_repeat_rel_diff_full_state;
summary.global_minimum_return_crossings = best.return_crossings;
summary.global_minimum_return_rel_diff_full_state = best_distance;
summary.global_minimum_uses_reflection = best.best_uses_reflection;
summary.interpretation = classification;
writetable(struct2table(summary, 'AsArray', true), ...
    fullfile(output_dir, 'recurrence_summary.csv'));

results = struct();
results.input_file = input_file;
results.forced_section_period = section.period;
results.forced_section_period_cv = section.period_cv;
results.best_return = best;
results.first_close_return = first_close;
results.summary = summary;
results.classification = classification;
results.records = records;
save(fullfile(output_dir, 'multi_return_recurrence.mat'), 'results');

fig = figure('Color', 'w', 'Position', [100 100 900 500]);
semilogy([records.return_crossings], ...
    max([records.direct_median_rel_diff_full_state], eps), '-o', ...
    'LineWidth', 1.15, 'DisplayName', 'direct'); hold on;
semilogy([records.return_crossings], ...
    max([records.reflection_median_rel_diff_full_state], eps), '-o', ...
    'LineWidth', 1.15, 'DisplayName', 'reflection');
yline(parent.configuration.classification.recovery_threshold, '--', ...
    'recovery-scale threshold', 'HandleVisibility', 'off');
xline(first_close.return_crossings, '--', 'first close return', ...
    'HandleVisibility', 'off');
xlabel('Number of Poincare crossings');
ylabel('Median relative full-state return distance');
title(sprintf('Multi-return recurrence at frozen b=%.2f', ...
    parent.configuration.confirmation.frozen_b));
legend('Location', 'best'); grid on; box on;
exportgraphics(fig, fullfile(output_dir, 'multi_return_recurrence.png'), ...
    'Resolution', 300);
close(fig);

analysis_cfg = parent.configuration;
analysis_cfg.experiment_name = 'full_state_poincare_recurrence_v3';
analysis_cfg.analysis = struct('input_file', input_file, ...
    'maximum_return_crossings', max_lag, ...
    'symmetries', 'identity, reflection');
metadata = brusselator_run_metadata(analysis_cfg, mfilename, toc(timer));
brusselator_write_metadata(output_dir, analysis_cfg, metadata);

fprintf('First close return: %d crossings (%.6g time), distance %.6g, reflection=%d\n', ...
    first_close.return_crossings, first_close.return_time, ...
    first_close.best_median_rel_diff_full_state, first_close.best_uses_reflection);
fprintf('Global minimum return: %d crossings (%.6g time), distance %.6g, reflection=%d\n', ...
    best.return_crossings, best.return_time, best.best_median_rel_diff_full_state, ...
    best.best_uses_reflection);
fprintf('Classification: %s\n', classification);
end

function distances = rowwise_relative_distance(A, reference)
numerator = sqrt(sum((A - reference).^2, 2));
denominator = sqrt(sum(reference.^2, 2));
distances = numerator ./ max(denominator, eps);
end

function reflected = reflect_full_states(states)
N = size(states,2) / 2;
if N ~= floor(N)
    error('analyze_full_state_poincare_recurrence:InvalidStateDimension', ...
        'Full states must concatenate equally sized u and v fields.');
end
reflected = [fliplr(states(:,1:N)), fliplr(states(:,N+1:end))];
end

function label = classify_recurrence(first_close, classification)
if first_close.best_median_rel_diff_full_state < classification.recovery_threshold
    if first_close.best_uses_reflection && ...
            first_close.direct_repeat_rel_diff_full_state < classification.recovery_threshold
        label = 'candidate_spatiotemporally_reflection_symmetric_periodic_state';
    else
        label = 'candidate_periodic_full_state_recurrence';
    end
elseif isnan(first_close.best_median_rel_diff_full_state)
    label = 'no_close_full_state_recurrence_in_observed_window';
else
    label = 'intermediate_or_unresolved_full_state_recurrence';
end
end
