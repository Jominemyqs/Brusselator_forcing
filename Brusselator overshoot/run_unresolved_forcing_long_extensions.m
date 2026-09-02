function results = run_unresolved_forcing_long_extensions()
%RUN_UNRESOLVED_FORCING_LONG_EXTENSIONS Follow all unresolved cases to T=1000.
%   Collects every unresolved trajectory from the forcing-reachability pilot
%   and ramp refinement, continues it at frozen b=10 to a common 1000-unit
%   post-forcing horizon, and records windowed A/B/C distances, full-state
%   Poincare recurrence, and a mean-field spectral diagnostic.
%
%   A recurrent unresolved trajectory is reported as an intermediate
%   candidate only. Reproducibility across protocols/solvers and invariant-set
%   computation remain separate requirements.

study_id = 'unresolved_forcing_long_extension_v1';
pilot_file = fullfile('experiment_outputs', ...
    'forcing_to_c_reachability_pilot_v1', 'forcing_to_c_reachability.mat');
refinement_file = fullfile('experiment_outputs', ...
    'ramp_transition_refinement_v1', 'ramp_transition_refinement.mat');
outdir = fullfile('experiment_outputs', study_id);
progress_file = fullfile(outdir, 'progress_long_extensions.mat');
if ~isfile(pilot_file) || ~isfile(refinement_file)
    error('run_unresolved_forcing_long_extensions:MissingInput', ...
        'Both the forcing pilot and ramp-refinement results are required.');
end
pilot_data = load(pilot_file, 'results');
refinement_data = load(refinement_file, 'results');
pilot = pilot_data.results;
refinement = refinement_data.results;
[refs, reference_cfg] = brusselator_three_way_references();
cases = collect_unresolved_cases(pilot, refinement);
if isempty(cases)
    error('run_unresolved_forcing_long_extensions:NoCases', ...
        'No unresolved forcing cases were found.');
end

cfg = brusselator_default_config(struct( ...
    'experiment_name', study_id, 'model', reference_cfg.model, ...
    'grid', reference_cfg.grid, 'initial', reference_cfg.initial, ...
    'solver', pilot.configuration.reachability.discovery_solver));
cfg.reachability = pilot.configuration.reachability;
cfg.long_extension = struct( ...
    'pilot_file', pilot_file, 'refinement_file', refinement_file, ...
    'number_cases', numel(cases), 'target_post_forcing_duration', 1000, ...
    'window_post_forcing_times', [80, 160, 240, 400, 600, 800, 1000], ...
    'recurrence_window', 400, 'spectrum_window', 400, ...
    'recurrence_threshold', cfg.reachability.median_threshold, ...
    'interpretation', ['common-horizon classification and intermediate-dynamics ', ...
        'screen for all previously unresolved forcing trajectories']);

resume = isfolder(outdir) && isfile(progress_file);
if isfolder(outdir) && ~resume
    error('run_unresolved_forcing_long_extensions:OutputExists', ...
        ['Refusing to use an existing output directory without a valid ', ...
        'checkpoint: %s'], outdir);
end
record_template = empty_record();
if resume
    checkpoint = load(progress_file, 'cfg', 'records', ...
        'continuation_files', 'window_tables');
    if numel(checkpoint.records) ~= numel(cases) || ...
            checkpoint.cfg.long_extension.target_post_forcing_duration ~= ...
            cfg.long_extension.target_post_forcing_duration
        error('run_unresolved_forcing_long_extensions:CheckpointMismatch', ...
            'The checkpoint does not match this long-extension study.');
    end
    cfg = checkpoint.cfg;
    records = checkpoint.records;
    continuation_files = checkpoint.continuation_files;
    window_tables = checkpoint.window_tables;
    fprintf('Resuming %s from its saved checkpoint.\n', study_id);
else
    mkdir(outdir);
    mkdir(fullfile(outdir, 'continuations'));
    brusselator_write_metadata(outdir, cfg, ...
        brusselator_run_metadata(cfg, mfilename, NaN));
    records = repmat(record_template, numel(cases), 1);
    continuation_files = cell(numel(cases), 1);
    window_tables = cell(numel(cases), 1);
end

for k = 1:numel(cases)
    item = cases(k);
    if ~strcmp(records(k).final_outcome, 'unclassified') && ...
            ~isempty(continuation_files{k}) && isfile(continuation_files{k})
        fprintf('Skipping completed long extension %d/%d: %s.\n', ...
            k, numel(cases), item.protocol_id);
        continue;
    end
    source_data = load(item.source_raw_file, 'case_result');
    source = source_data.case_result.discovery;
    source_post = source.post_forcing_duration;
    target_post = cfg.long_extension.target_post_forcing_duration;
    if source_post >= target_post
        error('run_unresolved_forcing_long_extensions:UnexpectedSourceHorizon', ...
            'Source %s already reaches the target horizon.', item.protocol_id);
    end

    fprintf('Long extension %d/%d: %s, post %.0f -> %.0f ...\n', ...
        k, numel(cases), item.protocol_id, source_post, target_post);
    case_timer = tic;
    frozen_cfg = source.configuration;
    frozen_cfg.solver = cfg.solver;
    frozen_par = brusselator_make_parameters(frozen_cfg, ...
        @(t) source.configuration.forcing.b0);
    continuation_duration = target_post - source_post;
    [final_state, S2, V2, U2] = solve_brusselator_1d_forced( ...
        source.final_state, frozen_par, continuation_duration, 0);
    S = [source.S; source.S(end) + S2(2:end)];
    U = [source.U; U2(2:end,:)];
    V = [source.V; V2(2:end,:)];

    [final_classification, final_history] = classify(S, U, V, refs, cfg);
    recurrence_start = S(end) - cfg.long_extension.recurrence_window;
    recurrence = brusselator_periodic_recurrence(S, U, V, ...
        recurrence_start, cfg.long_extension.recurrence_threshold);
    spectrum = mean_field_spectrum(S, V, ...
        S(end) - cfg.long_extension.spectrum_window);
    windows = windowed_classification(S, U, V, source.forcing_end, refs, cfg);
    runtime = toc(case_timer);
    record = summarize_case(item, source, final_classification, recurrence, ...
        spectrum, runtime, cfg);
    records(k) = record;
    window_tables{k} = windows;

    continuation_file = fullfile(outdir, 'continuations', ...
        [item.protocol_id, '_continuation.mat']);
    continuation = struct('source_case', item, 'source_post_forcing', source_post, ...
        'target_post_forcing', target_post, 'S', S2, 'U', U2, 'V', V2, ...
        'final_state', final_state, 'final_classification', final_classification, ...
        'final_distance_history', final_history, 'recurrence', recurrence, ...
        'spectrum', spectrum, 'windowed_classification', windows, ...
        'record', record);
    save(continuation_file, 'continuation', '-v7.3');
    continuation_files{k} = continuation_file;
    writetable(struct2table(records(1:k), 'AsArray', true), ...
        fullfile(outdir, 'progress_summary.csv'));
    save(progress_file, 'cfg', 'records', 'continuation_files', ...
        'window_tables', '-v7.3');
    fprintf(['  %s; dynamics=%s, d=[%.3e %.3e %.3e], ', ...
        'period=%.4g, recurrence=%.3e (%.1f s).\n'], ...
        record.final_outcome, record.dynamics_classification, ...
        record.final_median_distance_A, record.final_median_distance_B, ...
        record.final_median_distance_C, record.section_period, ...
        record.first_close_recurrence_distance, runtime);
end

summary_table = struct2table(records, 'AsArray', true);
writetable(summary_table, fullfile(outdir, 'long_extension_summary.csv'));
all_windows = combine_window_tables(cases, window_tables);
writetable(all_windows, fullfile(outdir, 'windowed_distance_summary.csv'));
summary = aggregate_summary(records);
writetable(struct2table(summary, 'AsArray', true), ...
    fullfile(outdir, 'long_extension_counts.csv'));
make_figures(outdir, records, all_windows);
results = struct('configuration', cfg, 'cases', cases, 'records', records, ...
    'summary', summary, 'continuation_files', {continuation_files}, ...
    'window_tables', {window_tables});
save(fullfile(outdir, 'unresolved_forcing_long_extensions.mat'), ...
    'results', '-v7.3');
metadata = brusselator_run_metadata(cfg, mfilename, sum([records.runtime_seconds]));
brusselator_write_metadata(outdir, cfg, metadata);
disp(summary_table);
disp(struct2table(summary, 'AsArray', true));
fprintf('Unresolved forcing extensions complete. Results saved in: %s\n', outdir);
end

function cases = collect_unresolved_cases(pilot, refinement)
template = struct('source_study', '', 'protocol_id', '', ...
    'source_raw_file', '', 'bmax', NaN, 'ramp_time', NaN, ...
    'Tup', NaN, 'Thold', NaN, 'Tdown', NaN);
cases = template([]);
sources = {pilot, refinement};
names = {'forcing_to_c_reachability_pilot_v1', 'ramp_transition_refinement_v1'};
for q = 1:numel(sources)
    source = sources{q};
    unresolved = find(strcmp({source.records.outcome}, 'unresolved_within_window'));
    for j = 1:numel(unresolved)
        index = unresolved(j);
        record = source.records(index);
        item = template;
        item.source_study = names{q};
        item.protocol_id = record.protocol_id;
        item.source_raw_file = source.raw_files{index};
        item.bmax = record.bmax;
        if isfield(record, 'ramp_time')
            item.ramp_time = record.ramp_time;
            item.Tup = record.ramp_time;
            item.Thold = record.Thold;
            item.Tdown = record.ramp_time;
        else
            item.Tup = record.Tup;
            item.Thold = record.Thold;
            item.Tdown = record.Tdown;
            if record.Tup == record.Tdown
                item.ramp_time = record.Tup;
            end
        end
        cases(end+1,1) = item; %#ok<AGROW>
    end
end
end

function [classification, history] = classify(S, U, V, refs, cfg)
[classification, history] = brusselator_classify_three_way_trajectory( ...
    S, U, V, refs, 'late_window', cfg.reachability.late_window, ...
    'sample_interval', cfg.reachability.distance_sample_interval, ...
    'median_threshold', cfg.reachability.median_threshold, ...
    'maximum_threshold', cfg.reachability.maximum_threshold);
end

function windows = windowed_classification(S, U, V, forcing_end, refs, cfg)
post_times = cfg.long_extension.window_post_forcing_times(:);
template = struct('post_forcing_time', NaN, 'median_distance_A', NaN, ...
    'maximum_distance_A', NaN, 'median_distance_B', NaN, ...
    'maximum_distance_B', NaN, 'median_distance_C', NaN, ...
    'maximum_distance_C', NaN, 'outcome', 'unclassified');
records = repmat(template, numel(post_times), 1);
for k = 1:numel(post_times)
    target = forcing_end + post_times(k);
    available = S <= target + 1e-9;
    [classification, ~] = classify(S(available), U(available,:), ...
        V(available,:), refs, cfg);
    record = template;
    record.post_forcing_time = post_times(k);
    record.median_distance_A = classification.late_median_distance_A;
    record.maximum_distance_A = classification.late_maximum_distance_A;
    record.median_distance_B = classification.late_median_distance_B;
    record.maximum_distance_B = classification.late_maximum_distance_B;
    record.median_distance_C = classification.late_median_distance_C;
    record.maximum_distance_C = classification.late_maximum_distance_C;
    record.outcome = classification.outcome;
    records(k) = record;
end
windows = struct2table(records, 'AsArray', true);
end

function spectrum = mean_field_spectrum(S, V, start_time)
late = S >= start_time;
t = S(late);
signal = mean(V(late,:), 2);
signal = signal - mean(signal);
dt = median(diff(t));
n = numel(signal);
power = abs(fft(signal)).^2;
frequency = (0:n-1)' / (n * dt);
positive = 2:floor(n/2);
if isempty(positive) || sum(power(positive)) <= eps
    dominant_frequency = NaN;
    concentration = NaN;
else
    [peak, local_index] = max(power(positive));
    dominant_frequency = frequency(positive(local_index));
    concentration = peak / sum(power(positive));
end
spectrum = struct('window_start', t(1), 'window_end', t(end), ...
    'dominant_frequency', dominant_frequency, ...
    'dominant_period', 1 / dominant_frequency, ...
    'peak_power_fraction', concentration, ...
    'relative_temporal_variation', std(signal) / ...
        max(abs(mean(mean(V(late,:),2))), eps));
end

function record = empty_record()
record = struct('source_study', '', 'protocol_id', '', 'bmax', NaN, ...
    'ramp_time', NaN, 'source_post_forcing_duration', NaN, ...
    'target_post_forcing_duration', NaN, 'final_median_distance_A', NaN, ...
    'final_maximum_distance_A', NaN, 'final_median_distance_B', NaN, ...
    'final_maximum_distance_B', NaN, 'final_median_distance_C', NaN, ...
    'final_maximum_distance_C', NaN, 'final_outcome', 'unclassified', ...
    'number_section_crossings', NaN, 'section_period', NaN, ...
    'section_period_cv', NaN, 'first_close_return_crossings', NaN, ...
    'first_close_recurrence_distance', NaN, 'best_recurrence_distance', NaN, ...
    'dominant_mean_field_period', NaN, 'spectral_peak_power_fraction', NaN, ...
    'relative_temporal_variation', NaN, 'dynamics_classification', '', ...
    'runtime_seconds', NaN);
end

function record = summarize_case(item, source, classification, recurrence, ...
        spectrum, runtime, cfg)
record = empty_record();
record.source_study = item.source_study;
record.protocol_id = item.protocol_id;
record.bmax = item.bmax;
record.ramp_time = item.ramp_time;
record.source_post_forcing_duration = source.post_forcing_duration;
record.target_post_forcing_duration = ...
    cfg.long_extension.target_post_forcing_duration;
record.final_median_distance_A = classification.late_median_distance_A;
record.final_maximum_distance_A = classification.late_maximum_distance_A;
record.final_median_distance_B = classification.late_median_distance_B;
record.final_maximum_distance_B = classification.late_maximum_distance_B;
record.final_median_distance_C = classification.late_median_distance_C;
record.final_maximum_distance_C = classification.late_maximum_distance_C;
record.final_outcome = classification.outcome;
record.number_section_crossings = recurrence.section.crossing_count;
record.section_period = recurrence.section.period;
record.section_period_cv = recurrence.section.period_cv;
record.first_close_return_crossings = ...
    recurrence.recurrence.first_close.return_crossings;
record.first_close_recurrence_distance = ...
    recurrence.recurrence.first_close.best_median_rel_diff_full_state;
record.best_recurrence_distance = ...
    recurrence.recurrence.best.best_median_rel_diff_full_state;
record.dominant_mean_field_period = spectrum.dominant_period;
record.spectral_peak_power_fraction = spectrum.peak_power_fraction;
record.relative_temporal_variation = spectrum.relative_temporal_variation;
record.dynamics_classification = dynamics_label(record, cfg);
record.runtime_seconds = runtime;
end

function label = dynamics_label(record, cfg)
if strcmp(record.final_outcome, 'A_stationary_neighborhood')
    label = 'eventually_A';
elseif strcmp(record.final_outcome, 'B_periodic_reflection_neighborhood')
    label = 'eventually_B';
elseif strcmp(record.final_outcome, 'C_periodic_direct_neighborhood')
    label = 'eventually_C';
elseif isfinite(record.first_close_return_crossings) && ...
        record.first_close_recurrence_distance < ...
        cfg.long_extension.recurrence_threshold && record.section_period_cv < 0.02
    label = 'unresolved_recurrent_intermediate_candidate';
elseif isfinite(record.section_period_cv) && record.section_period_cv < 0.02
    label = 'unresolved_with_regular_section_timing';
else
    label = 'unresolved_nonrecurrent_within_1000';
end
end

function table_out = combine_window_tables(cases, tables)
table_out = table();
for k = 1:numel(tables)
    item = tables{k};
    item.source_study = repmat({cases(k).source_study}, height(item), 1);
    item.protocol_id = repmat({cases(k).protocol_id}, height(item), 1);
    table_out = [table_out; item]; %#ok<AGROW>
end
table_out = movevars(table_out, {'source_study','protocol_id'}, 'Before', 1);
end

function summary = aggregate_summary(records)
summary = struct();
summary.number_cases = numel(records);
summary.eventually_A = sum(strcmp({records.dynamics_classification}, 'eventually_A'));
summary.eventually_B = sum(strcmp({records.dynamics_classification}, 'eventually_B'));
summary.eventually_C = sum(strcmp({records.dynamics_classification}, 'eventually_C'));
summary.unresolved_recurrent = sum(strcmp({records.dynamics_classification}, ...
    'unresolved_recurrent_intermediate_candidate'));
summary.unresolved_regular_timing = sum(strcmp({records.dynamics_classification}, ...
    'unresolved_with_regular_section_timing'));
summary.unresolved_nonrecurrent = sum(strcmp({records.dynamics_classification}, ...
    'unresolved_nonrecurrent_within_1000'));
summary.interpretation = ['finite-time common-horizon resolution at 1000 ', ...
    'post-forcing units; recurrent unresolved cases remain candidates, ', ...
    'not established invariant attractors'];
end

function make_figures(outdir, records, windows)
names = unique(windows.protocol_id, 'stable');
fig = figure('Color', 'w', 'Position', [100 100 1100 720]);
tiledlayout(3,1,'TileSpacing','compact');
fields = {'median_distance_A','median_distance_B','median_distance_C'};
labels = {'d_A','d_B','d_C'};
for q = 1:3
    nexttile; hold on;
    for k = 1:numel(names)
        rows = strcmp(windows.protocol_id, names{k});
        semilogy(windows.post_forcing_time(rows), ...
            max(windows.(fields{q})(rows), eps), '-o', ...
            'DisplayName', names{k});
    end
    ylabel(labels{q}); grid on; box on;
    if q == 1
        title('Windowed distances for previously unresolved forcing cases');
        legend('Location','eastoutside','FontSize',7);
    end
end
xlabel('Post-forcing time');
exportgraphics(fig, fullfile(outdir, 'unresolved_windowed_distances.png'), ...
    'Resolution', 300);
close(fig);

fig = figure('Color', 'w', 'Position', [100 100 1050 460]);
bar([records.best_recurrence_distance]);
set(gca, 'YScale', 'log', 'XTick', 1:numel(records), ...
    'XTickLabel', {records.protocol_id}, 'XTickLabelRotation', 40);
ylabel('Best late Poincare recurrence distance');
title('Late recurrence of extended unresolved cases'); grid on; box on;
exportgraphics(fig, fullfile(outdir, 'unresolved_recurrence_comparison.png'), ...
    'Resolution', 300);
close(fig);
end
