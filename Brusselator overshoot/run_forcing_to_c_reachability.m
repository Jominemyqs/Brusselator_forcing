function results = run_forcing_to_c_reachability()
%RUN_FORCING_TO_C_REACHABILITY Test whether overshoots can prepare outcome C.
%   Samples a focused set of piecewise-linear overshoot protocols at the
%   validated N=400, L=40 problem. For every protocol, both fields are saved
%   from t=0 through at least 80 frozen post-forcing time units. Unresolved
%   cases are adaptively continued to 240 post-forcing units. Outcomes use
%   the full-state phase/reflection-aware A/B/C/U classifier.
%
%   Discovery runs use the moderate ode15s configuration under which C was
%   independently shown to persist. Any apparent C outcome is automatically
%   repeated with the tight reference solver before it is counted as a
%   confirmed forcing preparation.

study_id = 'forcing_to_c_reachability_pilot_v1';
gate_file = fullfile('experiment_outputs', ...
    'candidate_c_numerical_validation_v1', ...
    'candidate_c_numerical_validation.mat');
outdir = fullfile('experiment_outputs', study_id);
progress_file = fullfile(outdir, 'progress_reachability.mat');
if ~isfile(gate_file)
    error('run_forcing_to_c_reachability:MissingValidationGate', ...
        'Expected completed candidate-C validation at %s.', gate_file);
end
gate_data = load(gate_file, 'results');
if ~gate_data.results.summary.proceed_to_L40_forcing_scan
    error('run_forcing_to_c_reachability:ValidationGateFailed', ...
        'The fixed-L=40 candidate-C numerical gate did not pass.');
end

[refs, reference_cfg] = brusselator_three_way_references();
protocols = forcing_protocols();
moderate_solver = struct('name', 'ode15s', 'output_dt', 0.2, ...
    'RelTol', 1e-6, 'AbsTol', 1e-8, 'MaxStep', 0.2);
tight_solver = struct('name', 'ode15s', 'output_dt', 0.1, ...
    'RelTol', 1e-8, 'AbsTol', 1e-10, 'MaxStep', 0.1);
cfg = brusselator_default_config(struct( ...
    'experiment_name', study_id, ...
    'model', reference_cfg.model, ...
    'grid', reference_cfg.grid, ...
    'initial', reference_cfg.initial, ...
    'solver', moderate_solver));
cfg.reachability = struct( ...
    'validation_gate_file', gate_file, ...
    'reference_files', refs.info, ...
    'protocol_count', numel(protocols), ...
    'initial_post_forcing_duration', 80, ...
    'extended_post_forcing_duration', 240, ...
    'late_window', 20, ...
    'distance_sample_interval', 1, ...
    'median_threshold', reference_cfg.classification.recovery_threshold, ...
    'maximum_threshold', reference_cfg.basin_slice.late_max_threshold, ...
    'discovery_solver', moderate_solver, ...
    'C_confirmation_solver', tight_solver, ...
    'outcomes', {{'A', 'B', 'C', 'unresolved'}}, ...
    'interpretation', ['finite-time search for independent forcing preparation ', ...
        'of C at N=400 and L=40']);

resume = isfolder(outdir) && isfile(progress_file);
if isfolder(outdir) && ~resume
    error('run_forcing_to_c_reachability:OutputExists', ...
        ['Refusing to use an existing output directory without a valid ', ...
        'checkpoint: %s'], outdir);
end
record_template = empty_record();
if resume
    checkpoint = load(progress_file, 'cfg', 'records', 'raw_files');
    if numel(checkpoint.records) ~= numel(protocols) || ...
            ~strcmp(checkpoint.cfg.reachability.validation_gate_file, gate_file)
        error('run_forcing_to_c_reachability:CheckpointMismatch', ...
            'The existing checkpoint does not match this protocol matrix.');
    end
    cfg = checkpoint.cfg;
    records = checkpoint.records;
    raw_files = checkpoint.raw_files;
    fprintf('Resuming %s from its saved checkpoint.\n', study_id);
else
    mkdir(outdir);
    mkdir(fullfile(outdir, 'raw_cases'));
    brusselator_write_metadata(outdir, cfg, ...
        brusselator_run_metadata(cfg, mfilename, NaN));
    records = repmat(record_template, numel(protocols), 1);
    raw_files = cell(numel(protocols), 1);
end

[~, common_initial_state] = brusselator_initial_condition(cfg);
for k = 1:numel(protocols)
    protocol = protocols(k);
    if ~strcmp(records(k).outcome, 'unclassified') && ...
            ~isempty(raw_files{k}) && isfile(raw_files{k})
        fprintf('Skipping completed forcing protocol %d/%d: %s.\n', ...
            k, numel(protocols), protocol.id);
        continue;
    end

    fprintf(['Forcing reachability %d/%d: %s, bmax=%.2f, ', ...
        'Tup/hold/down=%.0f/%.0f/%.0f ...\n'], k, numel(protocols), ...
        protocol.id, protocol.bmax, protocol.Tup, protocol.Thold, protocol.Tdown);
    case_timer = tic;
    discovery = execute_protocol(common_initial_state, protocol, cfg, ...
        moderate_solver, refs);
    tight_confirmation = struct();
    if strcmp(discovery.classification.outcome, ...
            'C_periodic_direct_neighborhood')
        fprintf('  Candidate C hit; repeating with the tight solver ...\n');
        tight_confirmation = execute_protocol(common_initial_state, protocol, ...
            cfg, tight_solver, refs);
    end
    runtime = toc(case_timer);
    record = summarize_protocol(protocol, discovery, tight_confirmation, runtime);
    records(k) = record;

    raw_file = fullfile(outdir, 'raw_cases', [protocol.id, '.mat']);
    case_result = struct('protocol', protocol, 'discovery', discovery, ...
        'tight_confirmation', tight_confirmation, 'record', record);
    save(raw_file, 'case_result', '-v7.3');
    raw_files{k} = raw_file;
    writetable(struct2table(records(1:k), 'AsArray', true), ...
        fullfile(outdir, 'progress_summary.csv'));
    save(progress_file, 'cfg', 'records', 'raw_files', '-v7.3');
    fprintf('  %s; landing d=[%.3f %.3f %.3f], post=%.0f (%.1f s).\n', ...
        record.outcome, record.landing_distance_A, record.landing_distance_B, ...
        record.landing_distance_C, record.post_forcing_duration, runtime);
end

summary_table = struct2table(records, 'AsArray', true);
writetable(summary_table, fullfile(outdir, 'forcing_to_c_reachability_summary.csv'));
summary = aggregate_summary(records);
writetable(struct2table(summary, 'AsArray', true), ...
    fullfile(outdir, 'forcing_to_c_reachability_gate.csv'));
make_figures(outdir, records);
results = struct('configuration', cfg, 'protocols', protocols, ...
    'records', records, 'summary', summary, 'raw_files', {raw_files});
save(fullfile(outdir, 'forcing_to_c_reachability.mat'), 'results', '-v7.3');
metadata = brusselator_run_metadata(cfg, mfilename, sum([records.runtime_seconds]));
brusselator_write_metadata(outdir, cfg, metadata);
disp(summary_table);
disp(struct2table(summary, 'AsArray', true));
fprintf('Forcing-to-C reachability pilot complete. Results saved in: %s\n', outdir);
end

function protocols = forcing_protocols()
template = struct('id', '', 'family', '', 'bmax', NaN, 'Tup', NaN, ...
    'Thold', NaN, 'Tdown', NaN);
protocols = template([]);
amplitudes = [11.00, 11.12, 11.14, 11.20, 11.40, 11.60, 12.00, 12.50, 13.00];
for value = amplitudes
    protocols(end+1) = make_protocol(template, ...
        sprintf('amplitude_bmax_%05.2f', value), 'amplitude', value, 40, 80, 40); %#ok<AGROW>
end
for hold = [0, 20, 40, 160, 320]
    protocols(end+1) = make_protocol(template, ...
        sprintf('hold_%03d', hold), 'hold', 11.14, 40, hold, 40); %#ok<AGROW>
end
for ramp = [10, 20, 80, 160]
    protocols(end+1) = make_protocol(template, ...
        sprintf('symmetric_ramp_%03d', ramp), 'symmetric_ramp', ...
        11.14, ramp, 80, ramp); %#ok<AGROW>
end
protocols(end+1) = make_protocol(template, 'fast_up_slow_down', ...
    'asymmetric_ramp', 11.14, 20, 80, 80);
protocols(end+1) = make_protocol(template, 'slow_up_fast_down', ...
    'asymmetric_ramp', 11.14, 80, 80, 20);
end

function protocol = make_protocol(template, id, family, bmax, Tup, Thold, Tdown)
protocol = template;
protocol.id = id;
protocol.family = family;
protocol.bmax = bmax;
protocol.Tup = Tup;
protocol.Thold = Thold;
protocol.Tdown = Tdown;
end

function result = execute_protocol(initial_state, protocol, base_cfg, solver, refs)
initial_post = base_cfg.reachability.initial_post_forcing_duration;
extended_post = base_cfg.reachability.extended_post_forcing_duration;
forcing_end = protocol.Tup + protocol.Thold + protocol.Tdown;
case_cfg = brusselator_default_config(struct( ...
    'experiment_name', [base_cfg.experiment_name, '_', protocol.id], ...
    'model', base_cfg.model, 'grid', base_cfg.grid, 'initial', base_cfg.initial, ...
    'forcing', struct('b0', base_cfg.forcing.b0, 'bmax', protocol.bmax, ...
        'Tup', protocol.Tup, 'Thold', protocol.Thold, 'Tdown', protocol.Tdown), ...
    'time', struct('post_forcing', initial_post), 'solver', solver));
forced_par = brusselator_make_parameters(case_cfg, @(t) overshoot_B(t, ...
    case_cfg.forcing.b0, protocol.bmax, protocol.Tup, protocol.Thold, protocol.Tdown));
[final_state, S, V, U] = solve_brusselator_1d_forced( ...
    initial_state, forced_par, forcing_end + initial_post, 0);

landing_index = find(abs(S - forcing_end) <= max(1e-10, eps(forcing_end)), 1);
if isempty(landing_index)
    [~, landing_index] = min(abs(S - forcing_end));
end
landing_distances = brusselator_three_way_distances( ...
    U(landing_index,:), V(landing_index,:), refs);
[classification, history] = classify(S, U, V, refs, base_cfg);
post_duration = initial_post;
if strcmp(classification.outcome, 'unresolved_within_window')
    frozen_cfg = case_cfg;
    frozen_cfg.solver = solver;
    frozen_par = brusselator_make_parameters(frozen_cfg, @(t) case_cfg.forcing.b0);
    continuation_duration = extended_post - initial_post;
    [final_state, S2, V2, U2] = solve_brusselator_1d_forced( ...
        final_state, frozen_par, continuation_duration, 0);
    S = [S; S(end) + S2(2:end)];
    U = [U; U2(2:end,:)];
    V = [V; V2(2:end,:)];
    post_duration = extended_post;
    [classification, history] = classify(S, U, V, refs, base_cfg);
end

recurrence = struct();
if ismember(classification.outcome, ...
        {'B_periodic_reflection_neighborhood', 'C_periodic_direct_neighborhood'})
    recurrence = brusselator_periodic_recurrence(S, U, V, ...
        S(end) - 40, base_cfg.reachability.median_threshold);
end
result = struct('configuration', case_cfg, 'forcing_end', forcing_end, ...
    'post_forcing_duration', post_duration, 'landing_time', S(landing_index), ...
    'landing_state', [linspace(0, case_cfg.grid.Lx, case_cfg.grid.N)', ...
        U(landing_index,:)', V(landing_index,:)'], ...
    'landing_distances', landing_distances, 'S', S, 'U', U, 'V', V, ...
    'final_state', final_state, 'classification', classification, ...
    'distance_history', history, 'recurrence', recurrence);
end

function [classification, history] = classify(S, U, V, refs, cfg)
[classification, history] = brusselator_classify_three_way_trajectory( ...
    S, U, V, refs, ...
    'late_window', cfg.reachability.late_window, ...
    'sample_interval', cfg.reachability.distance_sample_interval, ...
    'median_threshold', cfg.reachability.median_threshold, ...
    'maximum_threshold', cfg.reachability.maximum_threshold);
end

function record = empty_record()
record = struct('protocol_id', '', 'family', '', 'bmax', NaN, 'Tup', NaN, ...
    'Thold', NaN, 'Tdown', NaN, 'forcing_end', NaN, ...
    'post_forcing_duration', NaN, 'landing_distance_A', NaN, ...
    'landing_distance_B', NaN, 'landing_distance_C', NaN, ...
    'late_median_distance_A', NaN, 'late_median_distance_B', NaN, ...
    'late_median_distance_C', NaN, 'outcome', 'unclassified', ...
    'period_if_recurrent', NaN, 'recurrence_distance', NaN, ...
    'tight_confirmation_outcome', 'not_run', ...
    'tight_confirmation_distance_C', NaN, 'confirmed_C', false, ...
    'runtime_seconds', NaN);
end

function record = summarize_protocol(protocol, discovery, confirmation, runtime)
record = empty_record();
record.protocol_id = protocol.id;
record.family = protocol.family;
record.bmax = protocol.bmax;
record.Tup = protocol.Tup;
record.Thold = protocol.Thold;
record.Tdown = protocol.Tdown;
record.forcing_end = discovery.forcing_end;
record.post_forcing_duration = discovery.post_forcing_duration;
record.landing_distance_A = discovery.landing_distances(1);
record.landing_distance_B = discovery.landing_distances(2);
record.landing_distance_C = discovery.landing_distances(3);
record.late_median_distance_A = discovery.classification.late_median_distance_A;
record.late_median_distance_B = discovery.classification.late_median_distance_B;
record.late_median_distance_C = discovery.classification.late_median_distance_C;
record.outcome = discovery.classification.outcome;
if ~isempty(fieldnames(discovery.recurrence))
    record.period_if_recurrent = discovery.recurrence.section.period;
    record.recurrence_distance = ...
        discovery.recurrence.recurrence.first_close.best_median_rel_diff_full_state;
end
if ~isempty(fieldnames(confirmation))
    record.tight_confirmation_outcome = confirmation.classification.outcome;
    record.tight_confirmation_distance_C = ...
        confirmation.classification.late_median_distance_C;
    record.confirmed_C = strcmp(record.outcome, 'C_periodic_direct_neighborhood') && ...
        strcmp(record.tight_confirmation_outcome, ...
        'C_periodic_direct_neighborhood');
end
record.runtime_seconds = runtime;
end

function summary = aggregate_summary(records)
summary = struct();
summary.number_protocols = numel(records);
summary.number_A = sum(strcmp({records.outcome}, 'A_stationary_neighborhood'));
summary.number_B = sum(strcmp({records.outcome}, ...
    'B_periodic_reflection_neighborhood'));
summary.number_C_discovery = sum(strcmp({records.outcome}, ...
    'C_periodic_direct_neighborhood'));
summary.number_unresolved = sum(strcmp({records.outcome}, ...
    'unresolved_within_window'));
summary.number_confirmed_C = sum([records.confirmed_C]);
summary.forcing_reaches_C = summary.number_confirmed_C > 0;
if summary.forcing_reaches_C
    summary.interpretation = ['at least one sampled overshoot independently ', ...
        'prepares C under both discovery and tight solver configurations'];
else
    summary.interpretation = ['no tight-confirmed C preparation in this focused ', ...
        'finite protocol sample; this does not prove C is unreachable'];
end
end

function make_figures(outdir, records)
fig = figure('Color', 'w', 'Position', [100 100 1150 560]);
x = 1:numel(records);
bar(x - 0.25, [records.landing_distance_A], 0.25, 'DisplayName', 'A'); hold on;
bar(x, [records.landing_distance_B], 0.25, 'DisplayName', 'B');
bar(x + 0.25, [records.landing_distance_C], 0.25, 'DisplayName', 'C');
set(gca, 'YScale', 'log', 'XTick', x, 'XTickLabel', {records.protocol_id}, ...
    'XTickLabelRotation', 45);
ylabel('Full-state distance at return to b=10');
title('Forcing landing states relative to A, B, and C');
legend('Location', 'eastoutside'); grid on; box on;
exportgraphics(fig, fullfile(outdir, 'forcing_landing_distances.png'), ...
    'Resolution', 300);
close(fig);

fig = figure('Color', 'w', 'Position', [100 100 1150 500]);
bar(x - 0.25, [records.late_median_distance_A], 0.25, 'DisplayName', 'A'); hold on;
bar(x, [records.late_median_distance_B], 0.25, 'DisplayName', 'B');
bar(x + 0.25, [records.late_median_distance_C], 0.25, 'DisplayName', 'C');
set(gca, 'YScale', 'log', 'XTick', x, 'XTickLabel', {records.protocol_id}, ...
    'XTickLabelRotation', 45);
ylabel('Late median phase-aware distance');
title('Frozen post-forcing outcomes');
legend('Location', 'eastoutside'); grid on; box on;
exportgraphics(fig, fullfile(outdir, 'forcing_outcome_distances.png'), ...
    'Resolution', 300);
close(fig);
end
