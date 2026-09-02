function results = run_recurrent_class_local_stability()
%RUN_RECURRENT_CLASS_LOCAL_STABILITY Probe attraction near R5 and R6.
%   Applies smooth Neumann-compatible perturbations at one standardized phase
%   of each tight N=400 representative. Two symmetry sectors and two relative
%   amplitudes are tested. Distances minimize over temporal phase and exact
%   reflection. This is a reduced finite-direction local screen, not a Floquet
%   calculation or proof of nonlinear stability.

study_id = 'recurrent_class_local_stability_v1';
parent_file = fullfile('experiment_outputs', ...
    'recurrent_class_validation_v1', 'recurrent_class_validation.mat');
long_file = fullfile('experiment_outputs', ...
    'unresolved_forcing_long_extension_v1', ...
    'unresolved_forcing_long_extensions.mat');
outdir = fullfile('experiment_outputs', study_id);
progress_file = fullfile(outdir, 'progress_recurrent_class_local_stability.mat');
if ~isfile(parent_file) || ~isfile(long_file)
    error('run_recurrent_class_local_stability:MissingParent', ...
        'The recurrent validation and long-extension studies are required.');
end
loaded = load(parent_file, 'results');
parent = loaded.results;
loaded = load(long_file, 'results');
long_parent = loaded.results;
[references, reference_files] = load_references(parent);
cases = stability_cases();

cfg = parent.configuration;
cfg.experiment_name = study_id;
cfg.solver = struct('name', 'ode15s', 'output_dt', 0.1, ...
    'RelTol', 1e-8, 'AbsTol', 1e-10, 'MaxStep', 0.1);
cfg.recurrent_class_local_stability = struct( ...
    'parent_file', parent_file, 'reference_files', {reference_files}, ...
    'frozen_b', cfg.forcing.b0, 'duration', 60, 'late_window', 15, ...
    'sample_interval', 1, 'phase_fraction', 0, ...
    'amplitudes_relative_full_state', [1e-4, 1e-3], ...
    'directions', {{'reflection_odd_smooth','reflection_even_smooth'}}, ...
    'late_median_threshold', cfg.classification.recovery_threshold, ...
    'late_maximum_threshold', 2e-2, ...
    'interpretation', ['reduced finite-time local-return screen; not a ', ...
        'Floquet or full-basin calculation']);

resume = isfolder(outdir) && isfile(progress_file);
if isfolder(outdir) && ~resume
    error('run_recurrent_class_local_stability:OutputExists', ...
        'Refusing existing directory without a valid checkpoint: %s', outdir);
end
if resume
    checkpoint = load(progress_file, 'cfg', 'records', 'raw_files');
    if numel(checkpoint.records) ~= numel(cases)
        error('run_recurrent_class_local_stability:CheckpointMismatch', ...
            'The checkpoint does not match this stability design.');
    end
    cfg = checkpoint.cfg;
    records = checkpoint.records;
    raw_files = checkpoint.raw_files;
    fprintf('Resuming %s from its saved checkpoint.\n', study_id);
else
    mkdir(outdir);
    mkdir(fullfile(outdir, 'raw_cases'));
    records = repmat(empty_record(), numel(cases), 1);
    raw_files = cell(numel(cases),1);
    brusselator_write_metadata(outdir, cfg, ...
        brusselator_run_metadata(cfg, mfilename, NaN));
end

directions = smooth_neumann_directions(cfg.grid.Lx, cfg.grid.N);
frozen_par = brusselator_make_parameters(cfg, ...
    @(t) cfg.recurrent_class_local_stability.frozen_b);
x = linspace(0, cfg.grid.Lx, cfg.grid.N)';
for k = 1:numel(cases)
    spec = cases(k);
    if ~strcmp(records(k).classification, 'unclassified') && ...
            ~isempty(raw_files{k}) && isfile(raw_files{k})
        fprintf('Skipping completed stability case %d/%d: %s.\n', ...
            k, numel(cases), spec.case_id);
        continue;
    end
    template = references.(spec.class_id);
    direction = directions(strcmp({directions.name}, spec.direction));
    base_u = template.U(1,:);
    base_v = template.V(1,:);
    [u0, v0] = apply_perturbation(base_u, base_v, direction, spec.amplitude);
    initial_state = [x, u0', v0'];
    fprintf('Stability %d/%d: %s, %s, amplitude %.0e ...\n', ...
        k, numel(cases), spec.class_id, spec.direction, spec.amplitude);
    timer = tic;
    [final_state, S, V, U] = solve_brusselator_1d_forced( ...
        initial_state, frozen_par, ...
        cfg.recurrent_class_local_stability.duration, 0);
    runtime = toc(timer);
    stride = max(1, round(cfg.recurrent_class_local_stability.sample_interval / ...
        cfg.solver.output_dt));
    indices = unique([1:stride:numel(S), numel(S)]);
    [distances, phases, reflections] = ...
        brusselator_states_to_orbit_distances(U(indices,:), V(indices,:), template);
    sample_times = S(indices);
    late = sample_times >= cfg.recurrent_class_local_stability.duration - ...
        cfg.recurrent_class_local_stability.late_window;
    record = summarize_case(spec, distances, phases, reflections, ...
        sample_times, late, runtime, cfg);
    records(k) = record;
    case_result = struct('configuration', cfg, 'specification', spec, ...
        'initial_state', initial_state, 'final_state', final_state, ...
        'S', S, 'U', U, 'V', V, 'sample_times', sample_times, ...
        'orbit_distances', distances, 'best_phase_fraction', phases, ...
        'best_reflection', reflections, 'record', record);
    case_file = fullfile(outdir, 'raw_cases', [spec.case_id, '.mat']);
    save(case_file, 'case_result', '-v7.3');
    raw_files{k} = case_file;
    writetable(struct2table(records, 'AsArray', true), ...
        fullfile(outdir, 'progress_summary.csv'));
    save(progress_file, 'cfg', 'records', 'raw_files', '-v7.3');
    fprintf('  %s; late median %.3e, max %.3e (%.1f s).\n', ...
        record.classification, record.late_median_orbit_distance, ...
        record.late_maximum_orbit_distance, runtime);
end

summary_table = struct2table(records, 'AsArray', true);
writetable(summary_table, fullfile(outdir, 'local_stability_summary.csv'));
gates = aggregate_gates(records, parent.gates, long_parent.records);
writetable(struct2table(gates, 'AsArray', true), ...
    fullfile(outdir, 'recurrent_class_final_evidence_gates.csv'));
make_figure(outdir, records, raw_files, cfg);
results = struct('configuration', cfg, 'parent_file', parent_file, ...
    'references', references, 'records', records, 'raw_files', {raw_files}, ...
    'gates', gates);
save(fullfile(outdir, 'recurrent_class_local_stability.mat'), ...
    'results', '-v7.3');
brusselator_write_metadata(outdir, cfg, brusselator_run_metadata( ...
    cfg, mfilename, sum([records.runtime_seconds])));
disp(summary_table);
disp(struct2table(gates, 'AsArray', true));
fprintf('Recurrent-class local-stability screen saved in: %s\n', outdir);
end

function [references, files] = load_references(parent)
references = struct();
classes = {'period_5p0585','period_6p5837'};
files = cell(2,1);
for k = 1:2
    index = find(strcmp({parent.validation_records.class_id}, classes{k}) & ...
        [parent.validation_records.N] == 400 & ...
        strcmp({parent.validation_records.factor}, 'tolerance'), 1);
    if isempty(index)
        error('run_recurrent_class_local_stability:MissingReference', ...
            'No tight N=400 reference exists for %s.', classes{k});
    end
    files{k} = parent.validation_raw_files{index};
    loaded = load(files{k}, 'case_result');
    references.(classes{k}) = loaded.case_result.template;
end
end

function cases = stability_cases()
classes = {'period_5p0585','period_6p5837'};
directions = {'reflection_odd_smooth','reflection_even_smooth'};
amplitudes = [1e-4, 1e-3];
template = struct('case_id', '', 'class_id', '', 'phase_fraction', 0, ...
    'direction', '', 'amplitude', NaN);
cases = repmat(template, numel(classes)*numel(directions)*numel(amplitudes), 1);
count = 0;
for c = 1:numel(classes)
    for d = 1:numel(directions)
        for a = 1:numel(amplitudes)
            count = count + 1;
            cases(count) = struct('case_id', sprintf('%s_%s_%0.0e', ...
                classes{c}, directions{d}, amplitudes(a)), ...
                'class_id', classes{c}, 'phase_fraction', 0, ...
                'direction', directions{d}, 'amplitude', amplitudes(a));
        end
    end
end
end

function directions = smooth_neumann_directions(Lx, N)
x = linspace(0, Lx, N);
directions(1) = struct('name', 'reflection_odd_smooth', ...
    'u', cos(pi*x/Lx) + 0.35*cos(3*pi*x/Lx), ...
    'v', -0.7*cos(pi*x/Lx) + 0.25*cos(3*pi*x/Lx));
directions(2) = struct('name', 'reflection_even_smooth', ...
    'u', cos(2*pi*x/Lx) - 0.4*cos(4*pi*x/Lx), ...
    'v', 0.6*cos(2*pi*x/Lx) + 0.2*cos(4*pi*x/Lx));
end

function [u, v] = apply_perturbation(base_u, base_v, direction, amplitude)
base_norm = sqrt(sum(base_u.^2) + sum(base_v.^2));
direction_norm = sqrt(sum(direction.u.^2) + sum(direction.v.^2));
scale = amplitude * base_norm / max(direction_norm, eps);
u = base_u + scale * direction.u;
v = base_v + scale * direction.v;
end

function record = empty_record()
record = struct('case_id', '', 'class_id', '', 'phase_fraction', NaN, ...
    'direction', '', 'amplitude_relative_full_state', NaN, ...
    'initial_orbit_distance', NaN, 'final_orbit_distance', NaN, ...
    'late_median_orbit_distance', NaN, 'late_maximum_orbit_distance', NaN, ...
    'late_to_initial_ratio', NaN, 'best_final_phase_fraction', NaN, ...
    'best_final_reflection', false, 'runtime_seconds', NaN, ...
    'classification', 'unclassified');
end

function record = summarize_case(spec, distances, phases, reflections, ...
        sample_times, late, runtime, cfg)
record = empty_record();
record.case_id = spec.case_id;
record.class_id = spec.class_id;
record.phase_fraction = spec.phase_fraction;
record.direction = spec.direction;
record.amplitude_relative_full_state = spec.amplitude;
record.initial_orbit_distance = distances(1);
record.final_orbit_distance = distances(end);
record.late_median_orbit_distance = median(distances(late));
record.late_maximum_orbit_distance = max(distances(late));
record.late_to_initial_ratio = record.late_median_orbit_distance / ...
    max(record.initial_orbit_distance, eps);
record.best_final_phase_fraction = phases(end);
record.best_final_reflection = reflections(end);
record.runtime_seconds = runtime;
if record.late_median_orbit_distance < ...
        cfg.recurrent_class_local_stability.late_median_threshold && ...
        record.late_maximum_orbit_distance < ...
        cfg.recurrent_class_local_stability.late_maximum_threshold
    record.classification = 'returns_to_representative_orbit_neighborhood';
elseif record.final_orbit_distance > cfg.classification.transition_threshold
    record.classification = 'does_not_return_within_test_window';
else
    record.classification = 'intermediate_or_inconclusive';
end
if sample_times(end) < cfg.recurrent_class_local_stability.duration
    error('run_recurrent_class_local_stability:IncompleteTrajectory', ...
        'A stability trajectory ended before its requested duration.');
end
end

function gates = aggregate_gates(records, numerical_gates, long_records)
classes = {'period_5p0585','period_6p5837'};
representatives = {'amplitude_bmax_12.50','symmetric_ramp_012'};
gates = repmat(struct('class_id', '', 'periodic_orbit_supported', false, ...
    'number_local_cases', NaN, 'number_local_returns', NaN, ...
    'local_attraction_supported', false, ...
    'cross_preparation_orbit_supported', false, ...
    'separated_from_other_period_class', false, ...
    'separated_from_A_B_C', false, ...
    'finite_time_distinct_attractor_supported', false, ...
    'interpretation', ['finite-time numerical evidence restricted to Lx=40; ', ...
        'not a Floquet proof or continuum theorem']), 2, 1);
for k = 1:2
    local = strcmp({records.class_id}, classes{k});
    returned = strcmp({records(local).classification}, ...
        'returns_to_representative_orbit_neighborhood');
    numerical = numerical_gates(strcmp({numerical_gates.class_id}, classes{k}));
    representative = long_records(strcmp({long_records.protocol_id}, ...
        representatives{k}));
    separated_ABC = min([representative.final_median_distance_A, ...
        representative.final_median_distance_B, ...
        representative.final_median_distance_C]) > 0.1;
    gates(k).class_id = classes{k};
    gates(k).periodic_orbit_supported = numerical.periodic_orbit_supported;
    gates(k).number_local_cases = sum(local);
    gates(k).number_local_returns = sum(returned);
    gates(k).local_attraction_supported = all(returned);
    gates(k).cross_preparation_orbit_supported = ...
        numerical.cross_preparation_orbit_supported;
    gates(k).separated_from_other_period_class = ...
        numerical.separated_from_other_period_class;
    gates(k).separated_from_A_B_C = separated_ABC;
    gates(k).finite_time_distinct_attractor_supported = ...
        gates(k).periodic_orbit_supported && ...
        gates(k).local_attraction_supported && ...
        gates(k).cross_preparation_orbit_supported && ...
        gates(k).separated_from_other_period_class && separated_ABC;
end
end

function make_figure(outdir, records, raw_files, cfg)
fig = figure('Color', 'w', 'Position', [100 100 1000 560]);
hold on;
for k = 1:numel(records)
    loaded = load(raw_files{k}, 'case_result');
    item = loaded.case_result;
    semilogy(item.sample_times, max(item.orbit_distances, eps), ...
        'LineWidth', 1.1, 'DisplayName', sprintf('%s, %s, %.0e', ...
        records(k).class_id, records(k).direction, ...
        records(k).amplitude_relative_full_state));
end
yline(cfg.recurrent_class_local_stability.late_median_threshold, '--', ...
    'median gate', 'HandleVisibility','off');
yline(cfg.recurrent_class_local_stability.late_maximum_threshold, ':', ...
    'maximum gate', 'HandleVisibility','off');
xlabel('Frozen b=10 integration time');
ylabel('Distance to representative orbit modulo phase/reflection');
title('Reduced local-attraction screen for recurrent classes');
legend('Location','eastoutside','FontSize',7); grid on; box on;
exportgraphics(fig, fullfile(outdir, 'local_orbit_distance_by_case.png'), ...
    'Resolution', 300);
close(fig);
end
