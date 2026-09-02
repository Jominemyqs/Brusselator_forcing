function results = run_candidate_c_numerical_validation()
%RUN_CANDIDATE_C_NUMERICAL_VALIDATION Test whether candidate C is robust.
%   Regrids a late state on the documented N=400 candidate-C orbit and
%   evolves it at frozen b=10 under changes in grid, solver, tolerance, and
%   domain length. Each case must show close full-state recurrence, regular
%   Poincare timing, a period near the source period, and a standardized
%   section profile near the source C profile to be labeled as supporting C.
%
%   Domain cases stretch/compress the source profile in normalized x/L.
%   They test persistence after a modest domain deformation, not independent
%   discovery of C on that domain. The later forcing-reachability experiment
%   is the intended independent preparation test.

study_id = 'candidate_c_numerical_validation_v1';
source_file = fullfile('experiment_outputs', ...
    'edge_candidate_frozen_confirmation_v1', ...
    'candidate_confirmation_trajectory.mat');
outdir = fullfile('experiment_outputs', study_id);
if ~isfile(source_file)
    error('run_candidate_c_numerical_validation:MissingSource', ...
        'Expected candidate-C source data at %s.', source_file);
end
progress_file = fullfile(outdir, 'progress_validation.mat');
resume = isfolder(outdir) && isfile(progress_file);
if isfolder(outdir) && ~resume
    error('run_candidate_c_numerical_validation:OutputExists', ...
        ['Refusing to use an existing experiment directory without a valid ', ...
        'checkpoint: %s'], outdir);
end

loaded = load(source_file, 'results');
source = loaded.results;
source_cfg = source.configuration;
source_period = source.section.period;
source_section_state = source.section.states(end,:);
source_initial_state = source.trajectory.final_state;
source_N = source_cfg.grid.N;
if numel(source_section_state) ~= 2 * source_N
    error('run_candidate_c_numerical_validation:SourceDimensionMismatch', ...
        'Candidate-C source section does not match its configured grid.');
end

cases = validation_cases();
cfg = source_cfg;
cfg.experiment_name = study_id;
cfg.candidate_c_validation = struct( ...
    'source_file', source_file, ...
    'source_grid_N', source_N, ...
    'source_domain_Lx', source_cfg.grid.Lx, ...
    'source_period', source_period, ...
    'state_transfer', 'shape-preserving interpolation in normalized coordinate x/L', ...
    'recurrence_threshold', source_cfg.classification.recovery_threshold, ...
    'maximum_period_cv', 0.02, ...
    'maximum_relative_period_error', 0.05, ...
    'maximum_standardized_profile_distance', 0.10, ...
    'interpretation', ['finite-time persistence and numerical-sensitivity ', ...
        'screen; not a Floquet calculation or continuum proof'], ...
    'cases', cases);

record_template = empty_record();
if resume
    checkpoint = load(progress_file, 'cfg', 'records', 'trajectory_files');
    if ~isfield(checkpoint.cfg, 'candidate_c_validation') || ...
            ~strcmp(checkpoint.cfg.candidate_c_validation.source_file, source_file) || ...
            numel(checkpoint.records) ~= numel(cases)
        error('run_candidate_c_numerical_validation:CheckpointMismatch', ...
            'The existing checkpoint does not match this validation matrix.');
    end
    cfg = checkpoint.cfg;
    records = checkpoint.records;
    trajectory_files = checkpoint.trajectory_files;
    fprintf('Resuming %s from its saved checkpoint.\n', study_id);
else
    mkdir(outdir);
    mkdir(fullfile(outdir, 'raw_cases'));
    brusselator_write_metadata(outdir, cfg, ...
        brusselator_run_metadata(cfg, mfilename, NaN));
    records = repmat(record_template, numel(cases), 1);
    trajectory_files = cell(numel(cases), 1);
end
for k = 1:numel(cases)
    spec = cases(k);
    if ~strcmp(records(k).classification, 'unclassified') && ...
            ~isempty(trajectory_files{k}) && isfile(trajectory_files{k})
        fprintf('Skipping completed validation case %d/%d: %s.\n', ...
            k, numel(cases), spec.id);
        continue;
    end
    case_cfg = brusselator_default_config(struct( ...
        'experiment_name', [study_id, '_', spec.id], ...
        'model', source_cfg.model, ...
        'grid', struct('N', spec.N, 'Lx', spec.Lx), ...
        'forcing', source_cfg.forcing, ...
        'initial', source_cfg.initial, ...
        'solver', spec.solver));
    initial_state = transfer_state(source_initial_state, ...
        source_cfg.grid.Lx, spec.N, spec.Lx);
    frozen_par = brusselator_make_parameters(case_cfg, ...
        @(t) source_cfg.candidate_confirmation.frozen_b);

    fprintf(['Candidate-C validation %d/%d: %s, N=%d, L=%.1f, ', ...
        '%s, T=%.0f ...\n'], k, numel(cases), spec.id, spec.N, ...
        spec.Lx, spec.solver.name, spec.duration);
    case_timer = tic;
    [final_state, S, V, U] = solve_brusselator_1d_forced( ...
        initial_state, frozen_par, spec.duration, 0);
    runtime = toc(case_timer);

    recurrence = brusselator_periodic_recurrence(S, U, V, ...
        spec.analysis_start, cfg.candidate_c_validation.recurrence_threshold);
    record = summarize_case(spec, recurrence, source_section_state, ...
        source_cfg.grid.Lx, source_period, cfg, runtime);
    records(k) = record;

    case_file = fullfile(outdir, 'raw_cases', [spec.id, '.mat']);
    case_result = struct('configuration', case_cfg, 'specification', spec, ...
        'initial_state', initial_state, 'final_state', final_state, ...
        'S', S, 'U', U, 'V', V, 'recurrence', recurrence, 'record', record);
    save(case_file, 'case_result', '-v7.3');
    trajectory_files{k} = case_file;
    writetable(struct2table(records(1:k), 'AsArray', true), ...
        fullfile(outdir, 'progress_summary.csv'));
    save(progress_file, 'cfg', 'records', ...
        'trajectory_files', '-v7.3');
    fprintf('  %s; period %.5f, recurrence %.3e, profile %.3e (%.1f s).\n', ...
        record.classification, record.section_period, ...
        record.first_close_recurrence_distance, ...
        record.standardized_profile_distance_to_source_C, runtime);
end

summary_table = struct2table(records, 'AsArray', true);
writetable(summary_table, fullfile(outdir, 'candidate_c_numerical_validation_summary.csv'));
make_figures(outdir, records, cases, cfg);
summary = aggregate_summary(records);
writetable(struct2table(summary, 'AsArray', true), ...
    fullfile(outdir, 'candidate_c_validation_gate.csv'));

results = struct('configuration', cfg, 'source_file', source_file, ...
    'records', records, 'summary', summary, ...
    'trajectory_files', {trajectory_files});
save(fullfile(outdir, 'candidate_c_numerical_validation.mat'), ...
    'results', '-v7.3');
metadata = brusselator_run_metadata(cfg, mfilename, sum([records.runtime_seconds]));
brusselator_write_metadata(outdir, cfg, metadata);
disp(summary_table);
disp(struct2table(summary, 'AsArray', true));
fprintf('Candidate-C numerical validation complete. Results saved in: %s\n', outdir);
end

function cases = validation_cases()
tight = struct('name', 'ode15s', 'output_dt', 0.1, ...
    'RelTol', 1e-8, 'AbsTol', 1e-10, 'MaxStep', 0.1);
moderate = struct('name', 'ode15s', 'output_dt', 0.1, ...
    'RelTol', 1e-6, 'AbsTol', 1e-8, 'MaxStep', 0.2);
explicit = struct('name', 'ode45', 'output_dt', 0.1, ...
    'RelTol', 1e-7, 'AbsTol', 1e-9, 'MaxStep', 0.05);
template = struct('id', '', 'factor', '', 'N', NaN, 'Lx', NaN, ...
    'solver', tight, 'duration', 80, 'analysis_start', 40);
cases = repmat(template, 7, 1);
cases(1) = set_case(template, 'grid_N200_tight', 'grid', 200, 40, tight);
cases(2) = set_case(template, 'grid_N600_tight', 'grid', 600, 40, tight);
cases(3) = set_case(template, 'grid_N800_tight', 'grid', 800, 40, tight);
cases(4) = set_case(template, 'solver_ode45_N200', 'solver', 200, 40, explicit);
cases(5) = set_case(template, 'tolerance_N400_moderate', ...
    'tolerance', 400, 40, moderate);
cases(6) = set_case(template, 'domain_L38_dxmatched', 'domain', 380, 38, tight);
cases(7) = set_case(template, 'domain_L42_dxmatched', 'domain', 420, 42, tight);
end

function spec = set_case(template, id, factor, N, Lx, solver)
spec = template;
spec.id = id;
spec.factor = factor;
spec.N = N;
spec.Lx = Lx;
spec.solver = solver;
end

function transferred = transfer_state(source_state, source_Lx, target_N, target_Lx)
source_coordinate = source_state(:,1) / source_Lx;
target_x = linspace(0, target_Lx, target_N)';
target_coordinate = target_x / target_Lx;
u = interp1(source_coordinate, source_state(:,2), target_coordinate, 'pchip');
v = interp1(source_coordinate, source_state(:,3), target_coordinate, 'pchip');
transferred = [target_x, u, v];
end

function record = empty_record()
record = struct('case_id', '', 'factor', '', 'N', NaN, 'Lx', NaN, ...
    'dx', NaN, 'solver_name', '', 'RelTol', NaN, 'AbsTol', NaN, ...
    'MaxStep', NaN, 'duration', NaN, 'analysis_start', NaN, ...
    'number_section_crossings', NaN, 'section_period', NaN, ...
    'section_period_cv', NaN, 'relative_period_error_from_source_C', NaN, ...
    'first_close_return_crossings', NaN, ...
    'first_close_recurrence_distance', NaN, ...
    'first_close_uses_reflection', false, ...
    'best_recurrence_distance', NaN, ...
    'standardized_profile_distance_to_source_C', NaN, ...
    'runtime_seconds', NaN, 'classification', 'unclassified');
end

function record = summarize_case(spec, analysis, source_section_state, ...
        source_Lx, source_period, cfg, runtime)
record = empty_record();
record.case_id = spec.id;
record.factor = spec.factor;
record.N = spec.N;
record.Lx = spec.Lx;
record.dx = spec.Lx / (spec.N - 1);
record.solver_name = spec.solver.name;
record.RelTol = spec.solver.RelTol;
record.AbsTol = spec.solver.AbsTol;
record.MaxStep = spec.solver.MaxStep;
record.duration = spec.duration;
record.analysis_start = spec.analysis_start;
record.number_section_crossings = analysis.section.crossing_count;
record.section_period = analysis.section.period;
record.section_period_cv = analysis.section.period_cv;
record.relative_period_error_from_source_C = ...
    abs(analysis.section.period - source_period) / source_period;
record.first_close_return_crossings = ...
    analysis.recurrence.first_close.return_crossings;
record.first_close_recurrence_distance = ...
    analysis.recurrence.first_close.best_median_rel_diff_full_state;
record.first_close_uses_reflection = ...
    analysis.recurrence.first_close.best_uses_reflection;
record.best_recurrence_distance = ...
    analysis.recurrence.best.best_median_rel_diff_full_state;
if analysis.section.crossing_count > 0
    record.standardized_profile_distance_to_source_C = section_profile_distance( ...
        analysis.section.states(end,:), spec.Lx, source_section_state, source_Lx);
end
record.runtime_seconds = runtime;

has_close_return = isfinite(record.first_close_return_crossings) && ...
    record.first_close_recurrence_distance < ...
    cfg.candidate_c_validation.recurrence_threshold;
regular_timing = record.section_period_cv < ...
    cfg.candidate_c_validation.maximum_period_cv;
matching_period = record.relative_period_error_from_source_C < ...
    cfg.candidate_c_validation.maximum_relative_period_error;
matching_profile = record.standardized_profile_distance_to_source_C < ...
    cfg.candidate_c_validation.maximum_standardized_profile_distance;
if has_close_return && regular_timing && matching_period && matching_profile
    record.classification = 'supports_candidate_C_persistence';
elseif has_close_return && regular_timing
    record.classification = 'recurrent_periodic_but_not_identified_as_C';
else
    record.classification = 'does_not_validate_C_within_test_window';
end
end

function distance = section_profile_distance(candidate, candidate_Lx, source, source_Lx)
candidate_N = numel(candidate) / 2;
source_N = numel(source) / 2;
common = linspace(0, 1, 801);
candidate_coordinate = linspace(0, candidate_Lx, candidate_N) / candidate_Lx;
source_coordinate = linspace(0, source_Lx, source_N) / source_Lx;
candidate_common = [interp1(candidate_coordinate, candidate(1:candidate_N), ...
    common, 'pchip'), interp1(candidate_coordinate, candidate(candidate_N+1:end), ...
    common, 'pchip')];
source_common = [interp1(source_coordinate, source(1:source_N), common, 'pchip'), ...
    interp1(source_coordinate, source(source_N+1:end), common, 'pchip')];
source_reflected = [fliplr(source_common(1:801)), ...
    fliplr(source_common(802:end))];
direct = norm(candidate_common - source_common) / max(norm(source_common), eps);
reflected = norm(candidate_common - source_reflected) / ...
    max(norm(source_reflected), eps);
distance = min(direct, reflected);
end

function summary = aggregate_summary(records)
summary = struct();
summary.number_cases = numel(records);
summary.number_supporting_C = sum(strcmp({records.classification}, ...
    'supports_candidate_C_persistence'));
summary.grid_cases_pass = all(strcmp({records(strcmp({records.factor}, 'grid')).classification}, ...
    'supports_candidate_C_persistence'));
summary.solver_case_pass = all(strcmp({records(strcmp({records.factor}, 'solver')).classification}, ...
    'supports_candidate_C_persistence'));
summary.tolerance_case_pass = all(strcmp({records(strcmp({records.factor}, 'tolerance')).classification}, ...
    'supports_candidate_C_persistence'));
summary.domain_cases_pass = all(strcmp({records(strcmp({records.factor}, 'domain')).classification}, ...
    'supports_candidate_C_persistence'));
summary.fixed_domain_numerical_gate_pass = summary.grid_cases_pass && ...
    summary.solver_case_pass && summary.tolerance_case_pass;
summary.proceed_to_L40_forcing_scan = summary.fixed_domain_numerical_gate_pass;
if summary.fixed_domain_numerical_gate_pass && summary.domain_cases_pass
    summary.interpretation = ['candidate C passes the fixed-domain numerical ', ...
        'screen and the sampled domain deformations; proceed to forcing'];
elseif summary.fixed_domain_numerical_gate_pass
    summary.interpretation = ['candidate C passes the fixed L=40 numerical ', ...
        'screen but is domain-sensitive; proceed to L=40 forcing with ', ...
        'domain-qualified claims'];
else
    summary.interpretation = ['candidate C fails or remains inconclusive under ', ...
        'a fixed-domain numerical factor; inspect before forcing reachability'];
end
end

function make_figures(outdir, records, cases, cfg)
fig = figure('Color', 'w', 'Position', [100 100 1100 520]);
tiledlayout(1,2,'TileSpacing','compact');
nexttile;
bar([records.relative_period_error_from_source_C]); hold on;
yline(cfg.candidate_c_validation.maximum_relative_period_error, '--', ...
    'identification threshold');
set(gca, 'XTick', 1:numel(records), 'XTickLabel', {cases.id}, ...
    'XTickLabelRotation', 35);
ylabel('Relative period error from source C'); grid on; box on;
nexttile;
bar([records.standardized_profile_distance_to_source_C]); hold on;
yline(cfg.candidate_c_validation.maximum_standardized_profile_distance, '--', ...
    'identification threshold');
set(gca, 'XTick', 1:numel(records), 'XTickLabel', {cases.id}, ...
    'XTickLabelRotation', 35);
ylabel('Standardized section-profile distance'); grid on; box on;
exportgraphics(fig, fullfile(outdir, 'candidate_c_validation_comparison.png'), ...
    'Resolution', 300);
close(fig);
end
