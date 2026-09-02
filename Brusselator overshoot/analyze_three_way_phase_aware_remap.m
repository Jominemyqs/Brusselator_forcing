function results = analyze_three_way_phase_aware_remap()
%ANALYZE_THREE_WAY_PHASE_AWARE_REMAP Reclassify saved outcomes as A, B, or C.
%   A is stationary; B and C are periodic reference templates. B/C distances
%   minimize over phase and spatial reflection. Frozen slice trajectories use
%   full saved (u,v) fields. The legacy long overshoot trajectories saved only
%   v(x,t), so their late-window labels are explicitly V-only and their saved
%   full final snapshot is classified separately.

study_id = 'three_way_phase_aware_remap_v1';
slice_file = fullfile('experiment_outputs', 'phase_aware_basin_slice_v1', ...
    'raw_phase_aware_basin_slice.mat');
extension_file = fullfile('experiment_outputs', ...
    'phase_aware_basin_slice_edge_extension_v1', ...
    'raw_edge_extension_trajectories.mat');
candidate_file = fullfile('experiment_outputs', ...
    'edge_candidate_frozen_confirmation_v1', ...
    'candidate_confirmation_trajectory.mat');
overshoot_file = fullfile('experiment_outputs', ...
    'extended_postforcing_reference_v2', 'raw_extended_trajectories.mat');
outdir = fullfile('experiment_outputs', study_id);
input_files = {slice_file, extension_file, candidate_file, overshoot_file};
for k = 1:numel(input_files)
    if ~isfile(input_files{k})
        error('analyze_three_way_phase_aware_remap:MissingInput', ...
            'Expected saved input at %s.', input_files{k});
    end
end
if isfolder(outdir)
    error('analyze_three_way_phase_aware_remap:OutputExists', ...
        'Refusing to overwrite existing analysis directory: %s', outdir);
end

timer = tic;
loaded = load(slice_file, 'results');
slice = loaded.results;
loaded = load(extension_file, 'results');
extension = loaded.results;
loaded = load(candidate_file, 'results');
candidate = loaded.results;
loaded = load(overshoot_file, 'raw');
overshoot = loaded.raw;

cfg = slice.configuration;
cfg.experiment_name = study_id;
cfg.three_way_remap = struct( ...
    'slice_file', slice_file, ...
    'extension_file', extension_file, ...
    'candidate_file', candidate_file, ...
    'overshoot_file', overshoot_file, ...
    'outcomes', {{'A_stationary', 'B_periodic_reflection', 'C_periodic_direct'}}, ...
    'full_state_metric', ['relative_l2 full state; A modulo reflection; ', ...
        'B/C modulo phase and reflection'], ...
    'overshoot_v_only_late_window', 1000, ...
    'overshoot_v_only_sample_interval', 5, ...
    'classification_threshold', cfg.classification.recovery_threshold, ...
    'late_maximum_threshold', cfg.basin_slice.late_max_threshold, ...
    'limitations', ['Long overshoot output lacks U(x,t); its full-state ', ...
        'classification is available only at the saved final snapshot.']);

mkdir(outdir);
brusselator_write_metadata(outdir, cfg, ...
    brusselator_run_metadata(cfg, mfilename, NaN));
refs = construct_references(slice, candidate, cfg);
[slice_records, slice_histories] = remap_slice(slice, extension, refs, cfg);
[overshoot_records, overshoot_histories] = remap_overshoot(overshoot, refs, cfg);

slice_table = struct2table(slice_records, 'AsArray', true);
overshoot_table = struct2table(overshoot_records, 'AsArray', true);
writetable(slice_table, fullfile(outdir, 'three_way_slice_remap.csv'));
writetable(overshoot_table, fullfile(outdir, 'three_way_overshoot_remap.csv'));
make_figures(outdir, slice_records, slice_histories, overshoot_records, cfg);

summary = struct();
summary.slice_records = numel(slice_records);
summary.slice_A = sum(strcmp({slice_records.outcome}, 'A_stationary_neighborhood'));
summary.slice_B = sum(strcmp({slice_records.outcome}, 'B_periodic_reflection_neighborhood'));
summary.slice_C = sum(strcmp({slice_records.outcome}, 'C_periodic_direct_neighborhood'));
summary.slice_unresolved = sum(strcmp({slice_records.outcome}, 'unresolved_within_window'));
summary.overshoot_records = numel(overshoot_records);
summary.overshoot_A = sum(strcmp({overshoot_records.combined_outcome}, 'A_stationary_neighborhood'));
summary.overshoot_B = sum(strcmp({overshoot_records.combined_outcome}, 'B_periodic_reflection_neighborhood'));
summary.overshoot_C = sum(strcmp({overshoot_records.combined_outcome}, 'C_periodic_direct_neighborhood'));
summary.overshoot_inconclusive = sum(strcmp({overshoot_records.combined_outcome}, ...
    'inconclusive_due_to_reduced_long_trajectory_output'));
writetable(struct2table(summary, 'AsArray', true), ...
    fullfile(outdir, 'three_way_remap_summary.csv'));

results = struct('configuration', cfg, 'references', refs, ...
    'slice_records', slice_records, 'slice_histories', {slice_histories}, ...
    'overshoot_records', overshoot_records, 'overshoot_histories', {overshoot_histories}, ...
    'summary', summary);
save(fullfile(outdir, 'three_way_phase_aware_remap.mat'), 'results', '-v7.3');
metadata = brusselator_run_metadata(cfg, mfilename, toc(timer));
brusselator_write_metadata(outdir, cfg, metadata);
disp(slice_table);
disp(overshoot_table);
disp(struct2table(summary, 'AsArray', true));
fprintf('Three-way phase-aware remap complete. Results saved in: %s\n', outdir);
end

function refs = construct_references(slice, candidate, cfg)
refs = struct('A', slice.endpoint_A, 'B', slice.template);
S = candidate.trajectory.S;
U = candidate.trajectory.U;
V = candidate.trajectory.V;
period = candidate.section.period;
times = candidate.section.times(:);
valid = find(times >= 0.55 * S(end) & times + period <= S(end));
if isempty(valid)
    valid = find(times + period <= S(end));
end
if isempty(valid)
    error('analyze_three_way_phase_aware_remap:TemplateWindow', ...
        'No complete candidate-C reference period is available.');
end
index = valid(1);
phase_times = (0:cfg.basin_slice.template_phase_step:period)';
if phase_times(end) < period
    phase_times(end+1,1) = period;
end
start_time = times(index);
state = candidate.section.states(index,:);
N = cfg.grid.N;
refs.C = struct('phase_times', phase_times, ...
    'U', interp1(S, U, start_time + phase_times, 'linear'), ...
    'V', interp1(S, V, start_time + phase_times, 'linear'));
refs.C.U(1,:) = state(1:N);
refs.C.V(1,:) = state(N+1:end);
refs.info = struct('B_period', refs.B.phase_times(end), ...
    'C_period', period, 'C_section_time', start_time);
end

function [records, histories] = remap_slice(slice, extension, refs, cfg)
template = slice_record_template();
records = template([]);
histories = cell(0,1);
for k = 1:numel(slice.records)
    data = slice.trajectory_data{k};
    [record, history] = classify_full(data.S, data.U, data.V, refs, ...
        cfg.basin_slice.late_window, cfg.basin_slice.distance_sample_interval, cfg, template);
    record.source = 'initial_slice';
    record.source_stage = slice.records(k).stage;
    record.lambda = slice.records(k).lambda;
    record.trajectory_end_time = data.S(end);
    records(end+1,1) = record; %#ok<AGROW>
    histories{end+1,1} = history; %#ok<AGROW>
end
for k = 1:numel(extension.records)
    data = extension.trajectory_data{k};
    [S, U, V] = join_extension(data);
    [record, history] = classify_full(S, U, V, refs, ...
        cfg.basin_slice.late_window, cfg.basin_slice.distance_sample_interval, cfg, template);
    record.source = 'edge_extension';
    record.source_stage = 'continued_unresolved_case';
    record.lambda = extension.records(k).lambda;
    record.trajectory_end_time = S(end);
    records(end+1,1) = record; %#ok<AGROW>
    histories{end+1,1} = history; %#ok<AGROW>
end
end

function [S, U, V] = join_extension(data)
S = [data.old_S; data.old_S(end) + data.continuation_S(2:end)];
U = [data.old_U; data.continuation_U(2:end,:)];
V = [data.old_V; data.continuation_V(2:end,:)];
if size(U,1) ~= numel(S) || size(V,1) ~= numel(S)
    error('analyze_three_way_phase_aware_remap:SegmentMismatch', ...
        'Saved edge-extension trajectory segments are inconsistent.');
end
end

function [record, history] = classify_full(S, U, V, refs, late_window, sample_interval, cfg, template)
dt = median(diff(S));
indices = unique([1:max(1,round(sample_interval/dt)):numel(S), numel(S)]);
times = S(indices);
dA = zeros(numel(indices),1); dB = dA; dC = dA;
for k = 1:numel(indices)
    dA(k) = distance_A(U(indices(k),:), V(indices(k),:), refs.A);
    dB(k) = distance_orbit(U(indices(k),:), V(indices(k),:), refs.B);
    dC(k) = distance_orbit(U(indices(k),:), V(indices(k),:), refs.C);
end
late = times >= S(end) - late_window;
if ~any(late)
    error('analyze_three_way_phase_aware_remap:LateWindowEmpty', ...
        'No sampled points occur in the requested late window.');
end
record = template;
record.metric = 'full_state_phase_reflection_aware';
record.late_window = late_window;
record.late_median_distance_A = median(dA(late)); record.late_maximum_distance_A = max(dA(late));
record.late_median_distance_B = median(dB(late)); record.late_maximum_distance_B = max(dB(late));
record.late_median_distance_C = median(dC(late)); record.late_maximum_distance_C = max(dC(late));
record.final_distance_A = dA(end); record.final_distance_B = dB(end); record.final_distance_C = dC(end);
record.outcome = classify_window([record.late_median_distance_A, record.late_median_distance_B, record.late_median_distance_C], ...
    [record.late_maximum_distance_A, record.late_maximum_distance_B, record.late_maximum_distance_C], cfg);
history = struct('sample_times', times, 'distance_A', dA, 'distance_B', dB, 'distance_C', dC);
end

function record = slice_record_template()
record = struct('source', '', 'source_stage', '', 'lambda', NaN, ...
    'trajectory_end_time', NaN, 'metric', '', 'late_window', NaN, ...
    'late_median_distance_A', NaN, 'late_maximum_distance_A', NaN, ...
    'late_median_distance_B', NaN, 'late_maximum_distance_B', NaN, ...
    'late_median_distance_C', NaN, 'late_maximum_distance_C', NaN, ...
    'final_distance_A', NaN, 'final_distance_B', NaN, 'final_distance_C', NaN, ...
    'outcome', 'unclassified');
end

function [records, histories] = remap_overshoot(raw, refs, cfg)
template = overshoot_record_template();
records = repmat(template, 1 + numel(raw.forced_trajectories), 1);
histories = cell(numel(records),1);
[records(1), histories{1}] = classify_overshoot_case('control', NaN, ...
    raw.control.S, raw.control.V, raw.control.ic_final, refs, cfg, template);
for k = 1:numel(raw.forced_trajectories)
    forced = raw.forced_trajectories{k};
    [records(k+1), histories{k+1}] = classify_overshoot_case('forced', forced.bmax, ...
        forced.S, forced.V, forced.ic_final, refs, cfg, template);
end
end

function [record, history] = classify_overshoot_case(kind, bmax, S, V, final_state, refs, cfg, template)
late_window = cfg.three_way_remap.overshoot_v_only_late_window;
sample_interval = cfg.three_way_remap.overshoot_v_only_sample_interval;
indices = find(S >= S(end) - late_window);
stride = max(1, round(sample_interval / median(diff(S))));
indices = unique([indices(1:stride:end); indices(end)]);
times = S(indices);
dA = zeros(numel(indices),1); dB = dA; dC = dA;
for k = 1:numel(indices)
    dA(k) = distance_A_v(V(indices(k),:), refs.A);
    dB(k) = distance_orbit_v(V(indices(k),:), refs.B);
    dC(k) = distance_orbit_v(V(indices(k),:), refs.C);
end
u_final = final_state(:,2)';
v_final = final_state(:,3)';
full_distances = [distance_A(u_final, v_final, refs.A), ...
    distance_orbit(u_final, v_final, refs.B), ...
    distance_orbit(u_final, v_final, refs.C)];

record = template;
record.kind = kind;
record.bmax = bmax;
record.trajectory_end_time = S(end);
record.v_only_late_window = late_window;
record.v_only_sample_interval = sample_interval;
record.late_median_v_distance_A = median(dA); record.late_maximum_v_distance_A = max(dA);
record.late_median_v_distance_B = median(dB); record.late_maximum_v_distance_B = max(dB);
record.late_median_v_distance_C = median(dC); record.late_maximum_v_distance_C = max(dC);
record.v_only_outcome = classify_window([record.late_median_v_distance_A, record.late_median_v_distance_B, record.late_median_v_distance_C], ...
    [record.late_maximum_v_distance_A, record.late_maximum_v_distance_B, record.late_maximum_v_distance_C], cfg);
record.final_full_distance_A = full_distances(1);
record.final_full_distance_B = full_distances(2);
record.final_full_distance_C = full_distances(3);
record.final_full_snapshot_outcome = classify_snapshot(full_distances, cfg);
record.combined_outcome = combine_overshoot_outcomes(record);
history = struct('sample_times', times, 'v_distance_A', dA, 'v_distance_B', dB, 'v_distance_C', dC);
end

function record = overshoot_record_template()
record = struct('kind', '', 'bmax', NaN, 'trajectory_end_time', NaN, ...
    'v_only_late_window', NaN, 'v_only_sample_interval', NaN, ...
    'late_median_v_distance_A', NaN, 'late_maximum_v_distance_A', NaN, ...
    'late_median_v_distance_B', NaN, 'late_maximum_v_distance_B', NaN, ...
    'late_median_v_distance_C', NaN, 'late_maximum_v_distance_C', NaN, ...
    'v_only_outcome', 'unclassified', ...
    'final_full_distance_A', NaN, 'final_full_distance_B', NaN, ...
    'final_full_distance_C', NaN, 'final_full_snapshot_outcome', 'unclassified', ...
    'combined_outcome', 'unclassified');
end

function label = classify_window(medians, maxima, cfg)
valid = medians < cfg.classification.recovery_threshold & ...
    maxima < cfg.three_way_remap.late_maximum_threshold;
label = endpoint_label(valid);
end

function label = classify_snapshot(distances, cfg)
label = endpoint_label(distances < cfg.classification.recovery_threshold);
end

function label = endpoint_label(valid)
names = {'A_stationary_neighborhood', 'B_periodic_reflection_neighborhood', ...
    'C_periodic_direct_neighborhood'};
if sum(valid) == 1
    label = names{find(valid,1)};
elseif sum(valid) > 1
    label = 'ambiguous_multiple_endpoint_neighborhoods';
else
    label = 'unresolved_within_window';
end
end

function label = combine_overshoot_outcomes(record)
if strcmp(record.v_only_outcome, record.final_full_snapshot_outcome) && ...
        ~strcmp(record.v_only_outcome, 'unresolved_within_window')
    label = record.v_only_outcome;
else
    label = 'inconclusive_due_to_reduced_long_trajectory_output';
end
end

function distance = distance_A(u, v, A)
direct = sqrt(sum((u - A.U).^2) + sum((v - A.V).^2));
reflected = sqrt(sum((u - fliplr(A.U)).^2) + sum((v - fliplr(A.V)).^2));
distance = min(direct, reflected) / max(sqrt(sum(A.U.^2) + sum(A.V.^2)), eps);
end

function distance = distance_A_v(v, A)
direct = norm(v - A.V) / max(norm(A.V), eps);
reflected = norm(v - fliplr(A.V)) / max(norm(A.V), eps);
distance = min(direct, reflected);
end

function distance = distance_orbit(u, v, template)
direct = sqrt(sum((template.U - u).^2 + (template.V - v).^2, 2)) ./ ...
    max(sqrt(sum(template.U.^2 + template.V.^2, 2)), eps);
reflected = sqrt(sum((fliplr(template.U) - u).^2 + (fliplr(template.V) - v).^2, 2)) ./ ...
    max(sqrt(sum(template.U.^2 + template.V.^2, 2)), eps);
distance = min([direct; reflected]);
end

function distance = distance_orbit_v(v, template)
direct = sqrt(sum((template.V - v).^2, 2)) ./ max(sqrt(sum(template.V.^2, 2)), eps);
reflected = sqrt(sum((fliplr(template.V) - v).^2, 2)) ./ ...
    max(sqrt(sum(template.V.^2, 2)), eps);
distance = min([direct; reflected]);
end

function make_figures(outdir, slice_records, slice_histories, overshoot_records, cfg)
fig = figure('Color', 'w', 'Position', [100 100 1000 550]);
hold on;
sources = {'initial_slice', 'edge_extension'};
markers = {'o', 's'};
for q = 1:numel(sources)
    idx = strcmp({slice_records.source}, sources{q});
    if ~any(idx)
        continue;
    end
    lambda = [slice_records(idx).lambda];
    semilogy(lambda, max([slice_records(idx).late_median_distance_A], eps), ...
        [markers{q}, '-'], 'LineWidth', 1.0, 'DisplayName', [sources{q}, ': A']);
    semilogy(lambda, max([slice_records(idx).late_median_distance_B], eps), ...
        [markers{q}, '--'], 'LineWidth', 1.0, 'DisplayName', [sources{q}, ': B']);
    semilogy(lambda, max([slice_records(idx).late_median_distance_C], eps), ...
        [markers{q}, ':'], 'LineWidth', 1.4, 'DisplayName', [sources{q}, ': C']);
end
yline(cfg.classification.recovery_threshold, '--', 'classification threshold', ...
    'HandleVisibility', 'off');
xlabel('\lambda in X_\lambda = (1-\lambda)A + \lambda B');
ylabel('Late median relative full-state distance');
title('Three-way phase-aware remap of saved frozen slice trajectories');
legend('Location', 'eastoutside', 'FontSize', 7); grid on; box on;
exportgraphics(fig, fullfile(outdir, 'three_way_slice_late_distances.png'), ...
    'Resolution', 300);
close(fig);

fig = figure('Color', 'w', 'Position', [100 100 1100 650]);
tiledlayout(3,1,'TileSpacing','compact');
labels = {'A', 'B', 'C'};
for q = 1:numel(labels)
    nexttile; hold on;
    field = ['distance_', labels{q}];
    for k = 1:numel(slice_histories)
        history = slice_histories{k};
        semilogy(history.sample_times, max(history.(field), eps), 'LineWidth', 0.8, ...
            'DisplayName', sprintf('\\lambda=%.3f, %s', ...
            slice_records(k).lambda, slice_records(k).source));
    end
    yline(cfg.classification.recovery_threshold, '--', 'HandleVisibility', 'off');
    ylabel(['d_', labels{q}]); grid on; box on;
    if q == 1
        title('Saved frozen slice trajectories under three-way distances');
        legend('Location', 'eastoutside', 'FontSize', 6);
    end
end
xlabel('Frozen b=10 integration time');
exportgraphics(fig, fullfile(outdir, 'three_way_slice_distance_histories.png'), ...
    'Resolution', 300);
close(fig);

fig = figure('Color', 'w', 'Position', [100 100 900 500]);
x = 1:numel(overshoot_records);
bar(x - 0.22, [overshoot_records.late_median_v_distance_A], 0.22, ...
    'DisplayName', 'A'); hold on;
bar(x, [overshoot_records.late_median_v_distance_B], 0.22, ...
    'DisplayName', 'B');
bar(x + 0.22, [overshoot_records.late_median_v_distance_C], 0.22, ...
    'DisplayName', 'C');
set(gca, 'YScale', 'log');
yline(cfg.classification.recovery_threshold, '--', 'threshold', 'HandleVisibility', 'off');
names = cell(numel(overshoot_records),1);
for k = 1:numel(overshoot_records)
    if strcmp(overshoot_records(k).kind, 'control')
        names{k} = 'control';
    else
        names{k} = sprintf('b_{max}=%.2f', overshoot_records(k).bmax);
    end
end
xticks(x); xticklabels(names);
ylabel('Late median V-only distance');
title('Saved long overshoot trajectories: three-way V-only remap');
legend('Location', 'best'); grid on; box on;
exportgraphics(fig, fullfile(outdir, 'three_way_overshoot_v_only_distances.png'), ...
    'Resolution', 300);
close(fig);
end
