function results = analyze_unresolved_recurrent_groups()
%ANALYZE_UNRESOLVED_RECURRENT_GROUPS Compare long recurrent forcing outcomes.
%   Builds one-period full-state templates from the late continuation of every
%   formerly unresolved trajectory. Pairwise distances minimize an RMS
%   full-state mismatch over temporal phase and exact spatial reflection.
%   Connected components below a documented distance threshold identify
%   recurrent groups reached by multiple forcing protocols.
%
%   These groups are candidate invariant outcomes. This analysis does not by
%   itself establish local attraction, solver robustness, or asymptotic
%   convergence.

study_id = 'unresolved_recurrent_group_analysis_v1';
parent_file = fullfile('experiment_outputs', ...
    'unresolved_forcing_long_extension_v1', ...
    'unresolved_forcing_long_extensions.mat');
outdir = fullfile('experiment_outputs', study_id);
if ~isfile(parent_file)
    error('analyze_unresolved_recurrent_groups:MissingParent', ...
        'Expected long-extension result at %s.', parent_file);
end
if isfolder(outdir)
    error('analyze_unresolved_recurrent_groups:OutputExists', ...
        'Refusing to overwrite existing analysis directory: %s', outdir);
end
parent_data = load(parent_file, 'results');
parent = parent_data.results;
cfg = parent.configuration;
cfg.experiment_name = study_id;
cfg.recurrent_group_analysis = struct( ...
    'parent_file', parent_file, 'phase_samples', 360, ...
    'same_orbit_threshold', 2e-2, ...
    'metric', ['symmetric relative RMS full-state distance minimized over ', ...
        'normalized temporal phase and spatial reflection'], ...
    'interpretation', ['cross-protocol recurrent-orbit grouping; not an ', ...
        'attractor or continuation computation']);
mkdir(outdir);
brusselator_write_metadata(outdir, cfg, ...
    brusselator_run_metadata(cfg, mfilename, NaN));

timer = tic;
templates = cell(numel(parent.records), 1);
for k = 1:numel(parent.records)
    saved = load(parent.continuation_files{k}, 'continuation');
    period = parent.records(k).section_period;
    templates{k} = late_period_template(saved.continuation, period, ...
        cfg.recurrent_group_analysis.phase_samples);
end
[distances, reflections, phase_shifts] = pairwise_distances(templates);
groups = connected_groups(distances, ...
    cfg.recurrent_group_analysis.same_orbit_threshold);
[membership, group_records] = summarize_groups(parent.records, distances, groups, cfg);

names = {parent.records.protocol_id};
distance_table = array2table(distances, ...
    'VariableNames', matlab.lang.makeValidName(names), 'RowNames', names);
writetable(distance_table, fullfile(outdir, 'pairwise_orbit_distances.csv'), ...
    'WriteRowNames', true);
reflection_table = array2table(reflections, ...
    'VariableNames', matlab.lang.makeValidName(names), 'RowNames', names);
writetable(reflection_table, fullfile(outdir, 'pairwise_best_reflections.csv'), ...
    'WriteRowNames', true);
shift_table = array2table(phase_shifts, ...
    'VariableNames', matlab.lang.makeValidName(names), 'RowNames', names);
writetable(shift_table, fullfile(outdir, 'pairwise_best_phase_fraction.csv'), ...
    'WriteRowNames', true);
writetable(struct2table(membership, 'AsArray', true), ...
    fullfile(outdir, 'recurrent_group_membership.csv'));
writetable(struct2table(group_records, 'AsArray', true), ...
    fullfile(outdir, 'recurrent_group_summary.csv'));
make_figures(outdir, distances, membership, names);

summary = struct('number_trajectories', numel(parent.records), ...
    'number_groups', numel(groups), ...
    'number_cross_protocol_groups', sum([group_records.number_protocols] >= 2), ...
    'maximum_within_group_distance', max([group_records.maximum_within_distance]), ...
    'minimum_between_group_distance', minimum_between_group_distance(distances, groups), ...
    'interpretation', ['phase/reflection-aware recurrent groups shared across ', ...
        'forcing protocols; candidates require independent numerical and ', ...
        'local-stability checks']);
writetable(struct2table(summary, 'AsArray', true), ...
    fullfile(outdir, 'recurrent_group_analysis_summary.csv'));
results = struct('configuration', cfg, 'parent_file', parent_file, ...
    'templates', {templates}, 'distance_matrix', distances, ...
    'best_reflection_matrix', reflections, ...
    'best_phase_fraction_matrix', phase_shifts, 'groups', {groups}, ...
    'membership', membership, 'group_records', group_records, ...
    'summary', summary);
save(fullfile(outdir, 'unresolved_recurrent_group_analysis.mat'), ...
    'results', '-v7.3');
metadata = brusselator_run_metadata(cfg, mfilename, toc(timer));
brusselator_write_metadata(outdir, cfg, metadata);
disp(struct2table(membership, 'AsArray', true));
disp(struct2table(group_records, 'AsArray', true));
disp(struct2table(summary, 'AsArray', true));
fprintf('Unresolved recurrent-group analysis saved in: %s\n', outdir);
end

function template = late_period_template(continuation, period, phase_samples)
S = continuation.S(:);
start_time = S(end) - period;
phase = (0:phase_samples-1)' / phase_samples;
query = start_time + phase * period;
template = struct('phase_fraction', phase, ...
    'U', interp1(S, continuation.U, query, 'linear'), ...
    'V', interp1(S, continuation.V, query, 'linear'), ...
    'period', period);
end

function [matrix, reflection_matrix, shift_matrix] = pairwise_distances(templates)
n = numel(templates);
matrix = zeros(n);
reflection_matrix = false(n);
shift_matrix = zeros(n);
for i = 1:n
    for j = i+1:n
        [distance, reflection, shift] = orbit_distance(templates{i}, templates{j});
        matrix(i,j) = distance;
        matrix(j,i) = distance;
        reflection_matrix(i,j) = reflection;
        reflection_matrix(j,i) = reflection;
        shift_matrix(i,j) = shift;
        shift_matrix(j,i) = mod(-shift, 1);
    end
end
end

function [best, uses_reflection, phase_fraction] = orbit_distance(A, B)
X = [A.U, A.V];
B_direct = [B.U, B.V];
B_reflected = [fliplr(B.U), fliplr(B.V)];
normalizer = sqrt(0.5 * (mean(sum(X.^2,2)) + mean(sum(B_direct.^2,2))));
n = size(X,1);
direct = zeros(n,1);
reflected = zeros(n,1);
for shift = 0:n-1
    direct(shift+1) = sqrt(mean(sum((X - circshift(B_direct, shift, 1)).^2, 2))) / ...
        max(normalizer, eps);
    reflected(shift+1) = sqrt(mean(sum((X - circshift(B_reflected, shift, 1)).^2, 2))) / ...
        max(normalizer, eps);
end
[direct_minimum, direct_index] = min(direct);
[reflection_minimum, reflection_index] = min(reflected);
if reflection_minimum < direct_minimum
    best = reflection_minimum;
    uses_reflection = true;
    phase_fraction = (reflection_index - 1) / n;
else
    best = direct_minimum;
    uses_reflection = false;
    phase_fraction = (direct_index - 1) / n;
end
end

function groups = connected_groups(distances, threshold)
n = size(distances,1);
unassigned = true(n,1);
groups = cell(0,1);
while any(unassigned)
    seed = find(unassigned, 1);
    members = seed;
    frontier = seed;
    unassigned(seed) = false;
    while ~isempty(frontier)
        neighbors = find(any(distances(frontier,:) < threshold, 1))';
        neighbors = neighbors(unassigned(neighbors));
        if isempty(neighbors)
            break;
        end
        members = unique([members; neighbors]);
        frontier = neighbors;
        unassigned(neighbors) = false;
    end
    groups{end+1,1} = members; %#ok<AGROW>
end
end

function [membership, summaries] = summarize_groups(records, distances, groups, cfg)
member_template = struct('protocol_id', '', 'group_id', NaN, ...
    'candidate_label', '', 'period', NaN, 'distance_to_group_representative', NaN);
membership = repmat(member_template, numel(records), 1);
summary_template = struct('group_id', NaN, 'candidate_label', '', ...
    'number_protocols', NaN, 'representative_protocol', '', ...
    'median_period', NaN, 'maximum_relative_period_spread', NaN, ...
    'maximum_within_distance', NaN, 'cross_protocol_reproducible', false, ...
    'interpretation', ['same recurrent orbit modulo sampled phase/reflection; ', ...
        'candidate only']);
summaries = repmat(summary_template, numel(groups), 1);
for g = 1:numel(groups)
    members = groups{g};
    label = sprintf('candidate_group_%c', char('D' + g - 1));
    representative = members(1);
    periods = [records(members).section_period];
    block = distances(members,members);
    summary = summary_template;
    summary.group_id = g;
    summary.candidate_label = label;
    summary.number_protocols = numel(members);
    summary.representative_protocol = records(representative).protocol_id;
    summary.median_period = median(periods);
    summary.maximum_relative_period_spread = ...
        max(abs(periods - median(periods))) / median(periods);
    summary.maximum_within_distance = max(block(:));
    summary.cross_protocol_reproducible = numel(members) >= 2 && ...
        summary.maximum_within_distance < ...
        cfg.recurrent_group_analysis.same_orbit_threshold;
    summaries(g) = summary;
    for index = members(:)'
        item = member_template;
        item.protocol_id = records(index).protocol_id;
        item.group_id = g;
        item.candidate_label = label;
        item.period = records(index).section_period;
        item.distance_to_group_representative = distances(index, representative);
        membership(index) = item;
    end
end
end

function value = minimum_between_group_distance(distances, groups)
value = Inf;
for i = 1:numel(groups)
    for j = i+1:numel(groups)
        block = distances(groups{i}, groups{j});
        value = min(value, min(block(:)));
    end
end
if isinf(value)
    value = NaN;
end
end

function make_figures(outdir, distances, membership, names)
fig = figure('Color', 'w', 'Position', [100 100 800 700]);
imagesc(distances); axis image; colorbar;
set(gca, 'XTick', 1:numel(names), 'XTickLabel', names, ...
    'XTickLabelRotation', 45, 'YTick', 1:numel(names), 'YTickLabel', names);
title('Phase/reflection-aware distance between recurrent outcomes');
exportgraphics(fig, fullfile(outdir, 'recurrent_orbit_distance_matrix.png'), ...
    'Resolution', 300);
close(fig);

fig = figure('Color', 'w', 'Position', [100 100 950 420]);
scatter(1:numel(membership), [membership.period], 70, ...
    [membership.group_id], 'filled');
set(gca, 'XTick', 1:numel(names), 'XTickLabel', names, ...
    'XTickLabelRotation', 40);
ylabel('Late recurrent period'); title('Cross-protocol recurrent groups');
grid on; box on;
exportgraphics(fig, fullfile(outdir, 'recurrent_group_periods.png'), ...
    'Resolution', 300);
close(fig);
end
