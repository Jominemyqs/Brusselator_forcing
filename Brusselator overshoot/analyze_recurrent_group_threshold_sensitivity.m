function results = analyze_recurrent_group_threshold_sensitivity()
%ANALYZE_RECURRENT_GROUP_THRESHOLD_SENSITIVITY Audit candidate group labels.
%   Reclusters the phase/reflection-aware pairwise orbit distances over a
%   documented range of cutoffs. This prevents a single exploratory distance
%   threshold from being mistaken for evidence of a definite attractor count.
%
%   The two period classes are also compared directly. All labels remain
%   finite-time recurrent-pattern candidates rather than verified attractors.

study_id = 'unresolved_recurrent_threshold_sensitivity_v1';
parent_file = fullfile('experiment_outputs', ...
    'unresolved_recurrent_group_analysis_v1', ...
    'unresolved_recurrent_group_analysis.mat');
outdir = fullfile('experiment_outputs', study_id);
if ~isfile(parent_file)
    error('analyze_recurrent_group_threshold_sensitivity:MissingParent', ...
        'Expected recurrent-group result at %s.', parent_file);
end
if isfolder(outdir)
    error('analyze_recurrent_group_threshold_sensitivity:OutputExists', ...
        'Refusing to overwrite existing analysis directory: %s', outdir);
end

loaded = load(parent_file, 'results');
parent = loaded.results;
distances = parent.distance_matrix;
records = parent.membership;
names = string({records.protocol_id})';
periods = [records.period]';
thresholds = [0.005; 0.01; 0.015; 0.02; 0.025; 0.03; 0.04; 0.05; 0.075; 0.1];

cfg = parent.configuration;
cfg.experiment_name = study_id;
cfg.threshold_sensitivity = struct( ...
    'parent_file', parent_file, 'thresholds', thresholds, ...
    'interpretation', ['exploratory connected-component sensitivity; ', ...
        'not a count of invariant attractors']);
mkdir(outdir);
brusselator_write_metadata(outdir, cfg, ...
    brusselator_run_metadata(cfg, mfilename, NaN));
timer = tic;

sweep_template = struct('distance_threshold', NaN, 'number_groups', NaN, ...
    'number_cross_protocol_groups', NaN, 'largest_group_size', NaN);
sweep = repmat(sweep_template, numel(thresholds), 1);
membership_rows = cell(0,4);
for k = 1:numel(thresholds)
    groups = connected_groups(distances, thresholds(k));
    sizes = cellfun(@numel, groups);
    sweep(k) = struct('distance_threshold', thresholds(k), ...
        'number_groups', numel(groups), ...
        'number_cross_protocol_groups', sum(sizes >= 2), ...
        'largest_group_size', max(sizes));
    for g = 1:numel(groups)
        members = groups{g};
        for index = members(:)'
            membership_rows(end+1,:) = {thresholds(k), names(index), g, numel(members)}; %#ok<AGROW>
        end
    end
end

short = periods < 5.8;
long = periods >= 5.8;
period_class = struct( ...
    'short_period_count', sum(short), ...
    'short_period_median', median(periods(short)), ...
    'short_class_maximum_distance', maximum_block_distance(distances, short), ...
    'long_period_count', sum(long), ...
    'long_period_median', median(periods(long)), ...
    'long_class_maximum_distance', maximum_block_distance(distances, long), ...
    'minimum_between_period_class_distance', min(distances(short,long), [], 'all'), ...
    'interpretation', ['two strongly separated recurrence-period classes; ', ...
        'within-class orbit multiplicity remains unresolved']);

sweep_table = struct2table(sweep, 'AsArray', true);
membership_table = cell2table(membership_rows, 'VariableNames', ...
    {'distance_threshold','protocol_id','group_id','group_size'});
period_class_table = struct2table(period_class, 'AsArray', true);
writetable(sweep_table, fullfile(outdir, 'threshold_sweep_summary.csv'));
writetable(membership_table, fullfile(outdir, 'threshold_sweep_membership.csv'));
writetable(period_class_table, fullfile(outdir, 'period_class_separation.csv'));

fig = figure('Color', 'w', 'Position', [100 100 760 430]);
yyaxis left;
stairs([sweep.distance_threshold], [sweep.number_groups], '-o', ...
    'LineWidth', 1.4, 'MarkerFaceColor', 'auto');
ylabel('Connected candidate groups');
yyaxis right;
stairs([sweep.distance_threshold], [sweep.largest_group_size], '-s', ...
    'LineWidth', 1.4, 'MarkerFaceColor', 'auto');
ylabel('Largest group size');
xlabel('Full-state orbit-distance threshold');
title('Sensitivity of recurrent-pattern grouping to cutoff');
grid on; box on;
exportgraphics(fig, fullfile(outdir, 'recurrent_group_threshold_sensitivity.png'), ...
    'Resolution', 300);
close(fig);

results = struct('configuration', cfg, 'parent_file', parent_file, ...
    'threshold_sweep', sweep, 'threshold_membership', membership_table, ...
    'period_class_summary', period_class);
save(fullfile(outdir, 'recurrent_group_threshold_sensitivity.mat'), ...
    'results', '-v7.3');
brusselator_write_metadata(outdir, cfg, ...
    brusselator_run_metadata(cfg, mfilename, toc(timer)));
disp(sweep_table);
disp(period_class_table);
fprintf('Recurrent-group threshold sensitivity saved in: %s\n', outdir);
end

function value = maximum_block_distance(distances, selected)
block = distances(selected, selected);
value = max(block, [], 'all');
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
