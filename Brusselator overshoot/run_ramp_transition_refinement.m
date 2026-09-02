function results = run_ramp_transition_refinement()
%RUN_RAMP_TRANSITION_REFINEMENT Resolve structured symmetric-ramp changes.
%   At bmax=11.14 and Thold=80, the pilot found unresolved behavior at ramp
%   time 10, B at 20 and 40, and A at 80. This study samples every integer
%   ramp time in (10,20) and every five units in (20,80), reusing the saved
%   ramp-40 result rather than rerunning it. Unresolved cases are followed to
%   240 post-forcing units; the subsequent long-extension study treats all of
%   them at a common 1000-unit horizon.

study_id = 'ramp_transition_refinement_v1';
parent_file = fullfile('experiment_outputs', ...
    'forcing_to_c_reachability_pilot_v1', 'forcing_to_c_reachability.mat');
outdir = fullfile('experiment_outputs', study_id);
progress_file = fullfile(outdir, 'progress_ramp_refinement.mat');
if ~isfile(parent_file)
    error('run_ramp_transition_refinement:MissingParent', ...
        'Expected forcing pilot at %s.', parent_file);
end
parent_data = load(parent_file, 'results');
parent = parent_data.results;
[refs, reference_cfg] = brusselator_three_way_references();
protocols = refined_protocols();
moderate_solver = parent.configuration.reachability.discovery_solver;
tight_solver = parent.configuration.reachability.C_confirmation_solver;
cfg = brusselator_default_config(struct( ...
    'experiment_name', study_id, 'model', reference_cfg.model, ...
    'grid', reference_cfg.grid, 'initial', reference_cfg.initial, ...
    'solver', moderate_solver));
cfg.reachability = parent.configuration.reachability;
cfg.reachability.parent_file = parent_file;
cfg.reachability.protocol_count = numel(protocols);
cfg.reachability.refined_intervals = [10, 20; 20, 80];
cfg.reachability.ramp_values = [protocols.Tup];
cfg.reachability.interpretation = ['categorical refinement of the physically ', ...
    'realized symmetric-ramp outcome changes at bmax=11.14, Thold=80'];

resume = isfolder(outdir) && isfile(progress_file);
if isfolder(outdir) && ~resume
    error('run_ramp_transition_refinement:OutputExists', ...
        ['Refusing to use an existing output directory without a valid ', ...
        'checkpoint: %s'], outdir);
end
record_template = empty_record();
if resume
    checkpoint = load(progress_file, 'cfg', 'records', 'raw_files');
    if numel(checkpoint.records) ~= numel(protocols) || ...
            ~strcmp(checkpoint.cfg.reachability.parent_file, parent_file)
        error('run_ramp_transition_refinement:CheckpointMismatch', ...
            'The checkpoint does not match this ramp-refinement matrix.');
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
        fprintf('Skipping completed ramp case %d/%d: T_ramp=%.0f.\n', ...
            k, numel(protocols), protocol.Tup);
        continue;
    end
    fprintf('Ramp refinement %d/%d: symmetric T_ramp=%.0f ...\n', ...
        k, numel(protocols), protocol.Tup);
    case_timer = tic;
    discovery = brusselator_execute_forcing_protocol(common_initial_state, ...
        protocol, cfg, moderate_solver, refs, ...
        cfg.reachability.initial_post_forcing_duration, ...
        cfg.reachability.extended_post_forcing_duration);
    tight_confirmation = struct();
    if strcmp(discovery.classification.outcome, ...
            'C_periodic_direct_neighborhood')
        fprintf('  Candidate C hit; repeating with the tight solver ...\n');
        tight_confirmation = brusselator_execute_forcing_protocol( ...
            common_initial_state, protocol, cfg, tight_solver, refs, ...
            cfg.reachability.initial_post_forcing_duration, ...
            cfg.reachability.extended_post_forcing_duration);
    end
    runtime = toc(case_timer);
    record = summarize(protocol, discovery, tight_confirmation, runtime);
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
writetable(summary_table, fullfile(outdir, 'ramp_refinement_summary.csv'));
[combined, interval_summary] = combine_with_parent(records, parent);
writetable(struct2table(combined, 'AsArray', true), ...
    fullfile(outdir, 'combined_symmetric_ramp_summary.csv'));
writetable(struct2table(interval_summary, 'AsArray', true), ...
    fullfile(outdir, 'ramp_transition_intervals.csv'));
make_figures(outdir, combined);
summary = aggregate(records);
writetable(struct2table(summary, 'AsArray', true), ...
    fullfile(outdir, 'ramp_refinement_counts.csv'));
results = struct('configuration', cfg, 'protocols', protocols, ...
    'records', records, 'combined_records', combined, ...
    'interval_summary', interval_summary, 'summary', summary, ...
    'raw_files', {raw_files});
save(fullfile(outdir, 'ramp_transition_refinement.mat'), 'results', '-v7.3');
metadata = brusselator_run_metadata(cfg, mfilename, sum([records.runtime_seconds]));
brusselator_write_metadata(outdir, cfg, metadata);
disp(struct2table(combined, 'AsArray', true));
disp(struct2table(interval_summary, 'AsArray', true));
fprintf('Ramp-transition refinement complete. Results saved in: %s\n', outdir);
end

function protocols = refined_protocols()
ramps = [11:19, 25:5:35, 45:5:75];
template = struct('id', '', 'family', 'symmetric_ramp_refinement', ...
    'bmax', 11.14, 'Tup', NaN, 'Thold', 80, 'Tdown', NaN);
protocols = repmat(template, numel(ramps), 1);
for k = 1:numel(ramps)
    protocols(k).id = sprintf('symmetric_ramp_%03d', ramps(k));
    protocols(k).Tup = ramps(k);
    protocols(k).Tdown = ramps(k);
end
end

function record = empty_record()
record = struct('source', 'refinement', 'protocol_id', '', ...
    'ramp_time', NaN, 'bmax', NaN, 'Thold', NaN, 'forcing_end', NaN, ...
    'post_forcing_duration', NaN, 'landing_distance_A', NaN, ...
    'landing_distance_B', NaN, 'landing_distance_C', NaN, ...
    'late_median_distance_A', NaN, 'late_median_distance_B', NaN, ...
    'late_median_distance_C', NaN, 'outcome', 'unclassified', ...
    'tight_confirmation_outcome', 'not_run', 'confirmed_C', false, ...
    'runtime_seconds', NaN);
end

function record = summarize(protocol, discovery, confirmation, runtime)
record = empty_record();
record.protocol_id = protocol.id;
record.ramp_time = protocol.Tup;
record.bmax = protocol.bmax;
record.Thold = protocol.Thold;
record.forcing_end = discovery.forcing_end;
record.post_forcing_duration = discovery.post_forcing_duration;
record.landing_distance_A = discovery.landing_distances(1);
record.landing_distance_B = discovery.landing_distances(2);
record.landing_distance_C = discovery.landing_distances(3);
record.late_median_distance_A = discovery.classification.late_median_distance_A;
record.late_median_distance_B = discovery.classification.late_median_distance_B;
record.late_median_distance_C = discovery.classification.late_median_distance_C;
record.outcome = discovery.classification.outcome;
if ~isempty(fieldnames(confirmation))
    record.tight_confirmation_outcome = confirmation.classification.outcome;
    record.confirmed_C = strcmp(record.outcome, ...
        'C_periodic_direct_neighborhood') && strcmp( ...
        record.tight_confirmation_outcome, 'C_periodic_direct_neighborhood');
end
record.runtime_seconds = runtime;
end

function [combined, intervals] = combine_with_parent(records, parent)
template = empty_record();
combined = records;
parent_ramps = [10, 20, 40, 80];
for ramp = parent_ramps
    if ramp == 40
        parent_id = 'amplitude_bmax_11.14';
    else
        parent_id = sprintf('symmetric_ramp_%03d', ramp);
    end
    index = find(strcmp({parent.records.protocol_id}, parent_id), 1);
    if isempty(index)
        error('run_ramp_transition_refinement:MissingParentRamp', ...
            'The parent study lacks symmetric ramp %d.', ramp);
    end
    old = parent.records(index);
    record = template;
    record.source = 'parent_pilot';
    record.protocol_id = old.protocol_id;
    record.ramp_time = old.Tup;
    record.bmax = old.bmax;
    record.Thold = old.Thold;
    record.forcing_end = old.forcing_end;
    record.post_forcing_duration = old.post_forcing_duration;
    record.landing_distance_A = old.landing_distance_A;
    record.landing_distance_B = old.landing_distance_B;
    record.landing_distance_C = old.landing_distance_C;
    record.late_median_distance_A = old.late_median_distance_A;
    record.late_median_distance_B = old.late_median_distance_B;
    record.late_median_distance_C = old.late_median_distance_C;
    record.outcome = old.outcome;
    record.tight_confirmation_outcome = old.tight_confirmation_outcome;
    record.confirmed_C = old.confirmed_C;
    record.runtime_seconds = old.runtime_seconds;
    combined(end+1,1) = record; %#ok<AGROW>
end
[~, order] = sort([combined.ramp_time]);
combined = combined(order);

interval_template = struct('left_ramp', NaN, 'left_outcome', '', ...
    'right_ramp', NaN, 'right_outcome', '', 'width', NaN, ...
    'interpretation', 'adjacent sampled protocols with different outcomes');
intervals = interval_template([]);
for k = 1:(numel(combined)-1)
    if ~strcmp(combined(k).outcome, combined(k+1).outcome)
        item = interval_template;
        item.left_ramp = combined(k).ramp_time;
        item.left_outcome = combined(k).outcome;
        item.right_ramp = combined(k+1).ramp_time;
        item.right_outcome = combined(k+1).outcome;
        item.width = item.right_ramp - item.left_ramp;
        intervals(end+1,1) = item; %#ok<AGROW>
    end
end
end

function summary = aggregate(records)
summary = struct('number_new_protocols', numel(records), ...
    'number_A', sum(strcmp({records.outcome}, 'A_stationary_neighborhood')), ...
    'number_B', sum(strcmp({records.outcome}, ...
        'B_periodic_reflection_neighborhood')), ...
    'number_C', sum(strcmp({records.outcome}, ...
        'C_periodic_direct_neighborhood')), ...
    'number_unresolved', sum(strcmp({records.outcome}, ...
        'unresolved_within_window')), ...
    'number_confirmed_C', sum([records.confirmed_C]));
end

function make_figures(outdir, records)
ramps = [records.ramp_time];
codes = zeros(size(ramps));
for k = 1:numel(records)
    if strcmp(records(k).outcome, 'A_stationary_neighborhood')
        codes(k) = 1;
    elseif strcmp(records(k).outcome, 'B_periodic_reflection_neighborhood')
        codes(k) = 2;
    elseif strcmp(records(k).outcome, 'C_periodic_direct_neighborhood')
        codes(k) = 3;
    else
        codes(k) = 4;
    end
end
fig = figure('Color', 'w', 'Position', [100 100 950 430]);
scatter(ramps, codes, 55, codes, 'filled');
yticks(1:4); yticklabels({'A','B','C','unresolved'});
xlabel('Symmetric ramp time'); ylabel('Finite-time outcome');
title('Refined symmetric-ramp transitions at b_{max}=11.14, T_{hold}=80');
grid on; box on;
exportgraphics(fig, fullfile(outdir, 'ramp_outcome_map.png'), 'Resolution', 300);
close(fig);

fig = figure('Color', 'w', 'Position', [100 100 950 480]);
semilogy(ramps, [records.landing_distance_A], '-o', 'DisplayName', 'A'); hold on;
semilogy(ramps, [records.landing_distance_B], '-o', 'DisplayName', 'B');
semilogy(ramps, [records.landing_distance_C], '-o', 'DisplayName', 'C');
xlabel('Symmetric ramp time'); ylabel('Distance at return to b=10');
title('Landing-state distances across refined ramp times');
legend('Location','best'); grid on; box on;
exportgraphics(fig, fullfile(outdir, 'ramp_landing_distances.png'), 'Resolution', 300);
close(fig);
end
