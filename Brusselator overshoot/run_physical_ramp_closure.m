function results = run_physical_ramp_closure()
%RUN_PHYSICAL_RAMP_CLOSURE Refine the physical T_ramp=40--45 transition.
%   The fixed budget consists of an eleven-point half-unit seed scan and at
%   most nine adaptive midpoint runs. Midpoints are forcing protocols, not
%   interpolated states. The adaptive rule is multiclass-safe and never
%   coerces an A/B/C/R5/R6/unresolved/ambiguous outcome into a binary label.

cfg = brusselator_physical_ramp_closure_config();
study_id = cfg.experiment_name;
outdir = fullfile(cfg.output.root,study_id);
checkpoint_file = fullfile(outdir,'progress_physical_ramp_closure.mat');
final_file = fullfile(outdir,'physical_ramp_closure.mat');
if isfile(final_file)
    error('run_physical_ramp_closure:OutputExists', ...
        'Refusing to overwrite completed output: %s',final_file);
end
[library,library_provenance] = brusselator_build_ml_outcome_library(cfg);
[~,common_initial_state] = brusselator_initial_condition(cfg);
seed_ramps = (cfg.physical_ramp.left_ramp:cfg.physical_ramp.seed_spacing: ...
    cfg.physical_ramp.right_ramp)';

if isfile(checkpoint_file)
    saved = load(checkpoint_file,'records','raw_files','status');
    records = saved.records; raw_files = saved.raw_files; status = saved.status;
    fprintf('Resuming %s with %d completed protocol runs.\n',study_id,numel(records));
else
    if isfolder(outdir)
        error('run_physical_ramp_closure:UnsafeResume', ...
            'Existing output directory lacks a recognized checkpoint: %s',outdir);
    end
    mkdir(outdir); mkdir(fullfile(outdir,'raw_cases'));
    records = empty_record([]); raw_files = {}; status = 'seed_scan';
    brusselator_write_metadata(outdir,cfg, ...
        brusselator_run_metadata(cfg,mfilename,NaN));
    write_checkpoint();
end
run_timer = tic;

while numel(records) < cfg.physical_ramp.maximum_protocol_runs
    [ramp,source,parent] = next_protocol(records,seed_ramps,cfg);
    if ~isfinite(ramp)
        status = 'no_sampled_label_change_remaining_to_refine';
        break;
    end
    if strcmp(source,'seed_grid')
        status = 'seed_scan';
    else
        status = 'adaptive_refinement';
    end
    case_id = ramp_case_id(ramp);
    fprintf('Physical forcing run %d/%d: T_ramp=%.12g (%s).\n', ...
        numel(records)+1,cfg.physical_ramp.maximum_protocol_runs,ramp,source);
    raw_file = fullfile(outdir,'raw_cases',[case_id,'.mat']);
    partial_file = fullfile(outdir,'raw_cases',[case_id,'_partial.mat']);
    if isfile(raw_file)
        loaded = load(raw_file,'case_result'); case_result = loaded.case_result;
        record = case_result.record;
        fprintf('  Recovered completed case from %s.\n',raw_file);
    else
        [record,case_result] = execute_protocol(case_id,ramp,source,parent, ...
            common_initial_state,library,cfg,partial_file);
        record.raw_file = raw_file; case_result.record = record;
        save(raw_file,'case_result','-v7.3');
        if isfile(partial_file), delete(partial_file); end
    end
    records(end+1,1) = record; %#ok<AGROW>
    raw_files{end+1,1} = raw_file; %#ok<AGROW>
    write_checkpoint();
end
if numel(records) >= cfg.physical_ramp.maximum_protocol_runs
    status = 'fixed_protocol_budget_exhausted';
end

records = sort_records(records);
brackets = adjacent_brackets(records,cfg);
summary = summarize(records,brackets,status,cfg);
writetable(struct2table(records,'AsArray',true), ...
    fullfile(outdir,'physical_ramp_outcomes.csv'));
writetable(struct2table(brackets,'AsArray',true), ...
    fullfile(outdir,'adjacent_outcome_brackets.csv'));
writetable(struct2table(summary,'AsArray',true), ...
    fullfile(outdir,'physical_ramp_closure_summary.csv'));
make_figure(outdir,records);
results = struct('configuration',cfg,'library_provenance',library_provenance, ...
    'records',records,'raw_files',{raw_files},'brackets',brackets, ...
    'summary',summary,'status',status);
save(final_file,'results','-v7.3');
brusselator_write_metadata(outdir,cfg, ...
    brusselator_run_metadata(cfg,mfilename,toc(run_timer)));
disp(struct2table(records,'AsArray',true));
disp(struct2table(brackets,'AsArray',true));
disp(struct2table(summary,'AsArray',true));
fprintf('Physical ramp closure saved in: %s\n',outdir);

    function write_checkpoint()
        save(checkpoint_file,'cfg','records','raw_files','status','seed_ramps', ...
            'library_provenance','-v7.3');
        if ~isempty(records)
            writetable(struct2table(sort_records(records),'AsArray',true), ...
                fullfile(outdir,'progress_outcomes.csv'));
        end
    end
end

function [ramp,source,parent] = next_protocol(records,seed_ramps,cfg)
ramp = NaN; source = ''; parent = empty_parent();
tested = [records.ramp_time];
for k = 1:numel(seed_ramps)
    if isempty(tested) || all(abs(tested-seed_ramps(k))>1e-12)
        ramp = seed_ramps(k); source = 'seed_grid'; return;
    end
end
ordered = sort_records(records);
candidate = struct('width',{},'left_ramp',{},'right_ramp',{}, ...
    'left_outcome',{},'right_outcome',{});
for k = 1:(numel(ordered)-1)
    width = ordered(k+1).ramp_time-ordered(k).ramp_time;
    if ~strcmp(ordered(k).outcome,ordered(k+1).outcome) && ...
            width>cfg.physical_ramp.minimum_refinement_width
        midpoint = 0.5*(ordered(k).ramp_time+ordered(k+1).ramp_time);
        if all(abs(tested-midpoint)>1e-12)
            item = struct('width',width,'left_ramp',ordered(k).ramp_time, ...
                'right_ramp',ordered(k+1).ramp_time, ...
                'left_outcome',ordered(k).outcome, ...
                'right_outcome',ordered(k+1).outcome);
            candidate(end+1,1) = item; %#ok<AGROW>
        end
    end
end
if isempty(candidate), return; end
keys = [[candidate.width]',-[candidate.left_ramp]'];
[~,order] = sortrows(keys,[-1,2]); chosen = candidate(order(1));
ramp = 0.5*(chosen.left_ramp+chosen.right_ramp);
source = 'adaptive_midpoint'; parent = chosen;
end

function [record,result] = execute_protocol(case_id,ramp,source,parent, ...
        initial_state,library,cfg,partial_file)
protocol = struct('id',case_id,'bmax',cfg.physical_ramp.bmax, ...
    'Tup',ramp,'Thold',cfg.physical_ramp.Thold,'Tdown',ramp);
forcing_end = protocol.Tup+protocol.Thold+protocol.Tdown;
forcing_cfg = cfg; forcing_cfg.forcing = struct('b0',cfg.forcing.b0, ...
    'bmax',protocol.bmax,'Tup',protocol.Tup,'Thold',protocol.Thold, ...
    'Tdown',protocol.Tdown);
forcing_cfg.time.forcing_end = forcing_end;
forcing_cfg.time.Tfinal = forcing_end;
par = brusselator_make_parameters(forcing_cfg,@(t)overshoot_B(t, ...
    forcing_cfg.forcing.b0,protocol.bmax,protocol.Tup,protocol.Thold,protocol.Tdown));
timer = tic;
[landing_state,S,V,U] = solve_brusselator_1d_forced(initial_state,par,forcing_end,0);
landing_distances = state_distances(landing_state,library,cfg.grid.Lx);
frozen = brusselator_evolve_and_classify_state( ...
    [case_id,'_frozen'],landing_state,library,cfg,partial_file);
confirmation = struct(); confirmed_C = false;
if strcmp(frozen.classification.outcome,'C')
    fprintf('  Forcing-generated C candidate: running tighter confirmation.\n');
    confirm_cfg = cfg; confirm_cfg.solver = cfg.physical_ramp.confirmation_solver;
    confirm_cfg.forcing = forcing_cfg.forcing;
    confirm_par = brusselator_make_parameters(confirm_cfg,@(t)overshoot_B(t, ...
        confirm_cfg.forcing.b0,protocol.bmax,protocol.Tup,protocol.Thold,protocol.Tdown));
    [confirm_landing,CS,CV,CU] = solve_brusselator_1d_forced( ...
        initial_state,confirm_par,forcing_end,0);
    confirm_partial = strrep(partial_file,'.mat','_confirmation.mat');
    confirmation = brusselator_evolve_and_classify_state( ...
        [case_id,'_tight_C_confirmation'],confirm_landing,library, ...
        confirm_cfg,confirm_partial);
    confirmation.forcing_S = CS; confirmation.forcing_U = CU; ...
        confirmation.forcing_V = CV;
    if isfile(confirm_partial), delete(confirm_partial); end
    confirmed_C = strcmp(confirmation.classification.outcome,'C');
end
runtime = toc(timer);
record = empty_record(); record.protocol_id = case_id; record.source = source;
record.ramp_time = ramp; record.bmax = protocol.bmax; record.Thold = protocol.Thold;
record.forcing_end = forcing_end; record.outcome = frozen.classification.outcome;
record.decision_time = frozen.decision_time; record.tested_frozen_duration = frozen.duration;
record.landing_distance_A = landing_distances(1);
record.landing_distance_B = landing_distances(2);
record.landing_distance_C = landing_distances(3);
record.landing_distance_R5 = landing_distances(4);
record.landing_distance_R6 = landing_distances(5);
record.late_median_distance_A = frozen.classification.late_median_distances(1);
record.late_median_distance_B = frozen.classification.late_median_distances(2);
record.late_median_distance_C = frozen.classification.late_median_distances(3);
record.late_median_distance_R5 = frozen.classification.late_median_distances(4);
record.late_median_distance_R6 = frozen.classification.late_median_distances(5);
record.confirmation_outcome = 'not_run'; record.confirmed_C = confirmed_C;
if ~isempty(fieldnames(confirmation))
    record.confirmation_outcome = confirmation.classification.outcome;
end
record.parent_left_ramp = parent.left_ramp; record.parent_right_ramp = parent.right_ramp;
record.parent_left_outcome = parent.left_outcome;
record.parent_right_outcome = parent.right_outcome;
record.runtime_seconds = runtime;
result = struct('protocol',protocol,'initial_state',initial_state, ...
    'forcing_S',S,'forcing_U',U,'forcing_V',V,'landing_state',landing_state, ...
    'landing_distances',landing_distances,'frozen',frozen, ...
    'tight_C_confirmation',confirmation,'record',record);
fprintf('  outcome=%s, frozen t=%.0f, landing d=[%s] (%.1f s).\n', ...
    record.outcome,record.tested_frozen_duration, ...
    sprintf(' %.3e',landing_distances),runtime);
end

function distances = state_distances(state,library,Lx)
distances = zeros(1,numel(library));
for k = 1:numel(library)
    item = brusselator_physical_orbit_distance(state(:,2)',state(:,3)', ...
        library(k).U,library(k).V,Lx);
    distances(k) = item.distance;
end
end

function record = empty_record(varargin)
record = struct('protocol_id','','source','','ramp_time',NaN,'bmax',NaN, ...
    'Thold',NaN,'forcing_end',NaN,'outcome','unclassified', ...
    'decision_time',NaN,'tested_frozen_duration',NaN, ...
    'landing_distance_A',NaN,'landing_distance_B',NaN,'landing_distance_C',NaN, ...
    'landing_distance_R5',NaN,'landing_distance_R6',NaN, ...
    'late_median_distance_A',NaN,'late_median_distance_B',NaN, ...
    'late_median_distance_C',NaN,'late_median_distance_R5',NaN, ...
    'late_median_distance_R6',NaN,'confirmation_outcome','not_run', ...
    'confirmed_C',false,'parent_left_ramp',NaN,'parent_right_ramp',NaN, ...
    'parent_left_outcome','','parent_right_outcome','', ...
    'runtime_seconds',NaN,'raw_file','');
if nargin>0, record = record([]); end
end

function parent = empty_parent()
parent = struct('width',NaN,'left_ramp',NaN,'right_ramp',NaN, ...
    'left_outcome','','right_outcome','');
end

function records = sort_records(records)
if isempty(records), return; end
[~,order] = sort([records.ramp_time]); records = records(order);
end

function brackets = adjacent_brackets(records,cfg)
template = struct('bracket_id','','left_ramp',NaN,'left_outcome','', ...
    'right_ramp',NaN,'right_outcome','','width',NaN, ...
    'left_source_file','','right_source_file','', ...
    'eligible_for_dynamic_tracking',false,'contains_confirmed_C',false, ...
    'interpretation','adjacent different labels in the bounded physical scan');
brackets = template([]); known = cfg.ml_labeling.raw_library_order;
records = sort_records(records);
for k = 1:(numel(records)-1)
    if strcmp(records(k).outcome,records(k+1).outcome), continue; end
    item = template;
    item.bracket_id = sprintf('physical_%s_%s_%s_%s', ...
        records(k).outcome,records(k+1).outcome, ...
        numeric_tag(records(k).ramp_time),numeric_tag(records(k+1).ramp_time));
    item.left_ramp = records(k).ramp_time; item.left_outcome = records(k).outcome;
    item.right_ramp = records(k+1).ramp_time; item.right_outcome = records(k+1).outcome;
    item.width = item.right_ramp-item.left_ramp;
    item.left_source_file = records(k).raw_file;
    item.right_source_file = records(k+1).raw_file;
    item.eligible_for_dynamic_tracking = ismember(item.left_outcome,known) && ...
        ismember(item.right_outcome,known);
    item.contains_confirmed_C = (strcmp(item.left_outcome,'C') && records(k).confirmed_C) || ...
        (strcmp(item.right_outcome,'C') && records(k+1).confirmed_C);
    brackets(end+1,1) = item; %#ok<AGROW>
end
end

function summary = summarize(records,brackets,status,cfg)
labels = cfg.ml_labeling.raw_library_order;
counts = zeros(size(labels));
for k = 1:numel(labels), counts(k) = sum(strcmp({records.outcome},labels{k})); end
summary = struct('status',status,'protocol_runs',numel(records), ...
    'seed_runs',sum(strcmp({records.source},'seed_grid')), ...
    'adaptive_runs',sum(strcmp({records.source},'adaptive_midpoint')), ...
    'count_A',counts(1),'count_B',counts(2),'count_C',counts(3), ...
    'count_R5',counts(4),'count_R6',counts(5), ...
    'count_unresolved_or_ambiguous',numel(records)-sum(counts), ...
    'confirmed_forcing_C',sum([records.confirmed_C]), ...
    'adjacent_different_label_intervals',numel(brackets), ...
    'eligible_dynamic_brackets',sum([brackets.eligible_for_dynamic_tracking]), ...
    'minimum_different_label_width',NaN, ...
    'claim_scope',cfg.physical_ramp.claim_scope);
if ~isempty(brackets), summary.minimum_different_label_width = min([brackets.width]); end
end

function make_figure(outdir,records)
codes = containers.Map({'A','B','C','R5','R6','unresolved', ...
    'ambiguous_multiple_outcomes'},{1,2,3,4,5,6,7});
y = zeros(size(records));
for k = 1:numel(records)
    if isKey(codes,records(k).outcome), y(k) = codes(records(k).outcome); else, y(k)=6; end
end
fig = figure('Visible','off','Color','w','Position',[100 100 980 430]);
scatter([records.ramp_time],y,65,y,'filled');
yticks(1:7); yticklabels({'A','B','C','R5','R6','unresolved','ambiguous'});
xlabel('Symmetric ramp time'); ylabel('Frozen outcome');
title('Bounded physical forcing refinement, b_{max}=11.14, T_{hold}=80');
grid on; exportgraphics(fig,fullfile(outdir,'physical_ramp_outcome_map.png'), ...
    'Resolution',300); close(fig);
end

function id = ramp_case_id(ramp)
id = ['ramp_',numeric_tag(ramp)];
end

function tag = numeric_tag(value)
tag = strrep(sprintf('%.12g',value),'.','p'); tag = strrep(tag,'-','m');
end
