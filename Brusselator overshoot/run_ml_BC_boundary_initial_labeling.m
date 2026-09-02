function results = run_ml_BC_boundary_initial_labeling(cfg)
%RUN_ML_BC_BOUNDARY_INITIAL_LABELING Adaptively label the frozen 25-point seed.
%   Each state is evolved under the actual frozen b=10 PDE. Unique B and C
%   matches retain their labels; A/R5/R6 map to other; unresolved or ambiguous
%   cases map to U. Per-case partial checkpoints allow long trajectories to be
%   resumed without changing the scientific classification definitions.

if nargin<1, cfg=brusselator_ml_boundary_labeling_config(); end
validate_ml_BC_boundary_labeling_setup(cfg);
run_timer=tic;
outdir=fullfile(cfg.output.root,cfg.experiment_name);
rawdir=fullfile(outdir,'raw_cases'); partialdir=fullfile(outdir,'partial_cases');
checkpoint=fullfile(outdir,'labeling_checkpoint.mat');
final_file=fullfile(outdir,'ml_BC_boundary_initial_labels.mat');
if isfile(final_file)
    error('run_ml_BC_boundary_initial_labeling:OutputExists', ...
        'Refusing to overwrite completed labels: %s',final_file);
end

geometry_loaded=load(cfg.ml_labeling.geometry_file,'results');
geometry=geometry_loaded.results;
[library,library_provenance]=brusselator_build_ml_outcome_library(cfg);
selected=find(strcmp({geometry.manifest.design_role}, ...
    cfg.ml_labeling.selected_design_role));
selected=sort(selected(:));
geometry_id=geometry.summary.geometry_id;

if isfile(checkpoint)
    saved=load(checkpoint,'records','raw_files','next_position', ...
        'saved_cfg','geometry_id_saved');
    if ~isequaln(saved.saved_cfg,cfg) || ~strcmp(saved.geometry_id_saved,geometry_id)
        error('run_ml_BC_boundary_initial_labeling:CheckpointMismatch', ...
            'Existing checkpoint uses another configuration or geometry.');
    end
    records=saved.records; raw_files=saved.raw_files;
    next_position=saved.next_position;
    fprintf('Resuming adaptive ML-plane labeling at seed %d of %d.\n', ...
        next_position,numel(selected));
else
    if isfolder(outdir)
        error('run_ml_BC_boundary_initial_labeling:IncompleteOutput', ...
            'Existing output directory has no valid checkpoint: %s',outdir);
    end
    mkdir(outdir); mkdir(rawdir); mkdir(partialdir);
    records=empty_record([]); raw_files={}; next_position=1;
    brusselator_write_metadata(outdir,cfg, ...
        brusselator_run_metadata(cfg,mfilename,NaN));
    write_library_provenance(outdir,library_provenance,cfg);
    save_checkpoint();
end

processed=0;
for position=next_position:numel(selected)
    if processed>=cfg.ml_labeling.maximum_cases_per_invocation
        fprintf('Invocation limit reached after %d cases; checkpoint retained.\n',processed);
        results=incomplete_result(); return;
    end
    state_index=selected(position); row=geometry.manifest(state_index);
    case_id=row.sample_id;
    initial_state=[geometry.center.x, ...
        geometry.states.raw_U(state_index,:)', ...
        geometry.states.raw_V(state_index,:)'];
    raw_file=fullfile(rawdir,[case_id,'.mat']);
    partial_file=fullfile(partialdir,[case_id,'_partial.mat']);
    fprintf('\nSeed %d/%d: %s at (alpha_1,alpha_2)=(%.6g,%.6g).\n', ...
        position,numel(selected),case_id,row.alpha_1,row.alpha_2);
    if isfile(raw_file)
        loaded=load(raw_file,'case_data'); case_data=loaded.case_data;
    else
        case_data=brusselator_evolve_and_classify_state( ...
            case_id,initial_state,library,cfg,partial_file);
        manifest_row=row;
        save(raw_file,'case_data','manifest_row','-v7.3');
        if isfile(partial_file), delete(partial_file); end
    end
    record=summarize_case(case_data,row,library,cfg,raw_file);
    records(end+1,1)=record; %#ok<AGROW>
    raw_files{end+1,1}=raw_file; %#ok<AGROW>
    next_position=position+1; processed=processed+1;
    writetable(struct2table(records,'AsArray',true), ...
        fullfile(outdir,'labeling_progress.csv'));
    save_checkpoint();
end

records=sort_records(records);
summary=summarize_study(records,cfg);
writetable(struct2table(records,'AsArray',true), ...
    fullfile(outdir,'initial_seed_labels.csv'));
writetable(struct2table(summary,'AsArray',true), ...
    fullfile(outdir,'initial_labeling_summary.csv'));
results=struct('configuration',cfg,'geometry_id',geometry_id, ...
    'selected_state_indices',selected,'library',library, ...
    'library_provenance',library_provenance,'records',records, ...
    'raw_files',{raw_files},'summary',summary);
save(final_file,'results','-v7.3');
plot_ml_BC_boundary_initial_labels(final_file);
brusselator_write_metadata(outdir,cfg, ...
    brusselator_run_metadata(cfg,mfilename,toc(run_timer)));
disp(struct2table(records,'AsArray',true));
disp(struct2table(summary,'AsArray',true));
fprintf('Adaptive initial-seed labels saved in: %s\n',outdir);

    function save_checkpoint()
        saved_cfg=cfg; geometry_id_saved=geometry_id;
        save(checkpoint,'records','raw_files','next_position','saved_cfg', ...
            'geometry_id_saved','-v7.3');
    end
    function partial=incomplete_result()
        partial=struct('status','incomplete_checkpointed','configuration',cfg, ...
            'geometry_id',geometry_id,'records',records,'raw_files',{raw_files}, ...
            'next_position',next_position,'selected_count',numel(selected));
        brusselator_write_metadata(outdir,cfg, ...
            brusselator_run_metadata(cfg,mfilename,toc(run_timer)));
    end
end

function record = summarize_case(data,row,library,cfg,raw_file)
raw=data.classification.outcome;
if strcmp(raw,'B') || strcmp(raw,'C')
    label=raw; status='resolved_primary'; other='';
elseif ismember(raw,cfg.ml_labeling.other_raw_outcomes)
    label='other'; status='resolved_other'; other=raw;
else
    label='U'; other='';
    if strcmp(raw,'ambiguous_multiple_outcomes')
        status='ambiguous_at_maximum_horizon';
    else
        status='unresolved_at_maximum_horizon';
    end
end
labels={library.label};
indices=zeros(1,5);
for k=1:5, indices(k)=find(strcmp(labels,cfg.ml_labeling.raw_library_order{k}),1); end
med=data.classification.late_median_distances(indices);
mx=data.classification.late_maximum_distances(indices);
fin=data.classification.final_distances(indices);
minimum_history=min(data.history.distances(:,indices),[],1);
recurrence_period=NaN; recurrence_distance=NaN;
if data.recurrence.available
    recurrence_period=data.recurrence.period;
    recurrence_distance=data.recurrence.distance;
end
record=empty_record();
record.state_index=row.state_index; record.sample_id=row.sample_id;
record.row_index=row.row_index; record.column_index=row.column_index;
record.alpha_1=row.alpha_1; record.alpha_2=row.alpha_2;
record.normalized_alpha_1=row.normalized_alpha_1;
record.normalized_alpha_2=row.normalized_alpha_2;
record.final_label=label; record.label_status=status;
record.raw_five_way_outcome=raw; record.other_known_outcome=other;
record.tested_duration=data.duration; record.decision_time=data.decision_time;
record.runtime_seconds=data.runtime_seconds;
record.late_median_A=med(1); record.late_median_B=med(2);
record.late_median_C=med(3); record.late_median_R5=med(4);
record.late_median_R6=med(5); record.late_maximum_A=mx(1);
record.late_maximum_B=mx(2); record.late_maximum_C=mx(3);
record.late_maximum_R5=mx(4); record.late_maximum_R6=mx(5);
record.final_distance_A=fin(1); record.final_distance_B=fin(2);
record.final_distance_C=fin(3); record.final_distance_R5=fin(4);
record.final_distance_R6=fin(5); record.minimum_history_A=minimum_history(1);
record.minimum_history_B=minimum_history(2); record.minimum_history_C=minimum_history(3);
record.minimum_history_R5=minimum_history(4); record.minimum_history_R6=minimum_history(5);
record.recurrence_available=data.recurrence.available;
record.recurrence_period=recurrence_period;
record.recurrence_distance=recurrence_distance; record.raw_file=raw_file;
end

function record = empty_record(varargin)
record=struct('state_index',NaN,'sample_id','','row_index',NaN, ...
    'column_index',NaN,'alpha_1',NaN,'alpha_2',NaN, ...
    'normalized_alpha_1',NaN,'normalized_alpha_2',NaN, ...
    'final_label','','label_status','','raw_five_way_outcome','', ...
    'other_known_outcome','','tested_duration',NaN,'decision_time',NaN, ...
    'runtime_seconds',NaN,'late_median_A',NaN,'late_median_B',NaN, ...
    'late_median_C',NaN,'late_median_R5',NaN,'late_median_R6',NaN, ...
    'late_maximum_A',NaN,'late_maximum_B',NaN,'late_maximum_C',NaN, ...
    'late_maximum_R5',NaN,'late_maximum_R6',NaN,'final_distance_A',NaN, ...
    'final_distance_B',NaN,'final_distance_C',NaN,'final_distance_R5',NaN, ...
    'final_distance_R6',NaN,'minimum_history_A',NaN,'minimum_history_B',NaN, ...
    'minimum_history_C',NaN,'minimum_history_R5',NaN,'minimum_history_R6',NaN, ...
    'recurrence_available',false,'recurrence_period',NaN, ...
    'recurrence_distance',NaN,'raw_file','');
if nargin>0, record=record([]); end
end

function records = sort_records(records)
[~,order]=sort([records.state_index]); records=records(order);
end

function summary = summarize_study(records,cfg)
labels={records.final_label};
summary=struct('geometry_id',cfg.ml_labeling.expected_geometry_id, ...
    'case_count',numel(records),'B_count',sum(strcmp(labels,'B')), ...
    'C_count',sum(strcmp(labels,'C')),'other_count',sum(strcmp(labels,'other')), ...
    'U_count',sum(strcmp(labels,'U')), ...
    'A_raw_count',sum(strcmp({records.raw_five_way_outcome},'A')), ...
    'R5_raw_count',sum(strcmp({records.raw_five_way_outcome},'R5')), ...
    'R6_raw_count',sum(strcmp({records.raw_five_way_outcome},'R6')), ...
    'maximum_tested_duration',max([records.tested_duration]), ...
    'total_runtime_seconds',sum([records.runtime_seconds]), ...
    'exact_edge_used_for_labeling',false, ...
    'held_out_CA_region_used_for_labeling',false, ...
    'interpretation',cfg.ml_labeling.claim_scope);
end

function write_library_provenance(outdir,provenance,cfg)
writetable(struct2table(provenance,'AsArray',true), ...
    fullfile(outdir,'outcome_library_provenance.csv'));
mapping=table(cfg.ml_labeling.raw_library_order(:), ...
    {'other';'B';'C';'other';'other'}, ...
    'VariableNames',{'raw_five_way_outcome','reported_label'});
mapping=[mapping;table({'unresolved';'ambiguous_multiple_outcomes'}, {'U';'U'}, ...
    'VariableNames',mapping.Properties.VariableNames)];
writetable(mapping,fullfile(outdir,'label_mapping.csv'));
end
