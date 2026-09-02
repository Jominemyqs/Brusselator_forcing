function results = run_ml_BC_subset_labeling(cfg)
%RUN_ML_BC_SUBSET_LABELING Adaptive frozen-PDE labels for a frozen subset.

if nargin<1, cfg=brusselator_ml_acquisition_labeling_config(); end
validate_ml_BC_boundary_labeling_setup(cfg);
if ~isfile(cfg.ml_subset.prerequisite_file)
    error('run_ml_BC_subset_labeling:MissingPrerequisite', ...
        'Required frozen prerequisite is missing: %s',cfg.ml_subset.prerequisite_file);
end
if cfg.ml_subset.require_frozen_model_manifest
    validate_frozen_model_manifest(cfg.ml_subset.prerequisite_file);
elseif isfield(cfg.ml_subset,'require_frozen_pair_manifest') && ...
        cfg.ml_subset.require_frozen_pair_manifest
    validate_frozen_pair_manifest(cfg.ml_subset.prerequisite_file);
elseif isfield(cfg.ml_subset,'require_frozen_matched_protocol') && ...
        cfg.ml_subset.require_frozen_matched_protocol
    validate_frozen_matched_protocol(cfg.ml_subset.prerequisite_file,cfg);
end
if isfield(cfg.ml_subset,'require_frozen_shard_plan') && ...
        cfg.ml_subset.require_frozen_shard_plan
    validate_frozen_shard_plan(cfg.ml_subset.shard_plan_file,cfg);
end
timer=tic; outdir=fullfile(cfg.output.root,cfg.experiment_name);
rawdir=fullfile(outdir,'raw_cases'); partialdir=fullfile(outdir,'partial_cases');
checkpoint=fullfile(outdir,'subset_labeling_checkpoint.mat');
final_file=fullfile(outdir,cfg.ml_subset.final_filename);
if isfile(final_file)
    error('run_ml_BC_subset_labeling:OutputExists', ...
        'Refusing to overwrite completed subset labels: %s',final_file);
end
loaded=load(cfg.ml_labeling.geometry_file,'results'); geometry=loaded.results;
[library,library_provenance]=brusselator_build_ml_outcome_library(cfg);
selection=select_subset(geometry,cfg);
selected=selection.state_index;
geometry_id=geometry.summary.geometry_id;

if isfile(checkpoint)
    saved=load(checkpoint,'records','raw_files','next_position', ...
        'saved_cfg','geometry_id_saved');
    if ~execution_compatible(saved.saved_cfg,cfg) || ...
            ~strcmp(saved.geometry_id_saved,geometry_id)
        error('run_ml_BC_subset_labeling:CheckpointMismatch', ...
            'Checkpoint configuration or geometry differs.');
    end
    records=saved.records; raw_files=saved.raw_files; next_position=saved.next_position;
    fprintf('Resuming %s at case %d/%d.\n',cfg.ml_subset.kind,next_position,numel(selected));
else
    if isfolder(outdir)
        error('run_ml_BC_subset_labeling:IncompleteOutput', ...
            'Existing output has no valid checkpoint: %s',outdir);
    end
    mkdir(outdir); mkdir(rawdir); mkdir(partialdir);
    records=empty_record([]); raw_files={}; next_position=1;
    brusselator_write_metadata(outdir,cfg,brusselator_run_metadata(cfg,mfilename,NaN));
    writetable(struct2table(library_provenance,'AsArray',true), ...
        fullfile(outdir,'outcome_library_provenance.csv'));
    writetable(selection,fullfile(outdir,'frozen_selection_manifest.csv'));
    save_checkpoint();
end

processed=0;
parallel_workers=1;parallel_batch_size=1;
if isfield(cfg.ml_subset,'parallel_workers')
    parallel_workers=cfg.ml_subset.parallel_workers;
    parallel_batch_size=cfg.ml_subset.parallel_batch_size;
end
if parallel_workers>1 && exist('gcp','file')~=2
    warning('run_ml_BC_subset_labeling:ParallelRuntimeUnavailable', ...
        ['Parallel execution was requested, but gcp/parpool are not installed ', ...
        'in this MATLAB runtime. Falling back to the unchanged sequential path.']);
    parallel_workers=1;parallel_batch_size=1;
end
if parallel_workers>1
    pool=gcp('nocreate');
    if isempty(pool)||pool.NumWorkers~=parallel_workers
        if ~isempty(pool),delete(pool);end
        parpool('Processes',parallel_workers);
    end
end
while next_position<=numel(selected)
    remaining_budget=cfg.ml_subset.maximum_cases_per_invocation-processed;
    if remaining_budget<=0
        results=incomplete_result(); return;
    end
    count=min([parallel_batch_size,numel(selected)-next_position+1,remaining_budget]);
    positions=next_position:(next_position+count-1);
    batch_records=cell(count,1);batch_files=cell(count,1);
    if parallel_workers>1
        parfor local_index=1:count
            [batch_records{local_index},batch_files{local_index}]=process_subset_case( ...
                positions(local_index),selected,geometry,selection,library,cfg, ...
                rawdir,partialdir);
        end
    else
        for local_index=1:count
            [batch_records{local_index},batch_files{local_index}]=process_subset_case( ...
                positions(local_index),selected,geometry,selection,library,cfg, ...
                rawdir,partialdir);
        end
    end
    for local_index=1:count
        records(end+1,1)=batch_records{local_index}; %#ok<AGROW>
        raw_files{end+1,1}=batch_files{local_index}; %#ok<AGROW>
    end
    next_position=positions(end)+1;processed=processed+count;
    writetable(struct2table(records,'AsArray',true),fullfile(outdir,'labeling_progress.csv'));
    save_checkpoint();
end

records=sort_records(records); summary=summarize_study(records,cfg);
writetable(struct2table(records,'AsArray',true),fullfile(outdir,cfg.ml_subset.final_csv));
writetable(struct2table(summary,'AsArray',true),fullfile(outdir,cfg.ml_subset.summary_csv));
results=struct('configuration',cfg,'geometry_id',geometry_id, ...
    'selection',selection,'library',library,'library_provenance',library_provenance, ...
    'records',records,'raw_files',{raw_files},'summary',summary);
save(final_file,'results','-v7.3'); make_figure(outdir,records,cfg);
brusselator_write_metadata(outdir,cfg,brusselator_run_metadata(cfg,mfilename,toc(timer)));
disp(struct2table(summary,'AsArray',true));
fprintf('%s labels saved in: %s\n',cfg.ml_subset.kind,outdir);

    function save_checkpoint()
        saved_cfg=cfg; geometry_id_saved=geometry_id;
        save(checkpoint,'records','raw_files','next_position','saved_cfg', ...
            'geometry_id_saved','-v7.3');
    end
    function value=incomplete_result()
        value=struct('status','incomplete_checkpointed','records',records, ...
            'next_position',next_position,'selected_count',numel(selected));
        brusselator_write_metadata(outdir,cfg, ...
            brusselator_run_metadata(cfg,mfilename,toc(timer)));
    end
end

function compatible=execution_compatible(saved_cfg,current_cfg)
fields={'parallel_workers','parallel_batch_size'};
for k=1:numel(fields)
    if isfield(saved_cfg.ml_subset,fields{k})
        saved_cfg.ml_subset=rmfield(saved_cfg.ml_subset,fields{k});
    end
    if isfield(current_cfg.ml_subset,fields{k})
        current_cfg.ml_subset=rmfield(current_cfg.ml_subset,fields{k});
    end
end
compatible=isequaln(saved_cfg,current_cfg);
end

function [record,raw_file]=process_subset_case( ...
        position,selected,geometry,selection,library,cfg,rawdir,partialdir)
state_index=selected(position);row=geometry.manifest(state_index);
annotation=selection(position,:);case_id=row.sample_id;
initial_state=[geometry.center.x,geometry.states.raw_U(state_index,:)', ...
    geometry.states.raw_V(state_index,:)'];
raw_file=fullfile(rawdir,[case_id,'.mat']);
partial_file=fullfile(partialdir,[case_id,'_partial.mat']);
fprintf('\n%s %d/%d: %s at (%.6g,%.6g).\n',cfg.ml_subset.kind, ...
    position,numel(selected),case_id,row.alpha_1,row.alpha_2);
if isfile(raw_file)
    item=load(raw_file,'case_data');case_data=item.case_data;
else
    case_data=brusselator_evolve_and_classify_state( ...
        case_id,initial_state,library,cfg,partial_file);
    manifest_row=row;selection_row=annotation;
    save(raw_file,'case_data','manifest_row','selection_row','-v7.3');
    if isfile(partial_file),delete(partial_file);end
end
record=summarize_case(case_data,row,annotation,library,cfg,raw_file);
end

function selection=select_subset(geometry,cfg)
switch cfg.ml_subset.selection_mode
    case 'proposal_union'
        A=readtable(cfg.ml_subset.active_proposals_file);
        B=readtable(cfg.ml_subset.baseline_proposals_file);
        if height(A)~=cfg.ml_subset.required_active_count || ...
                height(B)~=cfg.ml_subset.required_baseline_count
            error('run_ml_BC_subset_labeling:ProposalCount','Proposal count changed.');
        end
        overlap=intersect(A.state_index,B.state_index);
        if numel(overlap)~=cfg.ml_subset.required_overlap_count
            error('run_ml_BC_subset_labeling:ProposalOverlap','Proposal overlap changed.');
        end
        indices=sort(union(A.state_index,B.state_index));
        selection=table(indices,'VariableNames',{'state_index'});
        selection.selected_by_active=ismember(indices,A.state_index);
        selection.selected_by_baseline=ismember(indices,B.state_index);
        selection.active_rank=nan(size(indices)); selection.baseline_rank=nan(size(indices));
        [~,ia]=ismember(indices,A.state_index); mask=ia>0;
        selection.active_rank(mask)=A.acquisition_step(ia(mask));
        [~,ib]=ismember(indices,B.state_index); mask=ib>0;
        selection.baseline_rank(mask)=B.baseline_step(ib(mask));
    case 'design_role'
        indices=find(strcmp({geometry.manifest.design_role},cfg.ml_subset.selected_design_role));
        selection=table(indices(:),'VariableNames',{'state_index'});
        selection.selected_by_active=false(size(indices(:)));
        selection.selected_by_baseline=false(size(indices(:)));
        selection.active_rank=nan(size(indices(:))); selection.baseline_rank=nan(size(indices(:)));
    case 'edge_bracket_union'
        A=readtable(cfg.ml_subset.active_proposals_file);
        B=readtable(cfg.ml_subset.baseline_proposals_file);
        active_endpoints=unique([A.B_state_index;A.C_state_index]);
        baseline_endpoints=unique([B.B_state_index;B.C_state_index]);
        indices=sort(union(active_endpoints,baseline_endpoints));
        selection=table(indices,'VariableNames',{'state_index'});
        selection.selected_by_active=ismember(indices,active_endpoints);
        selection.selected_by_baseline=ismember(indices,baseline_endpoints);
        selection.active_rank=nan(size(indices)); selection.baseline_rank=nan(size(indices));
        for k=1:height(A)
            mask=ismember(indices,[A.B_state_index(k),A.C_state_index(k)]);
            selection.active_rank(mask)=A.proposal_rank(k);
        end
        for k=1:height(B)
            mask=ismember(indices,[B.B_state_index(k),B.C_state_index(k)]);
            selection.baseline_rank(mask)=B.proposal_rank(k);
        end
    case 'pair_batch_union'
        P=readtable(cfg.ml_subset.pair_batch_file);
        indices=sort(unique([P.first_state_index;P.second_state_index]));
        selection=table(indices,'VariableNames',{'state_index'});
        selection.selected_by_active=false(size(indices));
        selection.selected_by_baseline=false(size(indices));
        selection.active_rank=nan(size(indices)); selection.baseline_rank=nan(size(indices));
        for k=1:height(P)
            mask=ismember(indices,[P.first_state_index(k),P.second_state_index(k)]);
            selection.active_rank(mask)=P.batch_step(k);
        end
    case 'state_index_manifest'
        P=readtable(cfg.ml_subset.state_index_file);
        if ~ismember('state_index',P.Properties.VariableNames)
            error('run_ml_BC_subset_labeling:StateIndexManifest', ...
                'The frozen manifest must contain state_index.');
        end
        indices=sort(unique(P.state_index));
        selection=table(indices,'VariableNames',{'state_index'});
        selection.selected_by_active=false(size(indices));
        selection.selected_by_baseline=false(size(indices));
        selection.active_rank=nan(size(indices));
        selection.baseline_rank=nan(size(indices));
    otherwise
        error('run_ml_BC_subset_labeling:SelectionMode','Unknown selection mode.');
end
if height(selection)~=cfg.ml_subset.expected_unique_count
    error('run_ml_BC_subset_labeling:UniqueCount','Selected-state count changed.');
end
roles={geometry.manifest(selection.state_index).design_role};
if isfield(cfg.ml_subset,'required_design_roles')
    allowed=cfg.ml_subset.required_design_roles;
else
    allowed={cfg.ml_subset.required_design_role};
end
if ~all(ismember(roles,allowed))
    error('run_ml_BC_subset_labeling:RoleLeakage','Selection includes a forbidden role.');
end
end

function validate_frozen_model_manifest(file)
data=jsondecode(fileread(file));
required={'models_frozen','fixed_test_labels_used','active_model_hash','baseline_model_hash'};
has_required=all(isfield(data,required));
has_hashes=has_required && strlength(string(data.active_model_hash))==64 && ...
    strlength(string(data.baseline_model_hash))==64;
if ~has_required || ~data.models_frozen || data.fixed_test_labels_used || ~has_hashes
    error('run_ml_BC_subset_labeling:ModelsNotFrozen', ...
        'Blind test labels are forbidden until both comparison models are frozen.');
end
end

function validate_frozen_pair_manifest(file)
data=jsondecode(fileread(file));
required={'pair_batch_frozen','endpoint_outcomes_used','endpoint_count', ...
    'frozen_pair_batch_hash','frozen_endpoint_manifest_hash'};
has_required=all(isfield(data,required));
has_hashes=has_required && strlength(string(data.frozen_pair_batch_hash))==64 && ...
    strlength(string(data.frozen_endpoint_manifest_hash))==64;
if ~has_required || ~data.pair_batch_frozen || data.endpoint_outcomes_used || ...
        data.endpoint_count~=8 || ~has_hashes
    error('run_ml_BC_subset_labeling:PairBatchNotFrozen', ...
        'Pair-aware endpoint labels require an outcome-blind frozen eight-state batch.');
end
end

function validate_frozen_matched_protocol(file,cfg)
data=jsondecode(fileread(file));
required={'protocol_frozen','replicate_count','training_count_per_replicate', ...
    'pair_count_per_method','endpoint_budget_per_method', ...
    'missing_label_state_count','labels_used_for_design', ...
    'candidate_outcomes_used','matched_initial_designs_sha256', ...
    'missing_label_manifest_sha256'};
has_required=all(isfield(data,required));
has_hashes=has_required && ...
    strlength(string(data.matched_initial_designs_sha256))==64 && ...
    strlength(string(data.missing_label_manifest_sha256))==64;
if ~has_required || ~data.protocol_frozen || data.replicate_count~=12 || ...
        data.training_count_per_replicate~=25 || ...
        data.pair_count_per_method~=4 || data.endpoint_budget_per_method~=8 || ...
        data.missing_label_state_count~=parent_missing_count(cfg) || ...
        data.labels_used_for_design || data.candidate_outcomes_used || ~has_hashes
    error('run_ml_BC_subset_labeling:MatchedProtocolNotFrozen', ...
        'Matched label-bank generation requires the frozen outcome-blind protocol.');
end
end

function count=parent_missing_count(cfg)
if isfield(cfg.ml_subset,'parent_missing_label_state_count')
    count=cfg.ml_subset.parent_missing_label_state_count;
else
    count=cfg.ml_subset.expected_unique_count;
end
end

function validate_frozen_shard_plan(file,cfg)
data=jsondecode(fileread(file));
if ~isfield(data,'shards_frozen')||~data.shards_frozen|| ...
        ~isfield(data,'shard_count')||data.shard_count~=4|| ...
        ~isfield(data,'outcomes_used_for_sharding')||data.outcomes_used_for_sharding|| ...
        ~isfield(data,'parent_missing_label_state_count')|| ...
        data.parent_missing_label_state_count~=parent_missing_count(cfg)
    error('run_ml_BC_subset_labeling:ShardPlanNotFrozen', ...
        'Label-bank shard execution requires the frozen outcome-blind plan.');
end
end

function record=summarize_case(data,row,annotation,library,cfg,raw_file)
raw=data.classification.outcome;
if ismember(raw,{'B','C'})
    label=raw; status='resolved_primary'; other='';
elseif ismember(raw,cfg.ml_labeling.other_raw_outcomes)
    label='other'; status='resolved_other'; other=raw;
else
    label='U'; other=''; status='unresolved_or_ambiguous_at_maximum_horizon';
end
labels={library.label}; order=cellfun(@(x)find(strcmp(labels,x),1), ...
    cfg.ml_labeling.raw_library_order);
med=data.classification.late_median_distances(order);
mx=data.classification.late_maximum_distances(order);
record=empty_record(); record.state_index=row.state_index; record.sample_id=row.sample_id;
record.alpha_1=row.alpha_1; record.alpha_2=row.alpha_2;
record.normalized_alpha_1=row.normalized_alpha_1;
record.normalized_alpha_2=row.normalized_alpha_2;
record.selected_by_active=annotation.selected_by_active;
record.selected_by_baseline=annotation.selected_by_baseline;
record.active_rank=annotation.active_rank; record.baseline_rank=annotation.baseline_rank;
record.final_label=label; record.label_status=status;
record.raw_five_way_outcome=raw; record.other_known_outcome=other;
record.tested_duration=data.duration; record.decision_time=data.decision_time;
record.runtime_seconds=data.runtime_seconds;
record.late_median_A=med(1); record.late_median_B=med(2); record.late_median_C=med(3);
record.late_median_R5=med(4); record.late_median_R6=med(5);
record.late_maximum_A=mx(1); record.late_maximum_B=mx(2); record.late_maximum_C=mx(3);
record.late_maximum_R5=mx(4); record.late_maximum_R6=mx(5); record.raw_file=raw_file;
end

function record=empty_record(varargin)
record=struct('state_index',NaN,'sample_id','','alpha_1',NaN,'alpha_2',NaN, ...
    'normalized_alpha_1',NaN,'normalized_alpha_2',NaN, ...
    'selected_by_active',false,'selected_by_baseline',false, ...
    'active_rank',NaN,'baseline_rank',NaN,'final_label','','label_status','', ...
    'raw_five_way_outcome','','other_known_outcome','','tested_duration',NaN, ...
    'decision_time',NaN,'runtime_seconds',NaN,'late_median_A',NaN, ...
    'late_median_B',NaN,'late_median_C',NaN,'late_median_R5',NaN, ...
    'late_median_R6',NaN,'late_maximum_A',NaN,'late_maximum_B',NaN, ...
    'late_maximum_C',NaN,'late_maximum_R5',NaN,'late_maximum_R6',NaN,'raw_file','');
if nargin>0, record=record([]); end
end

function records=sort_records(records)
[~,order]=sort([records.state_index]); records=records(order);
end

function summary=summarize_study(records,cfg)
labels={records.final_label};
summary=struct('subset_kind',cfg.ml_subset.kind,'case_count',numel(records), ...
    'B_count',sum(strcmp(labels,'B')),'C_count',sum(strcmp(labels,'C')), ...
    'other_count',sum(strcmp(labels,'other')),'U_count',sum(strcmp(labels,'U')), ...
    'active_budget_count',sum([records.selected_by_active]), ...
    'baseline_budget_count',sum([records.selected_by_baseline]), ...
    'maximum_tested_duration',max([records.tested_duration]), ...
    'total_runtime_seconds',sum([records.runtime_seconds]), ...
    'fixed_test_outcomes_touched',strcmp(cfg.ml_subset.kind,'blind_fixed_test'), ...
    'interpretation',cfg.ml_subset.interpretation);
end

function make_figure(outdir,records,cfg)
fig=figure('Visible','off','Color','w'); hold on;
for item={{'B',[0 .45 .74]},{'C',[.85 .33 .1]},{'other',[.49 .18 .56]},{'U',[.3 .3 .3]}}
    selected=strcmp({records.final_label},item{1}{1});
    scatter([records(selected).alpha_1],[records(selected).alpha_2],70, ...
        item{1}{2},'filled','DisplayName',item{1}{1});
end
xlabel('\alpha_1'); ylabel('\alpha_2'); axis equal; grid on; legend('Location','best');
title(strrep(cfg.ml_subset.kind,'_',' '));
exportgraphics(fig,fullfile(outdir,cfg.ml_subset.figure_filename),'Resolution',300); close(fig);
end
