function results = run_ml_heldout_CBA_round_labeling(round_index)
%RUN_ML_HELDOUT_CBA_ROUND_LABELING Label one prospectively frozen pair.

cfg=brusselator_ml_heldout_CBA_labeling_config(round_index);
round_dir=fileparts(cfg.heldout_round.selection_file);
final_file=fullfile(round_dir,'endpoint_labels.mat');
final_csv=fullfile(round_dir,'endpoint_labels.csv');
if isfile(final_file)||isfile(final_csv)
    error('run_ml_heldout_CBA_round_labeling:OutputExists', ...
        'Refusing to overwrite completed round %d labels.',round_index);
end
if ~isfile(cfg.heldout_round.selection_file)|| ...
        ~isfile(cfg.heldout_round.selection_manifest)
    error('run_ml_heldout_CBA_round_labeling:SelectionMissing', ...
        'The prospectively frozen round selection is missing.');
end
manifest=jsondecode(fileread(cfg.heldout_round.selection_manifest));
if ~manifest.selection_frozen_before_round_labels || manifest.labels_opened
    error('run_ml_heldout_CBA_round_labeling:NotFrozen', ...
        'Round selection was not frozen in an unopened state.');
end
selection=readtable(cfg.heldout_round.selection_file,'TextType','string');
indices=[selection.first_state_index;selection.second_state_index];
if height(selection)~=1||numel(unique(indices))~=2|| ...
        ~isequal(indices(:),manifest.selected_state_indices(:))
    error('run_ml_heldout_CBA_round_labeling:SelectionMismatch', ...
        'Selection CSV and frozen manifest do not agree.');
end
loaded=load(cfg.ml_labeling.geometry_file,'results');geometry=loaded.results;
[library,library_provenance]=brusselator_build_ml_outcome_library(cfg);
records=repmat(empty_record(),2,1);raw_files=cell(2,1);timer=tic;
rawdir=fullfile(round_dir,'raw_cases');partialdir=fullfile(round_dir,'partial_cases');
if ~isfolder(rawdir),mkdir(rawdir);end
if ~isfolder(partialdir),mkdir(partialdir);end
for position=1:2
    index=indices(position);row=geometry.candidate_manifest(index,:);
    case_id=char(row.sample_id);
    state=[geometry.states.x,geometry.states.raw_U(index,:)', ...
        geometry.states.raw_V(index,:)'];
    raw_file=fullfile(rawdir,[case_id,'.mat']);
    partial_file=fullfile(partialdir,[case_id,'_partial.mat']);
    fprintf('\nHeld-out round %d endpoint %d/2: %s, lambda=%.10f.\n', ...
        round_index,position,case_id,row.lambda);
    if isfile(raw_file)
        saved=load(raw_file,'case_data');case_data=saved.case_data;
    else
        case_data=brusselator_evolve_and_classify_state( ...
            case_id,state,library,cfg,partial_file);
        manifest_row=row;save(raw_file,'case_data','manifest_row','-v7.3');
        if isfile(partial_file),delete(partial_file);end
    end
    records(position)=summarize(case_data,row,raw_file);
    raw_files{position}=raw_file;
end
writetable(struct2table(records,'AsArray',true),final_csv);
writetable(struct2table(library_provenance,'AsArray',true), ...
    fullfile(round_dir,'outcome_library_provenance.csv'));
results=struct('configuration',cfg,'selection',selection,'records',records, ...
    'raw_files',{raw_files},'library_provenance',library_provenance);
save(final_file,'results','-v7.3');
brusselator_write_metadata(round_dir,cfg, ...
    brusselator_run_metadata(cfg,mfilename,toc(timer)));
disp(struct2table(records,'AsArray',true));
end

function record=empty_record()
record=struct('state_index',NaN,'sample_id','','lambda',NaN, ...
    'normalized_coordinate',NaN,'final_label','','label_status','', ...
    'raw_five_way_outcome','','tested_duration',NaN,'decision_time',NaN, ...
    'runtime_seconds',NaN,'late_median_A',NaN,'late_median_B',NaN, ...
    'late_median_C',NaN,'late_median_R5',NaN,'late_median_R6',NaN, ...
    'raw_file','');
end

function record=summarize(data,row,raw_file)
record=empty_record();record.state_index=row.state_index;
record.sample_id=char(row.sample_id);record.lambda=row.lambda;
record.normalized_coordinate=row.normalized_coordinate;
raw=data.classification.outcome;record.raw_five_way_outcome=raw;
if ismember(raw,{'A','B','C'})
    record.final_label=raw;record.label_status='resolved_primary';
elseif ismember(raw,{'R5','R6'})
    record.final_label='other';record.label_status='resolved_novelty';
else
    record.final_label='U';record.label_status='unresolved_or_ambiguous';
end
record.tested_duration=data.duration;record.decision_time=data.decision_time;
record.runtime_seconds=data.runtime_seconds;
values=data.classification.late_median_distances;
record.late_median_A=values(1);record.late_median_B=values(2);
record.late_median_C=values(3);record.late_median_R5=values(4);
record.late_median_R6=values(5);record.raw_file=raw_file;
end
