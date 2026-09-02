function results = run_ml_heldout_CBA_postbudget_refinement()
%RUN_ML_HELDOUT_CBA_POSTBUDGET_REFINEMENT Guard the final ML C/A pair.

cfg=brusselator_ml_heldout_CBA_postbudget_config();
outdir=fullfile(cfg.output.root,cfg.experiment_name);
final_file=fullfile(outdir,'postbudget_refinement.mat');
if isfolder(outdir)
    error('run_ml_heldout_CBA_postbudget_refinement:OutputExists', ...
        'Refusing to overwrite post-budget diagnostic: %s',outdir);
end
required={cfg.postbudget.source_selection,cfg.postbudget.source_labels, ...
    cfg.postbudget.parent_edge_file};
if ~all(cellfun(@isfile,required))
    error('run_ml_heldout_CBA_postbudget_refinement:MissingSource', ...
        'A final-pair source is missing.');
end
selection=readtable(cfg.postbudget.source_selection,'TextType','string');
labels=readtable(cfg.postbudget.source_labels,'TextType','string');
if height(selection)~=1||height(labels)~=2
    error('run_ml_heldout_CBA_postbudget_refinement:SourceShape', ...
        'The final selected pair must have exactly two endpoint labels.');
end
left_lambda=selection.first_lambda;right_lambda=selection.second_lambda;
left_label=char(labels.final_label(1));right_label=char(labels.final_label(2));
if left_lambda>right_lambda
    [left_lambda,right_lambda]=deal(right_lambda,left_lambda);
    [left_label,right_label]=deal(right_label,left_label);
end
if ~strcmp(left_label,'C')||~strcmp(right_label,'A')
    error('run_ml_heldout_CBA_postbudget_refinement:NotCA', ...
        'The final ML pair is not ordered C/A.');
end
parent_loaded=load(cfg.postbudget.parent_edge_file,'results');
parent=parent_loaded.results;
[library,library_provenance]=brusselator_build_ml_outcome_library(cfg);
mkdir(outdir);mkdir(fullfile(outdir,'raw_cases'));
brusselator_write_metadata(outdir,cfg, ...
    brusselator_run_metadata(cfg,mfilename,NaN));
timer=tic;records=empty_record([]);status='maximum_tests_reached_without_B';
B_lambda=NaN;B_file='';
for step=1:cfg.postbudget.maximum_midpoint_tests
    lambda=0.5*(left_lambda+right_lambda);
    case_id=sprintf('postbudget_step_%02d',step);
    raw_file=fullfile(outdir,'raw_cases',[case_id,'.mat']);
    data=brusselator_execute_edge_slice_case(case_id,lambda, ...
        cfg.postbudget.classification_duration,parent.endpoint_B, ...
        parent.endpoint_A,library,cfg);
    case_data=data;save(raw_file,'case_data','-v7.3');
    record=empty_record();record.step=step;record.lambda=lambda;
    record.outcome=data.classification.outcome;record.left_before=left_lambda;
    record.right_before=right_lambda;record.runtime_seconds=data.runtime_seconds;
    record.tested_duration=data.duration;record.decision_time=data.decision_time;
    record.raw_file=raw_file;
    switch data.classification.outcome
        case 'C'
            left_lambda=lambda;record.update='C_endpoint';
        case 'A'
            right_lambda=lambda;record.update='A_endpoint';
        case 'B'
            B_lambda=lambda;B_file=raw_file;record.update='B_intrusion_found';
            status='B_intrusion_found';
        otherwise
            record.update='non_ABC_obstruction';
            status=['obstructed_by_',data.classification.outcome];
    end
    record.left_after=left_lambda;record.right_after=right_lambda;
    record.width_after=right_lambda-left_lambda;
    records(end+1,1)=record; %#ok<AGROW>
    if strcmp(status,'B_intrusion_found')||startsWith(status,'obstructed_by_')
        break;
    end
end
brackets=empty_bracket([]);
if strcmp(status,'B_intrusion_found')
    brackets(1,1)=make_bracket('CB',left_lambda,'C',B_lambda,'B', ...
        source_for_lambda(left_lambda,records,cfg),B_file,cfg);
    brackets(2,1)=make_bracket('BA',B_lambda,'B',right_lambda,'A', ...
        B_file,source_for_lambda(right_lambda,records,cfg),cfg);
end
writetable(struct2table(records,'AsArray',true), ...
    fullfile(outdir,'postbudget_midpoint_cases.csv'));
writetable(struct2table(brackets,'AsArray',true), ...
    fullfile(outdir,'discovered_adjacent_brackets.csv'));
writetable(struct2table(library_provenance,'AsArray',true), ...
    fullfile(outdir,'outcome_library_provenance.csv'));
results=struct('configuration',cfg,'status',status,'records',records, ...
    'brackets',brackets,'library_provenance',library_provenance, ...
    'fixed_budget_primary_success',false, ...
    'postbudget_B_discovered',strcmp(status,'B_intrusion_found'));
save(final_file,'results','-v7.3');
brusselator_write_metadata(outdir,cfg, ...
    brusselator_run_metadata(cfg,mfilename,toc(timer)));
disp(struct2table(records,'AsArray',true));
disp(struct2table(brackets,'AsArray',true));
fprintf('Post-budget CBA refinement saved in: %s\n',outdir);
end

function record=empty_record(varargin)
record=struct('step',NaN,'lambda',NaN,'outcome','','left_before',NaN, ...
    'right_before',NaN,'left_after',NaN,'right_after',NaN,'width_after',NaN, ...
    'update','','tested_duration',NaN,'decision_time',NaN, ...
    'runtime_seconds',NaN,'raw_file','');
if nargin>0,record=record([]);end
end

function bracket=empty_bracket(varargin)
bracket=struct('bracket_id','','left_lambda',NaN,'left_outcome','', ...
    'right_lambda',NaN,'right_outcome','','width',NaN, ...
    'left_source_file','','right_source_file','','interpretation','');
if nargin>0,bracket=bracket([]);end
end

function bracket=make_bracket(id,left,left_label,right,right_label, ...
        left_file,right_file,cfg)
bracket=empty_bracket();bracket.bracket_id=id;bracket.left_lambda=left;
bracket.left_outcome=left_label;bracket.right_lambda=right;
bracket.right_outcome=right_label;bracket.width=right-left;
bracket.left_source_file=left_file;bracket.right_source_file=right_file;
bracket.interpretation=cfg.postbudget.claim_scope;
end

function file=source_for_lambda(lambda,records,cfg)
index=find(abs([records.lambda]-lambda)<1e-12,1,'last');
if isempty(index)
    labels=readtable(cfg.postbudget.source_labels,'TextType','string');
    [distance,row]=min(abs(labels.lambda-lambda));
    if distance>1e-12,error('Endpoint source cannot be resolved.');end
    file=char(labels.raw_file(row));
else
    file=records(index).raw_file;
end
end
