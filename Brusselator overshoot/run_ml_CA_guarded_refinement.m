function results = run_ml_CA_guarded_refinement()
%RUN_ML_CA_GUARDED_REFINEMENT Refine the evolved-chord C--A candidate.

cfg=brusselator_ml_CA_closure_config();
outdir=fullfile(cfg.output.root,cfg.experiment_name);
checkpoint=fullfile(outdir,'progress_CA_refinement.mat');
final_file=fullfile(outdir,'CA_guarded_refinement.mat');
if isfile(final_file)
    error('run_ml_CA_guarded_refinement:OutputExists', ...
        'Refusing to overwrite completed refinement: %s',final_file);
end
if ~isfile(cfg.CA_closure.source_scan_file)
    error('run_ml_CA_guarded_refinement:MissingScan','Source scan is missing.');
end
loaded=load(cfg.CA_closure.source_scan_file,'results');scan=loaded.results;
[library,library_provenance]=brusselator_build_ml_outcome_library(cfg);
[left_alpha,right_alpha,left_file,right_file]=initialize_bracket(scan,cfg);
if isfile(checkpoint)
    saved=load(checkpoint,'tracker','records','raw_files');
    tracker=saved.tracker;records=saved.records;raw_files=saved.raw_files;
    validate_tracker(tracker,cfg);
    fprintf('Resuming C--A refinement after %d cases.\n',numel(records));
else
    if isfolder(outdir)
        error('run_ml_CA_guarded_refinement:UnsafeOutput', ...
            'Existing output has no recognized checkpoint: %s',outdir);
    end
    mkdir(outdir);mkdir(fullfile(outdir,'raw_cases'));
    tracker=struct('status','running','next_step',1, ...
        'left_alpha',left_alpha,'right_alpha',right_alpha, ...
        'left_outcome','C','right_outcome','A', ...
        'left_source_file',left_file,'right_source_file',right_file);
    records=empty_record([]);raw_files={};
    brusselator_write_metadata(outdir,cfg, ...
        brusselator_run_metadata(cfg,mfilename,NaN));
    write_checkpoint();
end
timer=tic;
while strcmp(tracker.status,'running')&& ...
        tracker.next_step<=cfg.CA_closure.maximum_bisection_steps
    step=tracker.next_step;before_left=tracker.left_alpha;
    before_right=tracker.right_alpha;alpha=0.5*(before_left+before_right);
    initial_state=interpolate_states(scan.endpoint_left,scan.endpoint_right,alpha);
    case_id=sprintf('CA_closure_step_%02d',step);
    raw_file=fullfile(outdir,'raw_cases',[case_id,'.mat']);
    partial_file=fullfile(outdir,'raw_cases',[case_id,'_partial.mat']);
    fprintf('C--A closure %d/%d: alpha=%.10f in [%.10f, %.10f].\n', ...
        step,cfg.CA_closure.maximum_bisection_steps,alpha,before_left,before_right);
    if isfile(raw_file)
        saved=load(raw_file,'closure_case');data=saved.closure_case;
    else
        data=brusselator_evolve_and_classify_state( ...
            case_id,initial_state,library,cfg,partial_file);
        closure_case=data;save(raw_file,'closure_case','-v7.3');
        if isfile(partial_file),delete(partial_file);end
    end
    outcome=data.classification.outcome;update='obstruction';
    if strcmp(outcome,'C')
        tracker.left_alpha=alpha;tracker.left_source_file=raw_file;
        update='left_C_endpoint';
    elseif strcmp(outcome,'A')
        tracker.right_alpha=alpha;tracker.right_source_file=raw_file;
        update='right_A_endpoint';
    else
        tracker.status=['obstructed_by_',outcome];
    end
    record=summarize_case(data,step,alpha,before_left,before_right, ...
        tracker.left_alpha,tracker.right_alpha,update,raw_file);
    records(end+1,1)=record;raw_files{end+1,1}=raw_file; %#ok<AGROW>
    tracker.next_step=step+1;
    if strcmp(tracker.status,'running')&& ...
            tracker.next_step>cfg.CA_closure.maximum_bisection_steps
        tracker.status='clean_CA_label_bracket_refined';
    end
    write_checkpoint();
end
summary=make_summary(tracker,records,cfg);
bracket=make_bracket(tracker,cfg);
writetable(struct2table(records,'AsArray',true), ...
    fullfile(outdir,'CA_refinement_cases.csv'));
writetable(struct2table(bracket,'AsArray',true), ...
    fullfile(outdir,'CA_refined_bracket.csv'));
writetable(struct2table(summary,'AsArray',true), ...
    fullfile(outdir,'CA_refinement_summary.csv'));
writetable(struct2table(library_provenance,'AsArray',true), ...
    fullfile(outdir,'outcome_library_provenance.csv'));
make_figure(outdir,scan,records,tracker);
results=struct('configuration',cfg,'source_scan_file', ...
    cfg.CA_closure.source_scan_file,'tracker',tracker,'records',records, ...
    'raw_files',{raw_files},'bracket',bracket,'summary',summary, ...
    'library_provenance',library_provenance);
save(final_file,'results','-v7.3');
brusselator_write_metadata(outdir,cfg, ...
    brusselator_run_metadata(cfg,mfilename,toc(timer)));
disp(struct2table(records,'AsArray',true));
disp(struct2table(bracket,'AsArray',true));
disp(struct2table(summary,'AsArray',true));
fprintf('C--A guarded refinement saved in: %s\n',outdir);

    function write_checkpoint()
        save(checkpoint,'cfg','tracker','records','raw_files','-v7.3');
        if ~isempty(records)
            writetable(struct2table(records,'AsArray',true), ...
                fullfile(outdir,'progress_summary.csv'));
        end
    end
end

function [left,right,left_file,right_file]=initialize_bracket(scan,cfg)
match=strcmp({scan.transition_records.left_outcome},'C')& ...
    strcmp({scan.transition_records.right_outcome},'A');
if sum(match)~=1,error('Source scan does not have exactly one C--A interval.');end
transition=scan.transition_records(match);left=transition.left_alpha;
right=transition.right_alpha;
if abs(left-cfg.CA_closure.expected_left_alpha)>1e-12|| ...
        abs(right-cfg.CA_closure.expected_right_alpha)>1e-12
    error('Source C--A interval differs from the frozen closure design.');
end
left_index=find(abs([scan.records.alpha]-left)<1e-12,1);
right_index=find(abs([scan.records.alpha]-right)<1e-12,1);
if isempty(left_index)||isempty(right_index)|| ...
        ~strcmp(scan.records(left_index).outcome,'C')|| ...
        ~strcmp(scan.records(right_index).outcome,'A')
    error('Source endpoint records do not provide the ordered C/A labels.');
end
left_file=scan.raw_files{left_index};right_file=scan.raw_files{right_index};
end

function validate_tracker(tracker,cfg)
if ~strcmp(tracker.left_outcome,cfg.CA_closure.left_outcome)|| ...
        ~strcmp(tracker.right_outcome,cfg.CA_closure.right_outcome)|| ...
        tracker.left_alpha>=tracker.right_alpha
    error('Saved C--A tracker is incompatible with the configuration.');
end
end

function state=interpolate_states(left,right,alpha)
state=left;state(:,2:3)=(1-alpha)*left(:,2:3)+alpha*right(:,2:3);
end

function record=empty_record(varargin)
record=struct('step',NaN,'alpha',NaN,'outcome','', ...
    'tested_duration',NaN,'decision_time',NaN,'runtime_seconds',NaN, ...
    'left_before',NaN,'right_before',NaN,'left_after',NaN,'right_after',NaN, ...
    'width_after',NaN,'update','','late_distance_A',NaN, ...
    'late_distance_B',NaN,'late_distance_C',NaN,'late_distance_R5',NaN, ...
    'late_distance_R6',NaN,'raw_file','');
if nargin>0,record=record([]);end
end

function record=summarize_case(data,step,alpha,before_left,before_right, ...
        after_left,after_right,update,raw_file)
record=empty_record();record.step=step;record.alpha=alpha;
record.outcome=data.classification.outcome;record.tested_duration=data.duration;
record.decision_time=data.decision_time;record.runtime_seconds=data.runtime_seconds;
record.left_before=before_left;record.right_before=before_right;
record.left_after=after_left;record.right_after=after_right;
record.width_after=after_right-after_left;record.update=update;
labels={'A','B','C','R5','R6'};
for k=1:numel(labels)
    record.(['late_distance_',labels{k}])= ...
        data.classification.late_median_distances(k);
end
record.raw_file=raw_file;
end

function bracket=make_bracket(tracker,cfg)
bracket=struct('bracket_id','CA','left_alpha',tracker.left_alpha, ...
    'left_outcome','C','right_alpha',tracker.right_alpha, ...
    'right_outcome','A','width',tracker.right_alpha-tracker.left_alpha, ...
    'left_source_file',tracker.left_source_file, ...
    'right_source_file',tracker.right_source_file, ...
    'status',tracker.status,'interpretation',cfg.CA_closure.claim_scope);
end

function summary=make_summary(tracker,records,cfg)
labels={records.outcome};summary=struct('status',tracker.status, ...
    'steps_completed',numel(records),'final_width', ...
    tracker.right_alpha-tracker.left_alpha,'target_width', ...
    cfg.CA_closure.target_width,'clean_CA', ...
    strcmp(tracker.status,'clean_CA_label_bracket_refined'), ...
    'number_C',sum(strcmp(labels,'C')),'number_A',sum(strcmp(labels,'A')), ...
    'number_other',sum(~strcmp(labels,'C')&~strcmp(labels,'A')), ...
    'interpretation',cfg.CA_closure.claim_scope);
end

function make_figure(outdir,scan,records,tracker)
fig=figure('Visible','off','Color','w','Position',[100 100 1000 500]);
ax=axes(fig);hold(ax,'on');order={'B','C','A'};
for k=1:numel(scan.records)
    y=find(strcmp(order,scan.records(k).outcome),1);
    scatter(ax,scan.records(k).alpha,y,45,[0.75 0.75 0.75],'filled');
end
for k=1:numel(records)
    y=find(strcmp(order,records(k).outcome),1);if isempty(y),y=NaN;end
    scatter(ax,records(k).alpha,y,85,outcome_color(records(k).outcome), ...
        'filled','MarkerEdgeColor','k');
end
xline(ax,tracker.left_alpha,'--','C endpoint');
xline(ax,tracker.right_alpha,'--','A endpoint');
yticks(ax,1:3);yticklabels(ax,order);ylim(ax,[0.7 3.3]);xlim(ax,[0.80 0.92]);
grid(ax,'on');box(ax,'on');xlabel(ax,'alpha on evolved B--A chord');
ylabel(ax,'five-way outcome');title(ax,'Guarded closure of the C--A interval');
exportgraphics(fig,fullfile(outdir,'CA_guarded_refinement.png'),'Resolution',300);
close(fig);
end

function color=outcome_color(label)
switch label
    case 'A',color=[0.10 0.35 0.80];
    case 'B',color=[0.90 0.35 0.10];
    case 'C',color=[0.85 0.65 0.05];
    case 'R5',color=[0.50 0.20 0.75];
    case 'R6',color=[0.10 0.60 0.30];
    otherwise,color=[0.45 0.45 0.45];
end
end
