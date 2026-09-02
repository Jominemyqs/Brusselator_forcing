function results = run_edge_slice_pairwise_refinement()
%RUN_EDGE_SLICE_PAIRWISE_REFINEMENT Refine B--C and C--A label transitions.
%   Every midpoint is classified against A/B/C/R5/R6. A bracket is updated
%   only for one of its endpoint labels; any other or unresolved result stops
%   that refinement and is saved rather than coerced into a binary outcome.

cfg=brusselator_edge_slice_refinement_config();study_id=cfg.experiment_name;
outdir=fullfile(cfg.output.root,study_id);
checkpoint=fullfile(outdir,'progress_pairwise_refinement.mat');
final_file=fullfile(outdir,'edge_slice_pairwise_refinement.mat');
if ~isfile(cfg.refinement.scan_file)
    error('run_edge_slice_pairwise_refinement:MissingScan', ...
        'Expected completed slice scan at %s.',cfg.refinement.scan_file);
end
if isfile(final_file)
    error('run_edge_slice_pairwise_refinement:OutputExists', ...
        'Refusing to overwrite completed refinement: %s',final_file);
end
loaded=load(cfg.refinement.scan_file,'results');scan=loaded.results;
parent_loaded=load(scan.parent_edge_file,'results');parent=parent_loaded.results;
validate_inputs(scan,parent,cfg);
if isfile(checkpoint)
    saved=load(checkpoint,'pair_states','records','raw_files');
    pair_states=saved.pair_states;records=saved.records;raw_files=saved.raw_files;
    fprintf('Resuming %s with %d saved midpoint cases.\n',study_id,numel(records));
else
    if isfolder(outdir)
        error('run_edge_slice_pairwise_refinement:IncompleteOutput', ...
            'Existing output lacks a valid checkpoint: %s',outdir);
    end
    mkdir(outdir);mkdir(fullfile(outdir,'raw_cases'));
    pair_states=initialize_pairs(scan,cfg);template=empty_record();
    records=template([]);raw_files={};
    brusselator_write_metadata(outdir,cfg, ...
        brusselator_run_metadata(cfg,mfilename,NaN));
    write_checkpoint(checkpoint,outdir,cfg,pair_states,records,raw_files);
end
timer=tic;
for p=1:numel(pair_states)
    while strcmp(pair_states(p).status,'running') && ...
            pair_states(p).next_step<=cfg.refinement.maximum_bisection_steps
        step=pair_states(p).next_step;
        before_left=pair_states(p).left_lambda;
        before_right=pair_states(p).right_lambda;
        lambda=0.5*(before_left+before_right);
        case_id=sprintf('%s_step_%02d',lower(pair_states(p).pair_id),step);
        fprintf('%s refinement %d/%d: lambda=%.10f in [%.10f, %.10f].\n', ...
            pair_states(p).pair_id,step,cfg.refinement.maximum_bisection_steps, ...
            lambda,before_left,before_right);
        raw_file=fullfile(outdir,'raw_cases',[case_id,'.mat']);
        if isfile(raw_file)
            saved=load(raw_file,'refinement_case');data=saved.refinement_case;
            fprintf('  Recovered completed midpoint %s.\n',case_id);
        else
            data=brusselator_execute_edge_slice_case(case_id,lambda, ...
                cfg.refinement.classification_duration,parent.endpoint_B, ...
                parent.endpoint_A,parent.library,cfg);
            refinement_case=data;save(raw_file,'refinement_case','-v7.3');
        end
        outcome=data.classification.outcome;
        update='obstruction';
        if strcmp(outcome,pair_states(p).left_outcome)
            pair_states(p).left_lambda=lambda;update='left_endpoint';
        elseif strcmp(outcome,pair_states(p).right_outcome)
            pair_states(p).right_lambda=lambda;update='right_endpoint';
        else
            pair_states(p).status=['obstructed_by_',outcome];
        end
        record=summarize_case(data,pair_states(p),step,before_left, ...
            before_right,update,raw_file);
        records(end+1,1)=record;raw_files{end+1,1}=raw_file; %#ok<AGROW>
        pair_states(p).next_step=step+1;
        pair_states(p).last_case_id=case_id;
        if strcmp(pair_states(p).status,'running') && ...
                pair_states(p).next_step>cfg.refinement.maximum_bisection_steps
            pair_states(p).status='finite_time_label_bracket_refined';
        end
        write_checkpoint(checkpoint,outdir,cfg,pair_states,records,raw_files);
    end
end

summaries=make_summaries(pair_states,cfg);
writetable(struct2table(records,'AsArray',true), ...
    fullfile(outdir,'pairwise_refinement_cases.csv'));
writetable(struct2table(summaries,'AsArray',true), ...
    fullfile(outdir,'pairwise_refined_brackets.csv'));
make_figure(outdir,scan,records,pair_states);
results=struct('configuration',cfg,'scan_file',cfg.refinement.scan_file, ...
    'pair_states',pair_states,'records',records,'raw_files',{raw_files}, ...
    'summaries',summaries);
save(final_file,'results','-v7.3');
brusselator_write_metadata(outdir,cfg, ...
    brusselator_run_metadata(cfg,mfilename,toc(timer)));
disp(struct2table(records,'AsArray',true));
disp(struct2table(summaries,'AsArray',true));
fprintf('Pairwise slice refinements saved in: %s\n',outdir);
end

function validate_inputs(scan,parent,cfg)
if scan.configuration.grid.N~=cfg.grid.N || parent.configuration.grid.N~=cfg.grid.N
    error('run_edge_slice_pairwise_refinement:GridMismatch', ...
        'Scan, parent endpoints, and refinement grid do not agree.');
end
if numel(scan.transition_records)~=numel(cfg.refinement.pairs)
    error('run_edge_slice_pairwise_refinement:TransitionCount', ...
        'The completed scan does not contain exactly the configured transitions.');
end
for k=1:numel(cfg.refinement.pairs)
    requested=cfg.refinement.pairs(k);
    match=strcmp({scan.transition_records.left_outcome},requested.left_outcome) & ...
        strcmp({scan.transition_records.right_outcome},requested.right_outcome);
    if sum(match)~=1
        error('run_edge_slice_pairwise_refinement:TransitionMismatch', ...
            'Expected one %s--%s interval in the scan.', ...
            requested.left_outcome,requested.right_outcome);
    end
end
if ~isequal({parent.library.label},{'A','B','C','R5','R6'})
    error('run_edge_slice_pairwise_refinement:LibraryMismatch', ...
        'The physical five-way library is not A/B/C/R5/R6.');
end
end

function states=initialize_pairs(scan,cfg)
template=struct('pair_id','','left_outcome','','right_outcome','', ...
    'initial_left_lambda',NaN,'initial_right_lambda',NaN, ...
    'left_lambda',NaN,'right_lambda',NaN,'left_source_file','', ...
    'right_source_file','','next_step',1,'status','running', ...
    'last_case_id','');
states=repmat(template,numel(cfg.refinement.pairs),1);
for k=1:numel(states)
    requested=cfg.refinement.pairs(k);
    index=find(strcmp({scan.transition_records.left_outcome}, ...
        requested.left_outcome) & strcmp({scan.transition_records.right_outcome}, ...
        requested.right_outcome),1);
    transition=scan.transition_records(index);
    left_index=find(abs([scan.records.lambda]-transition.left_lambda)<1e-12,1);
    right_index=find(abs([scan.records.lambda]-transition.right_lambda)<1e-12,1);
    states(k).pair_id=requested.pair_id;
    states(k).left_outcome=requested.left_outcome;
    states(k).right_outcome=requested.right_outcome;
    states(k).initial_left_lambda=transition.left_lambda;
    states(k).initial_right_lambda=transition.right_lambda;
    states(k).left_lambda=transition.left_lambda;
    states(k).right_lambda=transition.right_lambda;
    states(k).left_source_file=scan.raw_files{left_index};
    states(k).right_source_file=scan.raw_files{right_index};
end
end

function record=empty_record()
record=struct('pair_id','','step',NaN,'case_id','','lambda',NaN, ...
    'outcome','','duration',NaN,'was_extended',false,'decision_time',NaN, ...
    'runtime_seconds',NaN,'left_before',NaN,'right_before',NaN, ...
    'left_after',NaN,'right_after',NaN,'width_after',NaN,'update','', ...
    'late_median_distance_A',NaN,'late_median_distance_B',NaN, ...
    'late_median_distance_C',NaN,'late_median_distance_R5',NaN, ...
    'late_median_distance_R6',NaN,'final_distance_A',NaN, ...
    'final_distance_B',NaN,'final_distance_C',NaN,'final_distance_R5',NaN, ...
    'final_distance_R6',NaN,'recurrence_available',false, ...
    'recurrence_period',NaN,'recurrence_distance',NaN,'raw_file','');
end

function record=summarize_case(data,state,step,left_before,right_before,update,raw)
record=empty_record();record.pair_id=state.pair_id;record.step=step;
record.case_id=data.case_id;record.lambda=data.lambda;
record.outcome=data.classification.outcome;record.duration=data.duration;
record.was_extended=data.was_extended;record.decision_time=data.decision_time;
record.runtime_seconds=data.runtime_seconds;record.left_before=left_before;
record.right_before=right_before;record.left_after=state.left_lambda;
record.right_after=state.right_lambda;
record.width_after=state.right_lambda-state.left_lambda;record.update=update;
labels={'A','B','C','R5','R6'};
for k=1:numel(labels)
    record.(['late_median_distance_',labels{k}])= ...
        data.classification.late_median_distances(k);
    record.(['final_distance_',labels{k}])= ...
        data.classification.final_distances(k);
end
record.recurrence_available=data.recurrence.available;
record.recurrence_period=data.recurrence.period;
record.recurrence_distance=data.recurrence.distance;record.raw_file=raw;
end

function summaries=make_summaries(states,cfg)
template=struct('pair_id','','left_outcome','','right_outcome','', ...
    'initial_left_lambda',NaN,'initial_right_lambda',NaN, ...
    'final_left_lambda',NaN,'final_right_lambda',NaN,'final_width',NaN, ...
    'steps_completed',NaN,'target_width',cfg.refinement.target_width, ...
    'status','','interpretation',cfg.refinement.interpretation);
summaries=repmat(template,numel(states),1);
for k=1:numel(states)
    summaries(k).pair_id=states(k).pair_id;
    summaries(k).left_outcome=states(k).left_outcome;
    summaries(k).right_outcome=states(k).right_outcome;
    summaries(k).initial_left_lambda=states(k).initial_left_lambda;
    summaries(k).initial_right_lambda=states(k).initial_right_lambda;
    summaries(k).final_left_lambda=states(k).left_lambda;
    summaries(k).final_right_lambda=states(k).right_lambda;
    summaries(k).final_width=states(k).right_lambda-states(k).left_lambda;
    summaries(k).steps_completed=states(k).next_step-1;
    summaries(k).status=states(k).status;
end
end

function write_checkpoint(file,outdir,cfg,pair_states,records,raw_files)
save(file,'cfg','pair_states','records','raw_files','-v7.3');
if ~isempty(records)
    writetable(struct2table(records,'AsArray',true), ...
        fullfile(outdir,'progress_summary.csv'));
end
end

function make_figure(outdir,scan,records,states)
fig=figure('Color','w','Position',[100 100 1050 500]);ax=axes(fig);hold(ax,'on');
order={'B','C','A'};scan_y=nan(size(scan.records));
for k=1:numel(scan.records)
    scan_y(k)=find(strcmp(order,scan.records(k).outcome),1);
end
plot(ax,[scan.records.lambda],scan_y,'o-','Color',[0.65 0.65 0.65], ...
    'MarkerFaceColor',[0.85 0.85 0.85],'DisplayName','dense scan');
for k=1:numel(records)
    y=find(strcmp(order,records(k).outcome),1);
    if isempty(y),y=NaN;end
    scatter(ax,records(k).lambda,y,65,outcome_color(records(k).outcome), ...
        'filled','MarkerEdgeColor','k','HandleVisibility','off');
end
ylimits=[0.6 3.4];
for k=1:numel(states)
    xregion=[states(k).left_lambda states(k).right_lambda];
    patch(ax,[xregion fliplr(xregion)],[ylimits(1) ylimits(1) ylimits(2) ylimits(2)], ...
        [0.2 0.2 0.2],'FaceAlpha',0.12,'EdgeColor','none', ...
        'HandleVisibility','off');
end
yticks(ax,1:3);yticklabels(ax,order);ylim(ax,ylimits);xlim(ax,[0 0.5]);
xlabel(ax,'lambda');ylabel(ax,'five-way outcome');
title(ax,'Multiclass-guarded refinement of slice transitions');
style_axes(ax);legend(ax,'Location','best');
exportgraphics(fig,fullfile(outdir,'pairwise_refined_brackets.png'),'Resolution',300);
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

function style_axes(ax)
set(ax,'Color','w','XColor','k','YColor','k','GridColor',[0.75 0.75 0.75]);
grid(ax,'on');box(ax,'on');ax.Title.Color='k';ax.XLabel.Color='k';ax.YLabel.Color='k';
end
