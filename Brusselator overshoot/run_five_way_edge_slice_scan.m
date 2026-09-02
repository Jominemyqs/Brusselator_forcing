function results = run_five_way_edge_slice_scan()
%RUN_FIVE_WAY_EDGE_SLICE_SCAN Map outcomes before pairwise edge bisection.
%   Scans the reflection-aligned B-to-A landing-state segment with the
%   physical A/B/C/R5/R6 classifier. Existing lambda=0,0.25,0.5 edge cases
%   are reused exactly; all other trajectories are checkpointed separately.

cfg=brusselator_edge_slice_scan_config();study_id=cfg.experiment_name;
outdir=fullfile(cfg.output.root,study_id);
checkpoint=fullfile(outdir,'progress_edge_slice_scan.mat');
final_file=fullfile(outdir,'edge_slice_scan_results.mat');
if ~isfile(cfg.scan.parent_edge_file)
    error('run_five_way_edge_slice_scan:MissingParent', ...
        'Expected parent edge result at %s.',cfg.scan.parent_edge_file);
end
if isfile(final_file)
    error('run_five_way_edge_slice_scan:OutputExists', ...
        'Refusing to overwrite completed scan: %s',final_file);
end
parent=load(cfg.scan.parent_edge_file,'results');parent=parent.results;
validate_parent(parent,cfg);
lambdas=cfg.scan.lambda_values(:);
if isfile(checkpoint)
    saved=load(checkpoint,'records','raw_files');
    records=saved.records;raw_files=saved.raw_files;
    if numel(records)~=numel(lambdas)
        error('run_five_way_edge_slice_scan:CheckpointMismatch', ...
            'Saved checkpoint does not match the configured lambda grid.');
    end
    fprintf('Resuming %s with %d/%d completed cases.\n',study_id, ...
        sum([records.completed]),numel(records));
else
    if isfolder(outdir)
        raw_dir=fullfile(outdir,'raw_cases');
        if isfolder(raw_dir)&&~isempty(dir(fullfile(raw_dir,'*.mat')))
            error('run_five_way_edge_slice_scan:IncompleteOutput', ...
                ['Existing output has raw cases but no valid checkpoint; ', ...
                'manual inspection is required: %s'],outdir);
        end
        fprintf('Restarting scan initialization in metadata-only output: %s\n',outdir);
    else
        mkdir(outdir);
    end
    if ~isfolder(fullfile(outdir,'raw_cases')),mkdir(fullfile(outdir,'raw_cases'));end
    records=repmat(empty_record(),numel(lambdas),1);raw_files=cell(numel(lambdas),1);
    for k=1:numel(lambdas),records(k).scan_index=k;records(k).lambda=lambdas(k);end
    brusselator_write_metadata(outdir,cfg, ...
        brusselator_run_metadata(cfg,mfilename,NaN));
end
timer=tic;
for k=1:numel(lambdas)
    if records(k).completed && isfile(raw_files{k})
        fprintf('Skipping completed scan case %d/%d at lambda=%.3f.\n', ...
            k,numel(lambdas),lambdas(k));
        continue;
    end
    [reuse,parent_index]=find_parent_case(parent,lambdas(k),cfg);
    if reuse
        source_file=parent.raw_files{parent_index};
        loaded=load(source_file,'case_result');case_data=loaded.case_result;
        source_type='reused_parent_edge_case';
        source_case_id=parent.records(parent_index).case_id;
        fprintf('Reusing parent %s at lambda=%.3f.\n',source_case_id,lambdas(k));
    else
        case_id=sprintf('lambda_%05.3f',lambdas(k));
        source_file=fullfile(outdir,'raw_cases',[strrep(case_id,'.','p'),'.mat']);
        if isfile(source_file)
            loaded=load(source_file,'scan_case');case_data=loaded.scan_case;
            fprintf('Recovered raw scan case at lambda=%.3f.\n',lambdas(k));
        else
            case_data=brusselator_execute_edge_slice_case(case_id,lambdas(k), ...
                cfg.scan.classification_duration,parent.endpoint_B,parent.endpoint_A, ...
                parent.library,cfg);
            scan_case=case_data;save(source_file,'scan_case','-v7.3');
        end
        source_type='new_scan_case';source_case_id=case_data.case_id;
    end
    records(k)=summarize_case(records(k),case_data,source_type, ...
        source_case_id,source_file);
    raw_files{k}=source_file;
    write_progress(outdir,checkpoint,cfg,records,raw_files);
end

transition_records=find_transitions(records);
region_records=find_regions(records);
summary=summarize_scan(records,transition_records,region_records,cfg);
writetable(struct2table(records,'AsArray',true), ...
    fullfile(outdir,'five_way_slice_scan.csv'));
writetable(struct2table(transition_records,'AsArray',true), ...
    fullfile(outdir,'candidate_transition_intervals.csv'));
writetable(struct2table(region_records,'AsArray',true), ...
    fullfile(outdir,'contiguous_outcome_regions.csv'));
writetable(struct2table(summary,'AsArray',true), ...
    fullfile(outdir,'five_way_slice_scan_summary.csv'));
make_figures(outdir,records,parent.library,cfg);
results=struct('configuration',cfg,'parent_edge_file',cfg.scan.parent_edge_file, ...
    'records',records,'raw_files',{raw_files}, ...
    'transition_records',transition_records,'region_records',region_records, ...
    'summary',summary);
save(final_file,'results','-v7.3');
brusselator_write_metadata(outdir,cfg, ...
    brusselator_run_metadata(cfg,mfilename,toc(timer)));
disp(struct2table(records,'AsArray',true));
disp(struct2table(transition_records,'AsArray',true));
fprintf('Five-way edge slice scan saved in: %s\n',outdir);
end

function validate_parent(parent,cfg)
if parent.configuration.grid.N~=cfg.grid.N || ...
        abs(parent.configuration.grid.Lx-cfg.grid.Lx)>1e-12 || ...
        ~isequal({parent.library.label},{'A','B','C','R5','R6'})
    error('run_five_way_edge_slice_scan:ParentMismatch', ...
        'Parent grid or five-way library does not match the scan.');
end
end

function record=empty_record()
record=struct('scan_index',NaN,'lambda',NaN,'completed',false, ...
    'case_id','','source_type','','source_case_id','','source_file','', ...
    'duration',NaN,'was_extended',false,'outcome','unclassified', ...
    'decision_time',NaN,'runtime_seconds',NaN, ...
    'initial_distance_A',NaN,'initial_distance_B',NaN,'initial_distance_C',NaN, ...
    'initial_distance_R5',NaN,'initial_distance_R6',NaN, ...
    'late_median_distance_A',NaN,'late_median_distance_B',NaN, ...
    'late_median_distance_C',NaN,'late_median_distance_R5',NaN, ...
    'late_median_distance_R6',NaN,'late_maximum_distance_A',NaN, ...
    'late_maximum_distance_B',NaN,'late_maximum_distance_C',NaN, ...
    'late_maximum_distance_R5',NaN,'late_maximum_distance_R6',NaN, ...
    'final_distance_A',NaN,'final_distance_B',NaN,'final_distance_C',NaN, ...
    'final_distance_R5',NaN,'final_distance_R6',NaN, ...
    'recurrence_available',false,'recurrence_period',NaN, ...
    'recurrence_distance',NaN);
end

function [reuse,index]=find_parent_case(parent,lambda,cfg)
reuse=false;index=[];
if ~any(abs(cfg.scan.reuse_parent_lambdas-lambda)<1e-12),return;end
index=find(abs([parent.records.lambda]-lambda)<1e-12,1);
reuse=~isempty(index);
end

function record=summarize_case(record,data,source_type,source_case_id,source_file)
record.completed=true;record.case_id=data.case_id;record.source_type=source_type;
record.source_case_id=source_case_id;record.source_file=source_file;
record.duration=case_value(data,'duration');
record.was_extended=case_value(data,'was_extended');
record.outcome=data.classification.outcome;
record.decision_time=case_value(data,'decision_time');
record.runtime_seconds=case_value(data,'runtime_seconds');
names={'A','B','C','R5','R6'};
for j=1:numel(names)
    record.(['initial_distance_',names{j}])=data.history.distances(1,j);
    record.(['late_median_distance_',names{j}])= ...
        data.classification.late_median_distances(j);
    record.(['late_maximum_distance_',names{j}])= ...
        data.classification.late_maximum_distances(j);
    record.(['final_distance_',names{j}])=data.classification.final_distances(j);
end
if isfield(data,'recurrence')&&isfield(data.recurrence,'available')
    record.recurrence_available=data.recurrence.available;
    record.recurrence_period=data.recurrence.period;
    record.recurrence_distance=data.recurrence.distance;
end
end

function value=case_value(data,name)
% Legacy pilot files stored run bookkeeping in case_result.record.
if isfield(data,name)
    value=data.(name);
elseif isfield(data,'record')&&isfield(data.record,name)
    value=data.record.(name);
else
    error('run_five_way_edge_slice_scan:MissingCaseField', ...
        'Case %s is missing required field %s.',data.case_id,name);
end
end

function records=find_transitions(scan)
template=struct('left_lambda',NaN,'right_lambda',NaN,'width',NaN, ...
    'left_outcome','','right_outcome','','pair','','status','');
records=template([]);
known={'A','B','C','R5','R6'};
for k=1:numel(scan)-1
    if strcmp(scan(k).outcome,scan(k+1).outcome),continue;end
    item=template;item.left_lambda=scan(k).lambda;item.right_lambda=scan(k+1).lambda;
    item.width=scan(k+1).lambda-scan(k).lambda;
    item.left_outcome=scan(k).outcome;item.right_outcome=scan(k+1).outcome;
    item.pair=[scan(k).outcome,'--',scan(k+1).outcome];
    if ismember(scan(k).outcome,known)&&ismember(scan(k+1).outcome,known)
        item.status='candidate_adjacent_at_scan_resolution';
    else
        item.status='contains_unresolved_or_ambiguous_endpoint';
    end
    records(end+1,1)=item; %#ok<AGROW>
end
end

function regions=find_regions(scan)
template=struct('outcome','','lambda_minimum',NaN,'lambda_maximum',NaN, ...
    'sample_count',NaN);
regions=template([]);start=1;
for k=2:numel(scan)+1
    if k<=numel(scan)&&strcmp(scan(k).outcome,scan(start).outcome),continue;end
    item=template;item.outcome=scan(start).outcome;
    item.lambda_minimum=scan(start).lambda;item.lambda_maximum=scan(k-1).lambda;
    item.sample_count=k-start;regions(end+1,1)=item; %#ok<AGROW>
    start=k;
end
end

function summary=summarize_scan(records,transitions,regions,cfg)
labels={records.outcome};
summary=struct('number_samples',numel(records),'lambda_step',cfg.scan.lambda_step, ...
    'number_A',sum(strcmp(labels,'A')),'number_B',sum(strcmp(labels,'B')), ...
    'number_C',sum(strcmp(labels,'C')),'number_R5',sum(strcmp(labels,'R5')), ...
    'number_R6',sum(strcmp(labels,'R6')), ...
    'number_unresolved',sum(strcmp(labels,'unresolved')), ...
    'number_contiguous_regions',numel(regions), ...
    'number_candidate_transition_intervals',numel(transitions), ...
    'interpretation',cfg.scan.transition_rule);
end

function write_progress(outdir,file,cfg,records,raw_files)
completed=records([records.completed]);
writetable(struct2table(completed,'AsArray',true), ...
    fullfile(outdir,'progress_summary.csv'));
save(file,'cfg','records','raw_files','-v7.3');
end

function make_figures(outdir,records,library,cfg)
labels={records.outcome};order={'B','R5','C','R6','A','unresolved', ...
    'ambiguous_multiple_outcomes'};y=zeros(size(labels));
for k=1:numel(labels)
    index=find(strcmp(order,labels{k}),1);if isempty(index),index=numel(order);end
    y(k)=index;
end
fig=figure('Color','w','Position',[100 100 1000 500]);ax=axes(fig);hold(ax,'on');
for k=1:numel(records)
    scatter(ax,records(k).lambda,y(k),75,outcome_color(records(k).outcome), ...
        'filled','MarkerEdgeColor','k');
end
plot(ax,[records.lambda],y,'-k','LineWidth',0.8,'HandleVisibility','off');
yticks(ax,1:numel(order));yticklabels(ax,order);set(ax,'TickLabelInterpreter','none');
xlabel(ax,'lambda: 0 = B landing, 0.5 = interior A-basin sample');
ylabel(ax,'five-way outcome');title(ax,'Multiclass frozen-system landing-state slice');
style_axes(ax);exportgraphics(fig,fullfile(outdir,'five_way_outcome_by_lambda.png'), ...
    'Resolution',300);close(fig);

distances=[[records.late_median_distance_A]',[records.late_median_distance_B]', ...
    [records.late_median_distance_C]',[records.late_median_distance_R5]', ...
    [records.late_median_distance_R6]'];
fig=figure('Color','w','Position',[100 100 950 600]);ax=axes(fig);
imagesc(ax,[records.lambda],1:numel(library),distances');colorbar(ax);
yticks(ax,1:numel(library));yticklabels(ax,{library.label});
xlabel(ax,'lambda');ylabel(ax,'validated outcome');
title(ax,'Late median physical distances across the five-way slice');style_axes(ax);
exportgraphics(fig,fullfile(outdir,'five_way_late_distance_heatmap.png'),'Resolution',300);
close(fig);

lifetimes=[records.decision_time];
lifetimes(~isfinite(lifetimes))=[records(~isfinite(lifetimes)).duration];
fig=figure('Color','w','Position',[100 100 950 500]);ax=axes(fig);hold(ax,'on');
for k=1:numel(records)
    scatter(ax,records(k).lambda,lifetimes(k),65,outcome_color(records(k).outcome), ...
        'filled','MarkerEdgeColor','k');
end
xlabel(ax,'lambda');ylabel(ax,'decision time or tested duration');
title(ax,sprintf('Outcome commitment times, Delta lambda = %.3f',cfg.scan.lambda_step));
style_axes(ax);exportgraphics(fig,fullfile(outdir,'five_way_decision_time.png'), ...
    'Resolution',300);close(fig);
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
