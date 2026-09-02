function results = run_ml_evolved_BA_five_way_scan()
%RUN_ML_EVOLVED_BA_FIVE_WAY_SCAN Map the C obstruction on an evolved BA chord.

cfg=brusselator_ml_evolved_BA_scan_config();
outdir=fullfile(cfg.output.root,cfg.experiment_name);
checkpoint=fullfile(outdir,'progress_evolved_BA_scan.mat');
final_file=fullfile(outdir,'evolved_BA_scan.mat');
if isfile(final_file)
    error('run_ml_evolved_BA_five_way_scan:OutputExists', ...
        'Refusing to overwrite completed scan: %s',final_file);
end
[left,right,reuse_cases,library,parent]=load_sources(cfg);
alphas=cfg.evolved_scan.alpha_values(:);validate_sources(left,right,reuse_cases,cfg);
if isfile(checkpoint)
    saved=load(checkpoint,'records','raw_files');
    records=saved.records;raw_files=saved.raw_files;
    if numel(records)~=numel(alphas)||max(abs([records.alpha]'-alphas))>1e-12
        error('run_ml_evolved_BA_five_way_scan:CheckpointMismatch', ...
            'Checkpoint alpha grid does not match the configuration.');
    end
    fprintf('Resuming %s with %d/%d cases complete.\n', ...
        cfg.experiment_name,sum([records.completed]),numel(records));
else
    if isfolder(outdir)
        error('run_ml_evolved_BA_five_way_scan:UnsafeOutput', ...
            'Existing output has no recognized checkpoint: %s',outdir);
    end
    mkdir(outdir);mkdir(fullfile(outdir,'raw_cases'));
    records=repmat(empty_record(),numel(alphas),1);raw_files=cell(numel(alphas),1);
    for k=1:numel(alphas),records(k).scan_index=k;records(k).alpha=alphas(k);end
    brusselator_write_metadata(outdir,cfg, ...
        brusselator_run_metadata(cfg,mfilename,NaN));
    write_progress();
end
timer=tic;
for k=1:numel(alphas)
    if records(k).completed,continue;end
    alpha=alphas(k);initial_state=interpolate_states(left,right,alpha);
    [reuse,data,source_type,source_file]=reuse_case(alpha,initial_state, ...
        left,right,reuse_cases,cfg);
    if ~reuse
        case_id=sprintf('evolved_BA_alpha_%06.3f',alpha);
        source_file=fullfile(outdir,'raw_cases',[strrep(case_id,'.','p'),'.mat']);
        partial_file=fullfile(outdir,'raw_cases',[strrep(case_id,'.','p'),'_partial.mat']);
        if isfile(source_file)
            loaded=load(source_file,'scan_case');data=loaded.scan_case;
        else
            data=brusselator_evolve_and_classify_state( ...
                case_id,initial_state,library,cfg,partial_file);
            scan_case=data;save(source_file,'scan_case','-v7.3');
            if isfile(partial_file),delete(partial_file);end
        end
        source_type='new_frozen_PDE_label';
    end
    records(k)=summarize_case(records(k),data,source_type,source_file,library);
    raw_files{k}=source_file;write_progress();
end
transitions=find_transitions(records);regions=find_regions(records);
summary=make_summary(records,transitions,regions,cfg);
writetable(struct2table(records,'AsArray',true), ...
    fullfile(outdir,'evolved_BA_five_way_scan.csv'));
writetable(struct2table(transitions,'AsArray',true), ...
    fullfile(outdir,'candidate_transition_intervals.csv'));
writetable(struct2table(regions,'AsArray',true), ...
    fullfile(outdir,'contiguous_sampled_regions.csv'));
writetable(struct2table(summary,'AsArray',true), ...
    fullfile(outdir,'evolved_BA_scan_summary.csv'));
make_figures(outdir,records,cfg);
results=struct('configuration',cfg,'parent_result_file', ...
    cfg.evolved_scan.parent_result_file,'parent_global_time',parent.tracker.global_time, ...
    'endpoint_left',left,'endpoint_right',right,'records',records, ...
    'raw_files',{raw_files},'transition_records',transitions, ...
    'region_records',regions,'summary',summary);
save(final_file,'results','-v7.3');
brusselator_write_metadata(outdir,cfg, ...
    brusselator_run_metadata(cfg,mfilename,toc(timer)));
disp(struct2table(records,'AsArray',true));
disp(struct2table(transitions,'AsArray',true));
disp(struct2table(summary,'AsArray',true));
fprintf('Evolved B--A five-way scan saved in: %s\n',outdir);

    function write_progress()
        save(checkpoint,'cfg','records','raw_files','-v7.3');
        completed=records([records.completed]);
        if ~isempty(completed)
            writetable(struct2table(completed,'AsArray',true), ...
                fullfile(outdir,'progress_summary.csv'));
        end
    end
end

function [left,right,reuse,library,parent]=load_sources(cfg)
if ~isfile(cfg.evolved_scan.parent_result_file)|| ...
        ~isfile(cfg.evolved_scan.reuse_result_file)
    error('Configured evolved-edge source result is missing.');
end
loaded=load(cfg.evolved_scan.parent_result_file,'results');parent=loaded.results;
loaded=load(cfg.evolved_scan.reuse_result_file,'results');reuse=loaded.results;
left=parent.tracker.state_left;right=parent.tracker.state_right;
[library,~]=brusselator_build_ml_outcome_library(cfg);
end

function validate_sources(left,right,reuse,cfg)
if size(left,1)~=cfg.grid.N||size(right,1)~=cfg.grid.N|| ...
        max(abs(left(:,1)-right(:,1)))>1e-12
    error('Evolved endpoint grids do not match the configured grid.');
end
expected=[0.5,0.75,0.875];actual=[reuse.case_records.alpha];
if numel(actual)~=3||max(abs(actual-expected))>1e-12|| ...
        ~isequal({reuse.case_records.outcome},{'B','B','C'})
    error('Saved obstruction cases are not the expected B, B, C sequence.');
end
for k=1:numel(actual)
    loaded=load(reuse.raw_files{k},'edge_case');
    expected_state=interpolate_states(left,right,actual(k));
    if max(abs(loaded.edge_case.initial_state(:)-expected_state(:)))> ...
            cfg.evolved_scan.reuse_tolerance
        error('Reusable alpha %.3f does not lie on the configured chord.',actual(k));
    end
end
end

function [found,data,source_type,source_file]=reuse_case( ...
        alpha,state,left,right,reuse,cfg)
found=true;data=[];source_type='';source_file='';
if abs(alpha)<1e-12
    data=endpoint_case('evolved_BA_alpha_000p000',state,'B');
    source_type='invariant_propagated_B_endpoint';source_file=cfg.evolved_scan.parent_result_file;
elseif abs(alpha-1)<1e-12
    data=endpoint_case('evolved_BA_alpha_001p000',state,'A');
    source_type='invariant_propagated_A_endpoint';source_file=cfg.evolved_scan.parent_result_file;
else
    index=find(abs([reuse.case_records.alpha]-alpha)<1e-12,1);
    if isempty(index),found=false;return;end
    source_file=reuse.raw_files{index};loaded=load(source_file,'edge_case');data=loaded.edge_case;
    if max(abs(data.initial_state(:)-state(:)))>cfg.evolved_scan.reuse_tolerance
        error('Reusable state mismatch at alpha %.6f.',alpha);
    end
    source_type='reused_obstruction_label';
end
if alpha==0&&max(abs(state(:)-left(:)))>1e-12,error('Left endpoint mismatch.');end
if alpha==1&&max(abs(state(:)-right(:)))>1e-12,error('Right endpoint mismatch.');end
end

function data=endpoint_case(case_id,state,outcome)
classification=struct('outcome',outcome,'late_median_distances',nan(1,5), ...
    'late_maximum_distances',nan(1,5),'final_distances',nan(1,5));
data=struct('case_id',case_id,'initial_state',state,'duration',0, ...
    'runtime_seconds',0,'decision_time',0,'classification',classification);
end

function state=interpolate_states(left,right,alpha)
state=left;state(:,2:3)=(1-alpha)*left(:,2:3)+alpha*right(:,2:3);
end

function record=empty_record()
record=struct('scan_index',NaN,'alpha',NaN,'completed',false, ...
    'outcome','','source_type','','source_file','','duration',NaN, ...
    'decision_time',NaN,'runtime_seconds',NaN, ...
    'late_distance_A',NaN,'late_distance_B',NaN,'late_distance_C',NaN, ...
    'late_distance_R5',NaN,'late_distance_R6',NaN);
end

function record=summarize_case(record,data,source_type,source_file,library)
record.completed=true;record.outcome=data.classification.outcome;
record.source_type=source_type;record.source_file=source_file;
record.duration=data.duration;record.decision_time=data.decision_time;
record.runtime_seconds=data.runtime_seconds;
for j=1:numel(library)
    record.(['late_distance_',library(j).label])= ...
        data.classification.late_median_distances(j);
end
end

function records=find_transitions(scan)
template=struct('left_alpha',NaN,'right_alpha',NaN,'width',NaN, ...
    'left_outcome','','right_outcome','','pair','','status','');records=template([]);
known={'A','B','C','R5','R6'};
for k=1:numel(scan)-1
    if strcmp(scan(k).outcome,scan(k+1).outcome),continue;end
    item=template;item.left_alpha=scan(k).alpha;item.right_alpha=scan(k+1).alpha;
    item.width=item.right_alpha-item.left_alpha;item.left_outcome=scan(k).outcome;
    item.right_outcome=scan(k+1).outcome;item.pair=[item.left_outcome,'--',item.right_outcome];
    if ismember(item.left_outcome,known)&&ismember(item.right_outcome,known)
        item.status='candidate_adjacent_at_sample_resolution';
    else
        item.status='contains_unresolved_or_ambiguous_endpoint';
    end
    records(end+1,1)=item; %#ok<AGROW>
end
end

function regions=find_regions(scan)
template=struct('outcome','','alpha_minimum',NaN,'alpha_maximum',NaN, ...
    'sample_count',NaN);regions=template([]);start=1;
for k=2:numel(scan)+1
    if k<=numel(scan)&&strcmp(scan(k).outcome,scan(start).outcome),continue;end
    item=template;item.outcome=scan(start).outcome;
    item.alpha_minimum=scan(start).alpha;item.alpha_maximum=scan(k-1).alpha;
    item.sample_count=k-start;regions(end+1,1)=item; %#ok<AGROW>
    start=k;
end
end

function summary=make_summary(records,transitions,regions,cfg)
labels={records.outcome};new_count=sum(strcmp({records.source_type},'new_frozen_PDE_label'));
summary=struct('sample_count',numel(records),'new_PDE_label_count',new_count, ...
    'reused_count',numel(records)-new_count,'number_A',sum(strcmp(labels,'A')), ...
    'number_B',sum(strcmp(labels,'B')),'number_C',sum(strcmp(labels,'C')), ...
    'number_R5',sum(strcmp(labels,'R5')),'number_R6',sum(strcmp(labels,'R6')), ...
    'number_unresolved',sum(strcmp(labels,'unresolved')), ...
    'contiguous_region_count',numel(regions),'transition_interval_count', ...
    numel(transitions),'interpretation',cfg.evolved_scan.transition_rule);
end

function make_figures(outdir,records,cfg)
order={'B','R5','C','R6','A','unresolved','ambiguous_multiple_outcomes'};
labels={records.outcome};y=zeros(size(labels));
for k=1:numel(labels)
    index=find(strcmp(order,labels{k}),1);if isempty(index),index=numel(order);end
    y(k)=index;
end
fig=figure('Visible','off','Color','w','Position',[100 100 1050 520]);
ax=axes(fig);hold(ax,'on');
for k=1:numel(records)
    scatter(ax,records(k).alpha,y(k),75,outcome_color(records(k).outcome), ...
        'filled','MarkerEdgeColor','k');
end
plot(ax,[records.alpha],y,'-k','LineWidth',0.8);
yticks(ax,1:numel(order));yticklabels(ax,order);grid(ax,'on');box(ax,'on');
xlabel(ax,'alpha on evolved B--A chord');ylabel(ax,'five-way outcome');
title(ax,'Multiclass geometry after the first B--A divergence');
exportgraphics(fig,fullfile(outdir,'evolved_BA_outcome_scan.png'),'Resolution',300);close(fig);

times=[records.decision_time];times(~isfinite(times))=[records(~isfinite(times)).duration];
fig=figure('Visible','off','Color','w','Position',[100 100 1050 520]);
ax=axes(fig);hold(ax,'on');
for k=1:numel(records)
    scatter(ax,records(k).alpha,times(k),70,outcome_color(records(k).outcome), ...
        'filled','MarkerEdgeColor','k');
end
grid(ax,'on');box(ax,'on');xlabel(ax,'alpha on evolved B--A chord');
ylabel(ax,'decision time');title(ax,sprintf( ...
    'Outcome commitment; focused spacing %.3f for alpha >= %.2f', ...
    cfg.evolved_scan.focused_spacing,cfg.evolved_scan.focused_minimum));
exportgraphics(fig,fullfile(outdir,'evolved_BA_decision_times.png'),'Resolution',300);close(fig);
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
