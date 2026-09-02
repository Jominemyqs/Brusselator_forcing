function results = run_dynamic_edge_tracking_BC()
%RUN_DYNAMIC_EDGE_TRACKING_BC One-event dynamic B--C edge-tracking pilot.
%   Evolves a close B/C pair until it separates, rebisects the complete
%   evolved states using the physical five-way outcome classifier, and then
%   evolves the renewed pair to test whether boundary shadowing is prolonged.

cfg=brusselator_dynamic_edge_BC_config();outdir=fullfile(cfg.output.root, ...
    cfg.experiment_name);checkpoint=fullfile(outdir,'progress_dynamic_edge.mat');
final_file=fullfile(outdir,'dynamic_edge_tracking_BC.mat');
if isfile(final_file)
    error('run_dynamic_edge_tracking_BC:OutputExists', ...
        'Refusing to overwrite completed result: %s',final_file);
end
[source_B,source_C,library]=load_sources(cfg);validate_sources(source_B,source_C,cfg);
if isfile(checkpoint)
    saved=load(checkpoint,'tracker','records','raw_files');
    tracker=saved.tracker;records=saved.records;raw_files=saved.raw_files;
    fprintf('Resuming %s in phase %s.\n',cfg.experiment_name,tracker.phase);
else
    if isfolder(outdir)
        error('run_dynamic_edge_tracking_BC:IncompleteOutput', ...
            'Existing output lacks a valid checkpoint: %s',outdir);
    end
    mkdir(outdir);mkdir(fullfile(outdir,'raw_cases'));mkdir(fullfile(outdir,'segments'));
    template=empty_record();records=template([]);raw_files={};
    tracker=initialize_tracker(source_B,source_C,cfg,outdir,library);
    brusselator_write_metadata(outdir,cfg, ...
        brusselator_run_metadata(cfg,mfilename,NaN));
    write_checkpoint(checkpoint,outdir,cfg,tracker,records,raw_files);
end
timer=tic;

while strcmp(tracker.phase,'rebisect')
    if tracker.rebisection_step>cfg.dynamic_edge.maximum_rebisection_steps
        tracker.status='maximum_rebisection_steps_reached';
        tracker.phase='complete';break;
    end
    step=tracker.rebisection_step;
    alpha=0.5*(tracker.alpha_B+tracker.alpha_C);
    midpoint=interpolate_states(tracker.state_B,tracker.state_C,0.5);
    case_id=sprintf('event01_step%02d',step);
    raw_file=fullfile(outdir,'raw_cases',[case_id,'.mat']);
    before_B=tracker.alpha_B;before_C=tracker.alpha_C;
    fprintf('Dynamic B--C rebisection %d: alpha=%.10f, separation=%.4e.\n', ...
        step,alpha,pair_distance(tracker.state_B,tracker.state_C,cfg.grid.Lx));
    if isfile(raw_file)
        loaded=load(raw_file,'edge_case');edge_case=loaded.edge_case;
        fprintf('  Recovered saved case %s.\n',case_id);
    else
        edge_case=brusselator_evolve_and_classify_state( ...
            case_id,midpoint,library,cfg);
        save(raw_file,'edge_case','-v7.3');
    end
    outcome=edge_case.classification.outcome;update='obstruction';
    if strcmp(outcome,cfg.dynamic_edge.left_outcome)
        tracker.state_B=midpoint;tracker.alpha_B=alpha;update='B_endpoint';
    elseif strcmp(outcome,cfg.dynamic_edge.right_outcome)
        tracker.state_C=midpoint;tracker.alpha_C=alpha;update='C_endpoint';
    else
        tracker.status=['rebisection_obstructed_by_',outcome];
        tracker.phase='complete';
    end
    separation=pair_distance(tracker.state_B,tracker.state_C,cfg.grid.Lx);
    record=summarize_case(edge_case,step,alpha,before_B,before_C, ...
        tracker.alpha_B,tracker.alpha_C,separation,update,raw_file);
    records(end+1,1)=record;raw_files{end+1,1}=raw_file; %#ok<AGROW>
    tracker.rebisection_step=step+1;tracker.last_case_id=case_id;
    if strcmp(tracker.phase,'rebisect') && ...
            separation<=cfg.dynamic_edge.restart_separation_target
        tracker.status='first_dynamic_rebracketing_complete';
        if cfg.dynamic_edge.run_post_rebracketing_segment
            tracker.phase='post_segment';
        else
            tracker.phase='complete';
        end
    end
    write_checkpoint(checkpoint,outdir,cfg,tracker,records,raw_files);
end

if strcmp(tracker.phase,'post_segment')
    [tracker,segment_record]=run_post_segment(tracker,cfg,outdir,library);
    tracker.segment_records(end+1,1)=segment_record;
    tracker.phase='complete';tracker.status='one_rebracketing_cycle_complete';
    write_checkpoint(checkpoint,outdir,cfg,tracker,records,raw_files);
end

summary=make_summary(tracker,records,cfg);
writetable(struct2table(records,'AsArray',true), ...
    fullfile(outdir,'dynamic_rebisection_cases.csv'));
writetable(struct2table(tracker.segment_records,'AsArray',true), ...
    fullfile(outdir,'dynamic_edge_segments.csv'));
writetable(struct2table(summary,'AsArray',true), ...
    fullfile(outdir,'dynamic_edge_summary.csv'));
make_figures(outdir,tracker,records);
results=struct('configuration',cfg,'tracker',tracker,'records',records, ...
    'raw_files',{raw_files},'summary',summary);
save(final_file,'results','-v7.3');
brusselator_write_metadata(outdir,cfg, ...
    brusselator_run_metadata(cfg,mfilename,toc(timer)));
disp(struct2table(records,'AsArray',true));disp(struct2table(summary,'AsArray',true));
fprintf('Dynamic B--C edge pilot saved in: %s\n',outdir);
end

function [B,C,library]=load_sources(cfg)
if ~isfile(cfg.dynamic_edge.source_B_file)||~isfile(cfg.dynamic_edge.source_C_file)|| ...
        ~isfile(cfg.dynamic_edge.parent_library_file)
    error('run_dynamic_edge_tracking_BC:MissingSource','A configured source is missing.');
end
loaded=load(cfg.dynamic_edge.source_B_file,'refinement_case');B=loaded.refinement_case;
loaded=load(cfg.dynamic_edge.source_C_file,'results');C=loaded.results;
loaded=load(cfg.dynamic_edge.parent_library_file,'results');library=loaded.results.library;
end

function validate_sources(B,C,cfg)
if ~strcmp(B.classification.outcome,cfg.dynamic_edge.left_outcome)|| ...
        ~strcmp(C.classification.outcome,cfg.dynamic_edge.right_outcome)|| ...
        abs(B.lambda-cfg.dynamic_edge.initial_lambda_B)>1e-12|| ...
        abs(C.summary.lambda-cfg.dynamic_edge.initial_lambda_C)>1e-12
    error('run_dynamic_edge_tracking_BC:SourceLabels', ...
        'The saved initial pair does not match the configured B/C bracket.');
end
if max(abs(B.initial_state(:)-C.initial_state(:)))==0
    error('run_dynamic_edge_tracking_BC:IdenticalSources','B and C source states are identical.');
end
direct=pair_distance(B.initial_state,C.initial_state,cfg.grid.Lx);
reflected=pair_reflected_distance(B.initial_state,C.initial_state,cfg.grid.Lx);
if direct>=cfg.dynamic_edge.separation_threshold||direct>=reflected
    error('run_dynamic_edge_tracking_BC:SourceGeometry', ...
        'Initial pair is not close in its maintained common orientation.');
end
end

function tracker=initialize_tracker(B,C,cfg,outdir,library)
[segment_data,segment_record,state_B,state_C]=source_segment(B,C,cfg,library);
segment_file=fullfile(outdir,'segments','segment_00_source.mat');
save(segment_file,'segment_data','-v7.3');segment_record.raw_file=segment_file;
tracker=struct('phase','rebisect','status','source_pair_reached_divergence', ...
    'global_time',segment_record.global_end_time,'state_B',state_B, ...
    'state_C',state_C,'alpha_B',0,'alpha_C',1,'rebisection_step',1, ...
    'last_case_id','','segment_records',segment_record, ...
    'initial_pair_separation',segment_record.initial_separation, ...
    'initial_divergence_time',segment_record.global_end_time, ...
    'post_rebracketing_divergence_time',NaN);
end

function [data,record,state_B,state_C]=source_segment(B,C,cfg,library)
n=min(numel(B.S),numel(C.S));
if max(abs(B.S(1:n)-C.S(1:n)))>1e-12
    error('run_dynamic_edge_tracking_BC:SourceTimes','Source time grids differ.');
end
S=B.S(1:n);UB=B.U(1:n,:);VB=B.V(1:n,:);UC=C.U(1:n,:);VC=C.V(1:n,:);
[separation,reflected]=pair_distance_history(UB,VB,UC,VC,cfg.grid.Lx);
index=find(separation>=cfg.dynamic_edge.separation_threshold,1);
if isempty(index)
    error('run_dynamic_edge_tracking_BC:NoSourceDivergence', ...
        'Saved source pair never reaches the configured separation threshold.');
end
[data,record]=package_segment(0,'reused_saved_source',S(1:index), ...
    UB(1:index,:),VB(1:index,:),UC(1:index,:),VC(1:index,:), ...
    separation(1:index),reflected(1:index),0,cfg,library);
x=B.initial_state(:,1);state_B=[x,UB(index,:)',VB(index,:)'];
state_C=[x,UC(index,:)',VC(index,:)'];
end

function [tracker,record]=run_post_segment(tracker,cfg,outdir,library)
par=brusselator_make_parameters(cfg,@(t)cfg.edge_tracking.frozen_b);
fprintf('Evolving renewed B/C pair from global edge time %.3f.\n',tracker.global_time);
[~,SB,VB,UB]=solve_brusselator_1d_forced(tracker.state_B,par, ...
    cfg.dynamic_edge.segment_maximum_duration,0);
[~,SC,VC,UC]=solve_brusselator_1d_forced(tracker.state_C,par, ...
    cfg.dynamic_edge.segment_maximum_duration,0);
if numel(SB)~=numel(SC)||max(abs(SB-SC))>1e-12
    error('run_dynamic_edge_tracking_BC:PostTimes','Post-segment time grids differ.');
end
[separation,reflected]=pair_distance_history(UB,VB,UC,VC,cfg.grid.Lx);
index=find(separation>=cfg.dynamic_edge.separation_threshold,1);
threshold_reached=~isempty(index);if ~threshold_reached,index=numel(SB);end
[segment_data,record]=package_segment(1,'post_rebracketing',SB(1:index), ...
    UB(1:index,:),VB(1:index,:),UC(1:index,:),VC(1:index,:), ...
    separation(1:index),reflected(1:index),tracker.global_time,cfg,library);
record.threshold_reached=threshold_reached;
segment_file=fullfile(outdir,'segments','segment_01_post_rebracketing.mat');
save(segment_file,'segment_data','-v7.3');record.raw_file=segment_file;
tracker.global_time=record.global_end_time;
if threshold_reached,tracker.post_rebracketing_divergence_time=record.duration;end
end

function [data,record]=package_segment(index,source,S,UB,VB,UC,VC, ...
        separation,reflected,global_start,cfg,library)
shadow_U=0.5*(UB+UC);shadow_V=0.5*(VB+VC);
[shadow_classification,shadow_history]=brusselator_classify_five_way_trajectory( ...
    S,shadow_U,shadow_V,library,'late_window',min(cfg.edge_tracking.late_window,S(end)), ...
    'sample_interval',cfg.edge_tracking.distance_sample_interval, ...
    'median_threshold',cfg.edge_tracking.median_distance_threshold, ...
    'maximum_threshold',cfg.edge_tracking.maximum_distance_threshold,'Lx',cfg.grid.Lx);
data=struct('segment_index',index,'source',source,'S',S,'U_B',UB,'V_B',VB, ...
    'U_C',UC,'V_C',VC,'shadow_U',shadow_U,'shadow_V',shadow_V, ...
    'separation',separation,'reflected_separation',reflected, ...
    'shadow_classification_diagnostic',shadow_classification, ...
    'shadow_distance_history',shadow_history, ...
    'warning',['the arithmetic midpoint shadow is not an exact PDE trajectory; ', ...
        'its orbit distances are visualization diagnostics only']);
record=empty_segment_record();record.segment_index=index;record.source=source;
record.global_start_time=global_start;record.global_end_time=global_start+S(end);
record.duration=S(end);record.initial_separation=separation(1);
record.final_separation=separation(end);record.maximum_separation=max(separation);
record.initial_reflected_separation=reflected(1);
record.minimum_reflected_separation=min(reflected);
record.threshold_reached=separation(end)>=cfg.dynamic_edge.separation_threshold;
record.raw_file='';
end

function record=empty_segment_record()
record=struct('segment_index',NaN,'source','','global_start_time',NaN, ...
    'global_end_time',NaN,'duration',NaN,'initial_separation',NaN, ...
    'final_separation',NaN,'maximum_separation',NaN, ...
    'initial_reflected_separation',NaN,'minimum_reflected_separation',NaN, ...
    'threshold_reached',false,'raw_file','');
end

function state=interpolate_states(B,C,alpha)
state=B;state(:,2:3)=(1-alpha)*B(:,2:3)+alpha*C(:,2:3);
end

function record=empty_record()
record=struct('event_index',1,'step',NaN,'case_id','','alpha',NaN, ...
    'outcome','','tested_duration',NaN,'decision_time',NaN, ...
    'runtime_seconds',NaN,'alpha_B_before',NaN,'alpha_C_before',NaN, ...
    'alpha_B_after',NaN,'alpha_C_after',NaN,'state_separation_after',NaN, ...
    'update','','late_median_distance_A',NaN,'late_median_distance_B',NaN, ...
    'late_median_distance_C',NaN,'late_median_distance_R5',NaN, ...
    'late_median_distance_R6',NaN,'recurrence_available',false, ...
    'recurrence_period',NaN,'recurrence_distance',NaN,'raw_file','');
end

function record=summarize_case(data,step,alpha,bB,bC,aB,aC,separation,update,raw)
record=empty_record();record.step=step;record.case_id=data.case_id;
record.alpha=alpha;record.outcome=data.classification.outcome;
record.tested_duration=data.duration;record.decision_time=data.decision_time;
record.runtime_seconds=data.runtime_seconds;record.alpha_B_before=bB;
record.alpha_C_before=bC;record.alpha_B_after=aB;record.alpha_C_after=aC;
record.state_separation_after=separation;record.update=update;
labels={'A','B','C','R5','R6'};
for k=1:numel(labels)
    record.(['late_median_distance_',labels{k}])= ...
        data.classification.late_median_distances(k);
end
record.recurrence_available=data.recurrence.available;
record.recurrence_period=data.recurrence.period;
record.recurrence_distance=data.recurrence.distance;record.raw_file=raw;
end

function summary=make_summary(tracker,records,cfg)
summary=struct('status',tracker.status,'initial_pair_separation', ...
    tracker.initial_pair_separation,'separation_threshold', ...
    cfg.dynamic_edge.separation_threshold,'initial_divergence_time', ...
    tracker.initial_divergence_time,'rebisection_steps',numel(records), ...
    'renewed_pair_separation',pair_distance(tracker.state_B,tracker.state_C, ...
    cfg.grid.Lx),'post_rebracketing_divergence_time', ...
    tracker.post_rebracketing_divergence_time,'total_shadow_time', ...
    tracker.global_time,'interpretation',cfg.dynamic_edge.interpretation);
end

function write_checkpoint(file,outdir,cfg,tracker,records,raw_files)
save(file,'cfg','tracker','records','raw_files','-v7.3');
if ~isempty(records)
    writetable(struct2table(records,'AsArray',true), ...
        fullfile(outdir,'progress_rebisection.csv'));
end
end

function make_figures(outdir,tracker,records)
fig=figure('Color','w','Position',[100 100 1000 560]);ax=axes(fig);hold(ax,'on');
for k=1:numel(tracker.segment_records)
    loaded=load(tracker.segment_records(k).raw_file,'segment_data');d=loaded.segment_data;
    plot(ax,tracker.segment_records(k).global_start_time+d.S,d.separation, ...
        'LineWidth',1.5,'DisplayName',sprintf('segment %d',k-1));
end
yline(ax,0.02,'--k','rebracketing threshold','HandleVisibility','off');
set(ax,'YScale','log');xlabel(ax,'assembled shadow time');
ylabel(ax,'B/C bracket separation');title(ax,'Dynamic B--C edge-shadow segments');
style_axes(ax);legend(ax,'Location','best');
exportgraphics(fig,fullfile(outdir,'dynamic_edge_separation.png'),'Resolution',300);
close(fig);

fig=figure('Color','w','Position',[100 100 950 520]);ax=axes(fig);hold(ax,'on');
for k=1:numel(records)
    scatter(ax,records(k).step,records(k).state_separation_after,75, ...
        outcome_color(records(k).outcome),'filled','MarkerEdgeColor','k');
end
yline(ax,2.5e-4,'--k','restart target','HandleVisibility','off');
set(ax,'YScale','log');xlabel(ax,'rebisection step');
ylabel(ax,'renewed bracket separation');title(ax,'Multiclass-safe dynamic rebisection');
style_axes(ax);exportgraphics(fig,fullfile(outdir,'dynamic_rebisection.png'), ...
    'Resolution',300);close(fig);
end

function [distance,reflected]=pair_distance_history(UB,VB,UC,VC,Lx)
N=size(UB,2);dx=Lx/(N-1);w=dx*ones(1,N);w([1,end])=0.5*dx;
numerator=sum(((UB-UC).^2+(VB-VC).^2).*w,2);
energy_B=sum((UB.^2+VB.^2).*w,2);energy_C=sum((UC.^2+VC.^2).*w,2);
distance=sqrt(numerator./max(0.5*(energy_B+energy_C),eps));
ref_U=fliplr(UC);ref_V=fliplr(VC);
ref_num=sum(((UB-ref_U).^2+(VB-ref_V).^2).*w,2);
reflected=sqrt(ref_num./max(0.5*(energy_B+energy_C),eps));
end

function distance=pair_distance(B,C,Lx)
[distance,~]=pair_distance_history(B(:,2)',B(:,3)',C(:,2)',C(:,3)',Lx);
end

function distance=pair_reflected_distance(B,C,Lx)
[~,distance]=pair_distance_history(B(:,2)',B(:,3)',C(:,2)',C(:,3)',Lx);
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
