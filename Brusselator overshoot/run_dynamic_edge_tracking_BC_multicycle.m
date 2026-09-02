function results = run_dynamic_edge_tracking_BC_multicycle()
%RUN_DYNAMIC_EDGE_TRACKING_BC_MULTICYCLE Add three B--C rebracketing cycles.
%   Continues the completed one-event pilot. Each renewed segment yields a
%   Poincare-phase-fixed template, recurrence period, consecutive-template
%   distance, and direct one-period flow residual for the Newton-seed gate.

cfg=brusselator_dynamic_edge_multicycle_config();outdir=fullfile(cfg.output.root, ...
    cfg.experiment_name);checkpoint=fullfile(outdir,'progress_multicycle.mat');
final_file=fullfile(outdir,'dynamic_edge_tracking_BC_multicycle.mat');
if isfile(final_file)
    error('run_dynamic_edge_tracking_BC_multicycle:OutputExists', ...
        'Refusing to overwrite completed result: %s',final_file);
end
[parent,parent_segment,library]=load_sources(cfg);validate_parent(parent,parent_segment,cfg);
if isfile(checkpoint)
    saved=load(checkpoint,'tracker','case_records','raw_files', ...
        'event_records','template_records','template_states');
    tracker=saved.tracker;case_records=saved.case_records;raw_files=saved.raw_files;
    event_records=saved.event_records;template_records=saved.template_records;
    template_states=saved.template_states;
    fprintf('Resuming multicycle edge tracking at event %d, phase %s.\n', ...
        tracker.event_index,tracker.phase);
else
    if isfolder(outdir)
        error('run_dynamic_edge_tracking_BC_multicycle:IncompleteOutput', ...
            'Existing output lacks a valid checkpoint: %s',outdir);
    end
    mkdir(outdir);mkdir(fullfile(outdir,'raw_cases'));mkdir(fullfile(outdir,'segments'));
    case_template=empty_case_record();case_records=case_template([]);raw_files={};
    event_template=empty_event_record();event_records=event_template([]);
    [parent_template,parent_state]=extract_template(parent_segment,1,cfg,library);
    template_records=parent_template;template_states={parent_state};
    tracker=initialize_tracker(parent,parent_segment,cfg);
    brusselator_write_metadata(outdir,cfg, ...
        brusselator_run_metadata(cfg,mfilename,NaN));
    write_checkpoint(checkpoint,outdir,cfg,tracker,case_records,raw_files, ...
        event_records,template_records,template_states);
end
timer=tic;last_event=cfg.multicycle.first_new_event_index+ ...
    cfg.multicycle.additional_rebracketing_events-1;

while tracker.event_index<=last_event && ~strcmp(tracker.phase,'complete')
    if strcmp(tracker.phase,'rebisect')
        [tracker,case_records,raw_files]=rebisect_event(tracker,case_records, ...
            raw_files,library,cfg,outdir,checkpoint,event_records, ...
            template_records,template_states);
    end
    if strcmp(tracker.phase,'evolve_segment')
        [tracker,segment_data,event_record]=evolve_segment(tracker,cfg,library);
        segment_file=fullfile(outdir,'segments',sprintf('segment_%02d.mat', ...
            tracker.event_index));save(segment_file,'segment_data','-v7.3');
        event_record.raw_file=segment_file;event_records(end+1,1)=event_record; %#ok<AGROW>
        [template_record,template_state]=extract_template( ...
            segment_data,tracker.event_index,cfg,library);
        [distance,reflected,uses_reflection]=template_distance( ...
            template_states{end},template_state,cfg.grid.Lx);
        template_record.previous_template_distance=distance;
        template_record.previous_template_direct_distance=reflected(1);
        template_record.previous_template_reflected_distance=reflected(2);
        template_record.previous_template_uses_reflection=uses_reflection;
        template_records(end+1,1)=template_record; %#ok<AGROW>
        template_states{end+1,1}=template_state; %#ok<AGROW>
        tracker.global_time=event_record.global_end_time;
        tracker.state_B=state_from_row(segment_data.U_B(end,:), ...
            segment_data.V_B(end,:),tracker.state_B(:,1));
        tracker.state_C=state_from_row(segment_data.U_C(end,:), ...
            segment_data.V_C(end,:),tracker.state_C(:,1));
        if event_record.threshold_reached
            tracker.event_index=tracker.event_index+1;tracker.rebisection_step=1;
            tracker.alpha_B=0;tracker.alpha_C=1;
        end
        if ~event_record.threshold_reached
            tracker.phase='complete';
        elseif tracker.event_index>last_event
            tracker.phase='complete';tracker.status='requested_multicycle_run_complete';
        else
            tracker.phase='rebisect';tracker.status='next_divergence_ready';
        end
        write_checkpoint(checkpoint,outdir,cfg,tracker,case_records,raw_files, ...
            event_records,template_records,template_states);
    end
end

summary=make_summary(parent,tracker,event_records,template_records,cfg);
writetable(struct2table(case_records,'AsArray',true), ...
    fullfile(outdir,'multicycle_rebisection_cases.csv'));
writetable(struct2table(event_records,'AsArray',true), ...
    fullfile(outdir,'multicycle_events.csv'));
writetable(struct2table(template_records,'AsArray',true), ...
    fullfile(outdir,'phase_template_convergence.csv'));
writetable(struct2table(summary,'AsArray',true), ...
    fullfile(outdir,'multicycle_summary.csv'));
make_figures(outdir,parent,event_records,template_records,cfg);
results=struct('configuration',cfg,'parent_result_file', ...
    cfg.multicycle.parent_result_file,'tracker',tracker, ...
    'case_records',case_records,'raw_files',{raw_files}, ...
    'event_records',event_records,'template_records',template_records, ...
    'template_states',{template_states},'summary',summary);
save(final_file,'results','-v7.3');
brusselator_write_metadata(outdir,cfg, ...
    brusselator_run_metadata(cfg,mfilename,toc(timer)));
disp(struct2table(event_records,'AsArray',true));
disp(struct2table(template_records,'AsArray',true));disp(struct2table(summary,'AsArray',true));
fprintf('Multicycle B--C edge tracking saved in: %s\n',outdir);
end

function [parent,segment,library]=load_sources(cfg)
if ~isfile(cfg.multicycle.parent_result_file)|| ...
        ~isfile(cfg.multicycle.parent_post_segment_file)|| ...
        ~isfile(cfg.dynamic_edge.parent_library_file)
    error('run_dynamic_edge_tracking_BC_multicycle:MissingSource', ...
        'A configured parent source is missing.');
end
loaded=load(cfg.multicycle.parent_result_file,'results');parent=loaded.results;
loaded=load(cfg.multicycle.parent_post_segment_file,'segment_data');segment=loaded.segment_data;
loaded=load(cfg.dynamic_edge.parent_library_file,'results');library=loaded.results.library;
end

function validate_parent(parent,segment,cfg)
if ~strcmp(parent.summary.status,'one_rebracketing_cycle_complete')|| ...
        segment.segment_index~=1|| ...
        segment.separation(end)<cfg.dynamic_edge.separation_threshold
    error('run_dynamic_edge_tracking_BC_multicycle:ParentMismatch', ...
        'The configured parent is not the completed one-event B--C pilot.');
end
end

function tracker=initialize_tracker(parent,segment,cfg)
x=parent.tracker.state_B(:,1);
tracker=struct('phase','rebisect','status','parent_divergence_ready', ...
    'event_index',cfg.multicycle.first_new_event_index, ...
    'rebisection_step',1,'global_time',parent.summary.total_shadow_time, ...
    'state_B',state_from_row(segment.U_B(end,:),segment.V_B(end,:),x), ...
    'state_C',state_from_row(segment.U_C(end,:),segment.V_C(end,:),x), ...
    'alpha_B',0,'alpha_C',1,'last_case_id','');
end

function [tracker,records,raw_files]=rebisect_event(tracker,records,raw_files, ...
        library,cfg,outdir,checkpoint,event_records,template_records,template_states)
while strcmp(tracker.phase,'rebisect')
    step=tracker.rebisection_step;
    if step>cfg.dynamic_edge.maximum_rebisection_steps
        tracker.phase='complete';tracker.status='maximum_rebisection_steps_reached';break;
    end
    alpha=0.5*(tracker.alpha_B+tracker.alpha_C);
    midpoint=interpolate_states(tracker.state_B,tracker.state_C,0.5);
    case_id=sprintf('event%02d_step%02d',tracker.event_index,step);
    raw_file=fullfile(outdir,'raw_cases',[case_id,'.mat']);
    before_B=tracker.alpha_B;before_C=tracker.alpha_C;
    fprintf('Event %d rebisection %d: alpha=%.10f, separation=%.4e.\n', ...
        tracker.event_index,step,alpha, ...
        pair_distance(tracker.state_B,tracker.state_C,cfg.grid.Lx));
    if isfile(raw_file)
        loaded=load(raw_file,'edge_case');edge_case=loaded.edge_case;
        fprintf('  Recovered saved case %s.\n',case_id);
    else
        partial_file=fullfile(outdir,'raw_cases',[case_id,'_partial.mat']);
        edge_case=brusselator_evolve_and_classify_state( ...
            case_id,midpoint,library,cfg,partial_file);
        save(raw_file,'edge_case','-v7.3');
    end
    outcome=edge_case.classification.outcome;update='obstruction';
    if strcmp(outcome,cfg.dynamic_edge.left_outcome)
        tracker.state_B=midpoint;tracker.alpha_B=alpha;update='B_endpoint';
    elseif strcmp(outcome,cfg.dynamic_edge.right_outcome)
        tracker.state_C=midpoint;tracker.alpha_C=alpha;update='C_endpoint';
    else
        tracker.phase='complete';tracker.status=['event_obstructed_by_',outcome];
    end
    separation=pair_distance(tracker.state_B,tracker.state_C,cfg.grid.Lx);
    record=summarize_case(edge_case,tracker.event_index,step,alpha,before_B, ...
        before_C,tracker.alpha_B,tracker.alpha_C,separation,update,raw_file);
    records(end+1,1)=record;raw_files{end+1,1}=raw_file; %#ok<AGROW>
    tracker.rebisection_step=step+1;tracker.last_case_id=case_id;
    if strcmp(tracker.phase,'rebisect')&& ...
            separation<=cfg.dynamic_edge.restart_separation_target
        tracker.phase='evolve_segment';tracker.status='event_rebracketed';
    end
    write_checkpoint(checkpoint,outdir,cfg,tracker,records,raw_files, ...
        event_records,template_records,template_states);
end
end

function [tracker,data,record]=evolve_segment(tracker,cfg,library)
par=brusselator_make_parameters(cfg,@(t)cfg.edge_tracking.frozen_b);
initial_separation=pair_distance(tracker.state_B,tracker.state_C,cfg.grid.Lx);
fprintf('Evolving event-%d renewed pair at global time %.3f.\n', ...
    tracker.event_index,tracker.global_time);
[~,SB,VB,UB]=solve_brusselator_1d_forced(tracker.state_B,par, ...
    cfg.dynamic_edge.segment_maximum_duration,0);
[~,SC,VC,UC]=solve_brusselator_1d_forced(tracker.state_C,par, ...
    cfg.dynamic_edge.segment_maximum_duration,0);
if numel(SB)~=numel(SC)||max(abs(SB-SC))>1e-12
    error('run_dynamic_edge_tracking_BC_multicycle:SegmentTimes', ...
        'B/C segment time grids differ.');
end
[separation,reflected]=pair_distance_history(UB,VB,UC,VC,cfg.grid.Lx);
index=find(separation>=cfg.dynamic_edge.separation_threshold,1);
if isempty(index)
    tracker.phase='complete';tracker.status='segment_did_not_reach_threshold';
    index=numel(SB);threshold=false;
else
    threshold=true;
end
S=SB(1:index);UB=UB(1:index,:);VB=VB(1:index,:);
UC=UC(1:index,:);VC=VC(1:index,:);separation=separation(1:index);
reflected=reflected(1:index);shadow_U=0.5*(UB+UC);shadow_V=0.5*(VB+VC);
[shadow_classification,shadow_history]=brusselator_classify_five_way_trajectory( ...
    S,shadow_U,shadow_V,library,'late_window',min(cfg.edge_tracking.late_window,S(end)), ...
    'sample_interval',cfg.edge_tracking.distance_sample_interval, ...
    'median_threshold',cfg.edge_tracking.median_distance_threshold, ...
    'maximum_threshold',cfg.edge_tracking.maximum_distance_threshold,'Lx',cfg.grid.Lx);
data=struct('segment_index',tracker.event_index,'S',S,'U_B',UB,'V_B',VB, ...
    'U_C',UC,'V_C',VC,'shadow_U',shadow_U,'shadow_V',shadow_V, ...
    'separation',separation,'reflected_separation',reflected, ...
    'shadow_classification_diagnostic',shadow_classification, ...
    'shadow_distance_history',shadow_history, ...
    'warning',['the midpoint shadow is not an exact PDE trajectory; recurrence ', ...
        'and residuals use explicitly identified states as documented']);
record=empty_event_record();record.event_index=tracker.event_index;
record.global_start_time=tracker.global_time;
record.global_end_time=tracker.global_time+S(end);record.segment_duration=S(end);
record.rebisection_steps=tracker.rebisection_step-1;
record.initial_separation=initial_separation;record.final_separation=separation(end);
record.threshold_reached=threshold;record.raw_file='';
end

function [record,state]=extract_template(segment,event_index,cfg,library)
start=cfg.multicycle.phase_sample_delay;
analysis_B=brusselator_periodic_recurrence(segment.S,segment.U_B,segment.V_B, ...
    start,cfg.multicycle.phase_recurrence_threshold);
analysis_C=brusselator_periodic_recurrence(segment.S,segment.U_C,segment.V_C, ...
    start,cfg.multicycle.phase_recurrence_threshold);
analysis_shadow=brusselator_periodic_recurrence(segment.S,segment.shadow_U, ...
    segment.shadow_V,start,cfg.multicycle.phase_recurrence_threshold);
period=0.5*(analysis_B.section.period+analysis_C.section.period);
phase_vector=analysis_shadow.section.states(1,:);N=size(segment.U_B,2);
x=linspace(0,cfg.grid.Lx,N)';state=[x,phase_vector(1:N)',phase_vector(N+1:end)'];
phase_time=analysis_shadow.section.times(1);
[residual,reflected_residual]=flow_residual(state,period,cfg);
[classification,history]=brusselator_classify_five_way_trajectory( ...
    segment.S,segment.U_B,segment.V_B,library, ...
    'late_window',cfg.edge_tracking.late_window, ...
    'sample_interval',cfg.edge_tracking.distance_sample_interval, ...
    'median_threshold',cfg.edge_tracking.median_distance_threshold, ...
    'maximum_threshold',cfg.edge_tracking.maximum_distance_threshold,'Lx',cfg.grid.Lx);
record=empty_template_record();record.event_index=event_index;
record.phase_time=phase_time;record.period_B=analysis_B.section.period;
record.period_C=analysis_C.section.period;record.period_mean=period;
record.period_difference=abs(record.period_B-record.period_C);
record.period_cv_B=analysis_B.section.period_cv;
record.period_cv_C=analysis_C.section.period_cv;
record.crossing_count_B=analysis_B.section.crossing_count;
record.crossing_count_C=analysis_C.section.crossing_count;
record.flow_residual=residual;record.reflected_flow_residual=reflected_residual;
record.segment_B_outcome_diagnostic=classification.outcome;
record.segment_B_minimum_known_distance=min(history.distances,[],'all');
end

function [residual,reflected]=flow_residual(state,period,cfg)
par=brusselator_make_parameters(cfg,@(t)cfg.edge_tracking.frozen_b);
[final_state,~,~,~]=solve_brusselator_1d_forced(state,par,period,0);
residual=pair_distance(state,final_state,cfg.grid.Lx);
reflected=pair_reflected_distance(state,final_state,cfg.grid.Lx);
end

function [distance,components,uses_reflection]=template_distance(A,B,Lx)
direct=pair_distance(A,B,Lx);reflected=pair_reflected_distance(A,B,Lx);
components=[direct,reflected];uses_reflection=reflected<direct;
distance=min(components);
end

function record=empty_case_record()
record=struct('event_index',NaN,'step',NaN,'case_id','','alpha',NaN, ...
    'outcome','','tested_duration',NaN,'decision_time',NaN, ...
    'runtime_seconds',NaN,'alpha_B_before',NaN,'alpha_C_before',NaN, ...
    'alpha_B_after',NaN,'alpha_C_after',NaN,'state_separation_after',NaN, ...
    'update','','late_median_distance_A',NaN,'late_median_distance_B',NaN, ...
    'late_median_distance_C',NaN,'late_median_distance_R5',NaN, ...
    'late_median_distance_R6',NaN,'raw_file','');
end

function record=summarize_case(data,event,step,alpha,bB,bC,aB,aC,separation,update,raw)
record=empty_case_record();record.event_index=event;record.step=step;
record.case_id=data.case_id;record.alpha=alpha;
record.outcome=data.classification.outcome;record.tested_duration=data.duration;
record.decision_time=data.decision_time;record.runtime_seconds=data.runtime_seconds;
record.alpha_B_before=bB;record.alpha_C_before=bC;
record.alpha_B_after=aB;record.alpha_C_after=aC;
record.state_separation_after=separation;record.update=update;
labels={'A','B','C','R5','R6'};
for k=1:numel(labels)
    record.(['late_median_distance_',labels{k}])= ...
        data.classification.late_median_distances(k);
end
record.raw_file=raw;
end

function record=empty_event_record()
record=struct('event_index',NaN,'global_start_time',NaN,'global_end_time',NaN, ...
    'segment_duration',NaN,'rebisection_steps',NaN,'initial_separation',NaN, ...
    'final_separation',NaN,'threshold_reached',false,'raw_file','');
end

function record=empty_template_record()
record=struct('event_index',NaN,'phase_time',NaN,'period_B',NaN,'period_C',NaN, ...
    'period_mean',NaN,'period_difference',NaN,'period_cv_B',NaN,'period_cv_C',NaN, ...
    'crossing_count_B',NaN,'crossing_count_C',NaN,'flow_residual',NaN, ...
    'reflected_flow_residual',NaN,'previous_template_distance',NaN, ...
    'previous_template_direct_distance',NaN, ...
    'previous_template_reflected_distance',NaN, ...
    'previous_template_uses_reflection',false, ...
    'segment_B_outcome_diagnostic','','segment_B_minimum_known_distance',NaN);
end

function summary=make_summary(parent,tracker,events,templates,cfg)
later=templates(max(1,end-2):end);distances=[later.previous_template_distance];
summary=struct('status',tracker.status,'parent_shadow_time', ...
    parent.summary.total_shadow_time,'new_event_count',numel(events), ...
    'total_shadow_time',tracker.global_time,'final_period',templates(end).period_mean, ...
    'last_period_change',abs(templates(end).period_mean-templates(end-1).period_mean), ...
    'last_template_distance',templates(end).previous_template_distance, ...
    'last_flow_residual',templates(end).flow_residual, ...
    'maximum_last_three_template_distance',max(distances), ...
    'newton_seed_gate',cfg.multicycle.newton_seed_gate);
end

function write_checkpoint(file,outdir,cfg,tracker,case_records,raw_files, ...
        event_records,template_records,template_states)
save(file,'cfg','tracker','case_records','raw_files','event_records', ...
    'template_records','template_states','-v7.3');
if ~isempty(case_records)
    writetable(struct2table(case_records,'AsArray',true), ...
        fullfile(outdir,'progress_rebisection.csv'));
end
if ~isempty(template_records)
    writetable(struct2table(template_records,'AsArray',true), ...
        fullfile(outdir,'progress_template_convergence.csv'));
end
end

function make_figures(outdir,parent,events,templates,cfg)
cycles=[templates.event_index];periods=[templates.period_mean];
fig=figure('Color','w','Position',[100 100 1000 750]);layout=tiledlayout(fig,3,1);
ax=nexttile(layout);plot(ax,cycles,periods,'o-','LineWidth',1.4,'MarkerFaceColor',[0.2 0.5 0.8]);
ylabel(ax,'period');title(ax,'Phase-fixed edge-template convergence');style_axes(ax);
ax=nexttile(layout);semilogy(ax,cycles(2:end),[templates(2:end).previous_template_distance], ...
    'o-','LineWidth',1.4,'MarkerFaceColor',[0.85 0.45 0.15]);
ylabel(ax,'template distance');style_axes(ax);
ax=nexttile(layout);semilogy(ax,cycles,[templates.flow_residual],'o-', ...
    'LineWidth',1.4,'MarkerFaceColor',[0.25 0.65 0.35]);
xlabel(ax,'post-rebracketing cycle');ylabel(ax,'one-period residual');style_axes(ax);
exportgraphics(fig,fullfile(outdir,'phase_template_convergence.png'),'Resolution',300);
close(fig);

fig=figure('Color','w','Position',[100 100 900 500]);ax=axes(fig);hold(ax,'on');
times=[parent.summary.initial_divergence_time, ...
    parent.summary.post_rebracketing_divergence_time,[events.segment_duration]];
bar(ax,0:numel(times)-1,times,'FaceColor',[0.2 0.5 0.8]);
xlabel(ax,'edge-shadow segment');ylabel(ax,'time to separation threshold');
title(ax,sprintf('Dynamic edge segments, threshold %.3g', ...
    cfg.dynamic_edge.separation_threshold));style_axes(ax);
exportgraphics(fig,fullfile(outdir,'multicycle_shadow_durations.png'),'Resolution',300);
close(fig);
end

function state=interpolate_states(B,C,alpha)
state=B;state(:,2:3)=(1-alpha)*B(:,2:3)+alpha*C(:,2:3);
end

function state=state_from_row(u,v,x)
state=[x,u(:),v(:)];
end

function [distance,reflected]=pair_distance_history(UB,VB,UC,VC,Lx)
N=size(UB,2);dx=Lx/(N-1);w=dx*ones(1,N);w([1,end])=0.5*dx;
numerator=sum(((UB-UC).^2+(VB-VC).^2).*w,2);
energy_B=sum((UB.^2+VB.^2).*w,2);energy_C=sum((UC.^2+VC.^2).*w,2);
denominator=max(0.5*(energy_B+energy_C),eps);
distance=sqrt(numerator./denominator);
ref_num=sum(((UB-fliplr(UC)).^2+(VB-fliplr(VC)).^2).*w,2);
reflected=sqrt(ref_num./denominator);
end

function distance=pair_distance(B,C,Lx)
[distance,~]=pair_distance_history(B(:,2)',B(:,3)',C(:,2)',C(:,3)',Lx);
end

function distance=pair_reflected_distance(B,C,Lx)
[~,distance]=pair_distance_history(B(:,2)',B(:,3)',C(:,2)',C(:,3)',Lx);
end

function style_axes(ax)
set(ax,'Color','w','XColor','k','YColor','k','GridColor',[0.75 0.75 0.75]);
grid(ax,'on');box(ax,'on');ax.Title.Color='k';ax.XLabel.Color='k';ax.YLabel.Color='k';
end
