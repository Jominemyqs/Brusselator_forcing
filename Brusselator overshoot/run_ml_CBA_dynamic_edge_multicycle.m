function results = run_ml_CBA_dynamic_edge_multicycle(bracket_id)
%RUN_ML_CBA_DYNAMIC_EDGE_MULTICYCLE Continue three rebracketing events.

cfg=brusselator_ml_CBA_multicycle_config(bracket_id);
outdir=fullfile(cfg.output.root,cfg.experiment_name);
checkpoint=fullfile(outdir,'progress_multicycle_pair.mat');
final_file=fullfile(outdir,'dynamic_pair_multicycle.mat');
if isfile(final_file)
    error('run_ml_CBA_dynamic_edge_multicycle:OutputExists', ...
        'Refusing to overwrite completed multicycle result: %s',final_file);
end
[parent,parent_segment,library]=load_parent(cfg);validate_parent(parent,parent_segment);
if isfile(checkpoint)
    saved=load(checkpoint,'tracker','case_records','raw_files','event_records');
    tracker=saved.tracker;case_records=saved.case_records;
    raw_files=saved.raw_files;event_records=saved.event_records;
    fprintf('Resuming %s at event %d phase %s.\n', ...
        cfg.experiment_name,tracker.event_index,tracker.phase);
else
    if isfolder(outdir)
        error('run_ml_CBA_dynamic_edge_multicycle:UnsafeResume', ...
            'Existing output has no recognized checkpoint: %s',outdir);
    end
    mkdir(outdir);mkdir(fullfile(outdir,'raw_cases'));mkdir(fullfile(outdir,'segments'));
    case_records=empty_case_record([]);raw_files={};event_records=empty_event_record([]);
    tracker=struct('phase','rebisect','status','parent_divergence_ready', ...
        'event_index',cfg.multicycle_pair.first_new_event_index, ...
        'rebisection_step',1,'global_time',parent.tracker.global_time, ...
        'state_left',parent.tracker.state_left,'state_right',parent.tracker.state_right, ...
        'alpha_left',0,'alpha_right',1,'restart_separation',NaN);
    brusselator_write_metadata(outdir,cfg, ...
        brusselator_run_metadata(cfg,mfilename,NaN));
    write_checkpoint();
end
timer=tic;last_event=cfg.multicycle_pair.first_new_event_index+ ...
    cfg.multicycle_pair.additional_events-1;
while tracker.event_index<=last_event && ~strcmp(tracker.phase,'complete')
    while strcmp(tracker.phase,'rebisect')
        step=tracker.rebisection_step;
        if step>cfg.dynamic_pair.maximum_rebisection_steps
            tracker.phase='complete';tracker.status='maximum_rebisection_steps_reached';
            break;
        end
        before_left=tracker.alpha_left;before_right=tracker.alpha_right;
        alpha=0.5*(before_left+before_right);
        midpoint=interpolate_states(tracker.state_left,tracker.state_right,0.5);
        case_id=sprintf('event%02d_step%02d',tracker.event_index,step);
        raw_file=fullfile(outdir,'raw_cases',[case_id,'.mat']);
        partial_file=fullfile(outdir,'raw_cases',[case_id,'_partial.mat']);
        fprintf('%s event %d rebisection %d: separation %.4e.\n', ...
            cfg.dynamic_pair.bracket_id,tracker.event_index,step, ...
            pair_distance(tracker.state_left,tracker.state_right,cfg.grid.Lx));
        if isfile(raw_file)
            loaded=load(raw_file,'edge_case');edge_case=loaded.edge_case;
        else
            edge_case=brusselator_evolve_and_classify_state( ...
                case_id,midpoint,library,cfg,partial_file);
            save(raw_file,'edge_case','-v7.3');
            if isfile(partial_file),delete(partial_file);end
        end
        outcome=edge_case.classification.outcome;update='obstruction';
        if strcmp(outcome,cfg.dynamic_pair.left_outcome)
            tracker.state_left=midpoint;tracker.alpha_left=alpha;update='left_endpoint';
        elseif strcmp(outcome,cfg.dynamic_pair.right_outcome)
            tracker.state_right=midpoint;tracker.alpha_right=alpha;update='right_endpoint';
        else
            tracker.phase='complete';tracker.status=['event_obstructed_by_',outcome];
        end
        separation=pair_distance(tracker.state_left,tracker.state_right,cfg.grid.Lx);
        record=summarize_case(edge_case,tracker.event_index,step,alpha, ...
            before_left,before_right,tracker.alpha_left,tracker.alpha_right, ...
            separation,update,raw_file);
        case_records(end+1,1)=record;raw_files{end+1,1}=raw_file; %#ok<AGROW>
        tracker.rebisection_step=step+1;
        if strcmp(tracker.phase,'rebisect')&& ...
                separation<=cfg.dynamic_pair.restart_separation_target
            tracker.restart_separation=separation;tracker.phase='evolve_segment';
            tracker.status='event_rebracketed';
        end
        write_checkpoint();
    end
    if strcmp(tracker.phase,'evolve_segment')
        [segment_data,event_record,left_final,right_final]=evolve_segment( ...
            tracker.state_left,tracker.state_right,tracker.event_index, ...
            tracker.global_time,library,cfg);
        segment_file=fullfile(outdir,'segments',sprintf('segment_%02d.mat', ...
            tracker.event_index));save(segment_file,'segment_data','-v7.3');
        event_record.raw_file=segment_file;
        [minimum_exact,maximum_shadow]=exact_diagnostic(segment_data,cfg);
        event_record.minimum_exact_edge_distance=minimum_exact;
        event_record.maximum_exact_edge_shadow_duration=maximum_shadow;
        event_records(end+1,1)=event_record; %#ok<AGROW>
        tracker.global_time=event_record.global_end_time;
        tracker.state_left=left_final;tracker.state_right=right_final;
        if ~event_record.threshold_reached
            tracker.phase='complete';tracker.status='segment_did_not_reach_threshold';
        elseif tracker.event_index>=last_event
            tracker.phase='complete';tracker.status='requested_multicycle_run_complete';
        else
            tracker.event_index=tracker.event_index+1;tracker.rebisection_step=1;
            tracker.alpha_left=0;tracker.alpha_right=1;tracker.phase='rebisect';
            tracker.status='next_divergence_ready';
        end
        write_checkpoint();
    end
end
summary=make_summary(parent,tracker,case_records,event_records,cfg);
writetable(struct2table(case_records,'AsArray',true), ...
    fullfile(outdir,'multicycle_rebisection_cases.csv'));
writetable(struct2table(event_records,'AsArray',true), ...
    fullfile(outdir,'multicycle_events.csv'));
writetable(struct2table(summary,'AsArray',true), ...
    fullfile(outdir,'multicycle_summary.csv'));
make_figures(outdir,parent,event_records,cfg);
results=struct('configuration',cfg,'parent_result_file', ...
    cfg.multicycle_pair.parent_result_file,'tracker',tracker, ...
    'case_records',case_records,'raw_files',{raw_files}, ...
    'event_records',event_records,'summary',summary);
save(final_file,'results','-v7.3');
brusselator_write_metadata(outdir,cfg, ...
    brusselator_run_metadata(cfg,mfilename,toc(timer)));
disp(struct2table(event_records,'AsArray',true));disp(struct2table(summary,'AsArray',true));
fprintf('Multicycle %s edge tracking saved in: %s\n', ...
    cfg.dynamic_pair.bracket_id,outdir);

    function write_checkpoint()
        save(checkpoint,'cfg','tracker','case_records','raw_files', ...
            'event_records','-v7.3');
        if ~isempty(case_records)
            writetable(struct2table(case_records,'AsArray',true), ...
                fullfile(outdir,'progress_rebisection.csv'));
        end
    end
end

function [parent,segment,library]=load_parent(cfg)
if ~isfile(cfg.multicycle_pair.parent_result_file)|| ...
        ~isfile(cfg.multicycle_pair.parent_post_segment_file)
    error('A completed one-event parent is missing.');
end
loaded=load(cfg.multicycle_pair.parent_result_file,'results');parent=loaded.results;
loaded=load(cfg.multicycle_pair.parent_post_segment_file,'segment_data');
segment=loaded.segment_data;[library,~]=brusselator_build_ml_outcome_library(cfg);
end

function validate_parent(parent,segment)
if ~strcmp(parent.summary.status,'one_rebracketing_cycle_complete')|| ...
        segment.segment_index~=1||~parent.segments(end).threshold_reached
    error('Configured parent is not a completed divergent one-event result.');
end
end

function [data,record,left_final,right_final]=evolve_segment( ...
        left_initial,right_initial,index,global_start,library,cfg)
par=brusselator_make_parameters(cfg,@(t)cfg.edge_tracking.frozen_b);
fprintf('Evolving %s event %d from global time %.3f.\n', ...
    cfg.dynamic_pair.bracket_id,index,global_start);
[~,SL,VL,UL]=solve_brusselator_1d_forced(left_initial,par, ...
    cfg.dynamic_pair.segment_maximum_duration,0);
[~,SR,VR,UR]=solve_brusselator_1d_forced(right_initial,par, ...
    cfg.dynamic_pair.segment_maximum_duration,0);
if numel(SL)~=numel(SR)||max(abs(SL-SR))>1e-12,error('Segment time grids differ.');end
[separation,reflected]=pair_distance_history(UL,VL,UR,VR,cfg.grid.Lx);
stop=find(separation>=cfg.dynamic_pair.separation_threshold,1);
threshold=~isempty(stop);if ~threshold,stop=numel(SL);end
S=SL(1:stop);UL=UL(1:stop,:);VL=VL(1:stop,:);
UR=UR(1:stop,:);VR=VR(1:stop,:);separation=separation(1:stop);
reflected=reflected(1:stop);shadow_U=0.5*(UL+UR);shadow_V=0.5*(VL+VR);
[classification,history]=brusselator_classify_five_way_trajectory( ...
    S,shadow_U,shadow_V,library,'late_window',min(cfg.edge_tracking.late_window,S(end)), ...
    'sample_interval',cfg.edge_tracking.distance_sample_interval, ...
    'median_threshold',cfg.edge_tracking.median_distance_threshold, ...
    'maximum_threshold',cfg.edge_tracking.maximum_distance_threshold,'Lx',cfg.grid.Lx);
data=struct('segment_index',index,'S',S,'U_left',UL,'V_left',VL, ...
    'U_right',UR,'V_right',VR,'shadow_U',shadow_U,'shadow_V',shadow_V, ...
    'separation',separation,'reflected_separation',reflected, ...
    'shadow_classification_diagnostic',classification, ...
    'shadow_distance_history',history, ...
    'warning','midpoint shadow is not an exact PDE trajectory');
record=empty_event_record();record.event_index=index;
record.global_start_time=global_start;record.global_end_time=global_start+S(end);
record.segment_duration=S(end);record.initial_separation=separation(1);
record.final_separation=separation(end);record.threshold_reached=threshold;
x=left_initial(:,1);left_final=[x,UL(end,:)',VL(end,:)'];
right_final=[x,UR(end,:)',VR(end,:)'];
end

function [minimum,shadow]=exact_diagnostic(segment,cfg)
if ~cfg.dynamic_pair.exact_edge_comparison,minimum=NaN;shadow=NaN;return;end
loaded=load(cfg.dynamic_pair.exact_edge_file,'results');edge=loaded.results;
minimum=inf;shadow=0;
for side={'left','right'}
    times=(segment.S(1):0.2:segment.S(end))';
    U=interp1(segment.S,segment.(['U_',side{1}]),times,'pchip');
    V=interp1(segment.S,segment.(['V_',side{1}]),times,'pchip');
    distance=zeros(size(times));
    for k=1:numel(times)
        value=brusselator_physical_orbit_distance(U(k,:),V(k,:), ...
            edge.orbit_U,edge.orbit_V,cfg.grid.Lx);distance(k)=value.distance;
    end
    minimum=min(minimum,min(distance));mask=distance< ...
        cfg.dynamic_pair.exact_edge_neighborhood_threshold;
    change=diff([false;mask;false]);starts=find(change==1);ends=find(change==-1)-1;
    if ~isempty(starts),shadow=max(shadow,max(times(ends)-times(starts)));end
end
end

function summary=make_summary(parent,tracker,cases,events,cfg)
parent_exact=parent.summary.minimum_exact_edge_distance;
values=[events.minimum_exact_edge_distance];finite=isfinite(values);
if any(finite),minimum_exact=min([parent_exact,values(finite)]);
else,minimum_exact=parent_exact;end
durations=[parent.summary.post_rebracketing_shadow_duration,[events.segment_duration]];
summary=struct('bracket_id',cfg.dynamic_pair.bracket_id,'status',tracker.status, ...
    'additional_event_count',numel(events),'additional_rebisection_cases',numel(cases), ...
    'assembled_global_time',tracker.global_time, ...
    'segment_durations',durations,'minimum_exact_edge_distance',minimum_exact, ...
    'known_exact_edge_recovered',isfinite(minimum_exact)&&minimum_exact< ...
        cfg.dynamic_pair.exact_edge_neighborhood_threshold, ...
    'interpretation',cfg.multicycle_pair.interpretation);
end

function record=empty_case_record(varargin)
record=struct('event_index',NaN,'step',NaN,'alpha',NaN,'outcome','', ...
    'tested_duration',NaN,'decision_time',NaN,'runtime_seconds',NaN, ...
    'alpha_left_before',NaN,'alpha_right_before',NaN, ...
    'alpha_left_after',NaN,'alpha_right_after',NaN, ...
    'state_separation_after',NaN,'update','','raw_file','');
if nargin>0,record=record([]);end
end

function record=summarize_case(data,event,step,alpha,bL,bR,aL,aR,separation,update,raw)
record=empty_case_record();record.event_index=event;record.step=step;
record.alpha=alpha;record.outcome=data.classification.outcome;
record.tested_duration=data.duration;record.decision_time=data.decision_time;
record.runtime_seconds=data.runtime_seconds;record.alpha_left_before=bL;
record.alpha_right_before=bR;record.alpha_left_after=aL;record.alpha_right_after=aR;
record.state_separation_after=separation;record.update=update;record.raw_file=raw;
end

function record=empty_event_record(varargin)
record=struct('event_index',NaN,'global_start_time',NaN,'global_end_time',NaN, ...
    'segment_duration',NaN,'initial_separation',NaN,'final_separation',NaN, ...
    'threshold_reached',false,'minimum_exact_edge_distance',NaN, ...
    'maximum_exact_edge_shadow_duration',NaN,'raw_file','');
if nargin>0,record=record([]);end
end

function state=interpolate_states(left,right,alpha)
state=left;state(:,2:3)=(1-alpha)*left(:,2:3)+alpha*right(:,2:3);
end

function [distance,reflected]=pair_distance_history(UL,VL,UR,VR,Lx)
N=size(UL,2);dx=Lx/(N-1);w=dx*ones(1,N);w([1,end])=0.5*dx;
energy_left=sum((UL.^2+VL.^2).*w,2);energy_right=sum((UR.^2+VR.^2).*w,2);
denominator=max(0.5*(energy_left+energy_right),eps);
distance=sqrt(sum(((UL-UR).^2+(VL-VR).^2).*w,2)./denominator);
reflected=sqrt(sum(((UL-fliplr(UR)).^2+(VL-fliplr(VR)).^2).*w,2)./denominator);
end

function distance=pair_distance(left,right,Lx)
[distance,~]=pair_distance_history(left(:,2)',left(:,3)',right(:,2)',right(:,3)',Lx);
end

function make_figures(outdir,parent,events,cfg)
durations=[parent.summary.post_rebracketing_shadow_duration,[events.segment_duration]];
fig=figure('Visible','off','Color','w');bar(1:numel(durations),durations);
xlabel('dynamic segment');ylabel('time to separation 0.02');grid on;
title(sprintf('%s repeated edge-shadow durations',cfg.dynamic_pair.bracket_id));
exportgraphics(fig,fullfile(outdir,'multicycle_shadow_durations.png'),'Resolution',300);close(fig);
if cfg.dynamic_pair.exact_edge_comparison
    values=[parent.summary.minimum_exact_edge_distance,[events.minimum_exact_edge_distance]];
    fig=figure('Visible','off','Color','w');semilogy(1:numel(values),values,'o-','LineWidth',1.4);
    yline(cfg.dynamic_pair.exact_edge_neighborhood_threshold,'--k','E_{BC} neighborhood');
    xlabel('dynamic segment');ylabel('minimum exact-edge distance');grid on;
    title('Tracked C--B boundary versus known E_{BC}');
    exportgraphics(fig,fullfile(outdir,'exact_edge_distance_by_cycle.png'),'Resolution',300);close(fig);
end
end
