function results = run_ml_CBA_dynamic_edge_tracking(bracket_id,cfg)
%RUN_ML_CBA_DYNAMIC_EDGE_TRACKING One-event five-way dynamic edge tracking.
%   An optional configuration permits the same audited tracking machinery
%   to operate on independently generated, schema-compatible brackets. The
%   original one-argument held-out CBA calls are unchanged.

if nargin<2 || isempty(cfg)
    cfg=brusselator_ml_CBA_edge_tracking_config(bracket_id);
end
outdir=fullfile(cfg.output.root,cfg.experiment_name);
checkpoint=fullfile(outdir,'progress_dynamic_pair.mat');
final_file=fullfile(outdir,'dynamic_pair_edge_tracking.mat');
if isfile(final_file)
    error('run_ml_CBA_dynamic_edge_tracking:OutputExists', ...
        'Refusing to overwrite completed tracking result: %s',final_file);
end
if ~isfile(cfg.dynamic_pair.bracket_file)
    error('run_ml_CBA_dynamic_edge_tracking:MissingBracket', ...
        'The discovered-bracket table is missing.');
end
[left_source,right_source,bracket,library]=load_sources(cfg);
validate_sources(left_source,right_source,bracket,cfg);
if isfile(checkpoint)
    saved=load(checkpoint,'tracker','case_records','raw_files','segments');
    tracker=saved.tracker;case_records=saved.case_records;
    raw_files=saved.raw_files;segments=saved.segments;
    fprintf('Resuming %s in phase %s.\n',cfg.experiment_name,tracker.phase);
else
    if isfolder(outdir)
        error('run_ml_CBA_dynamic_edge_tracking:UnsafeResume', ...
            'Existing output has no recognized checkpoint: %s',outdir);
    end
    mkdir(outdir);mkdir(fullfile(outdir,'raw_cases'));mkdir(fullfile(outdir,'segments'));
    case_records=empty_case_record([]);raw_files={};segments=empty_segment_record([]);
    tracker=struct('phase','source_segment','status','initialized', ...
        'global_time',0,'event_index',0,'rebisection_step',1, ...
        'state_left',left_source.initial_state, ...
        'state_right',right_source.initial_state, ...
        'alpha_left',0,'alpha_right',1,'restart_separation',NaN);
    brusselator_write_metadata(outdir,cfg, ...
        brusselator_run_metadata(cfg,mfilename,NaN));
    write_checkpoint();
end
timer=tic;

if strcmp(tracker.phase,'source_segment')
    [segment_data,segment_record,left_state,right_state]=evolve_segment( ...
        tracker.state_left,tracker.state_right,0,tracker.global_time,library,cfg);
    segment_file=fullfile(outdir,'segments','segment_00_source.mat');
    save(segment_file,'segment_data','-v7.3');segment_record.raw_file=segment_file;
    segments(end+1,1)=segment_record;tracker.global_time=segment_record.global_end_time;
    tracker.state_left=left_state;tracker.state_right=right_state;
    if segment_record.threshold_reached
        tracker.phase='rebisect';tracker.status='source_divergence_ready';
    else
        tracker.phase='complete';tracker.status='source_segment_did_not_diverge';
    end
    write_checkpoint();
end

while strcmp(tracker.phase,'rebisect')
    step=tracker.rebisection_step;
    if step>cfg.dynamic_pair.maximum_rebisection_steps
        tracker.phase='complete';tracker.status='maximum_rebisection_steps_reached';
        break;
    end
    before_left=tracker.alpha_left;before_right=tracker.alpha_right;
    alpha=0.5*(before_left+before_right);
    midpoint=interpolate_states(tracker.state_left,tracker.state_right,0.5);
    case_id=sprintf('event01_step%02d',step);
    raw_file=fullfile(outdir,'raw_cases',[case_id,'.mat']);
    partial_file=fullfile(outdir,'raw_cases',[case_id,'_partial.mat']);
    fprintf('%s rebisection %d: alpha=%.10f, separation=%.4e.\n', ...
        cfg.dynamic_pair.bracket_id,step,alpha, ...
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
        tracker.phase='complete';tracker.status=['rebisection_obstructed_by_',outcome];
    end
    separation=pair_distance(tracker.state_left,tracker.state_right,cfg.grid.Lx);
    record=summarize_case(edge_case,step,alpha,before_left,before_right, ...
        tracker.alpha_left,tracker.alpha_right,separation,update,raw_file);
    case_records(end+1,1)=record;raw_files{end+1,1}=raw_file; %#ok<AGROW>
    tracker.rebisection_step=step+1;
    if strcmp(tracker.phase,'rebisect')&& ...
            separation<=cfg.dynamic_pair.restart_separation_target
        tracker.restart_separation=separation;
        tracker.phase='post_segment';tracker.status='dynamic_rebracketing_complete';
    end
    write_checkpoint();
end

if strcmp(tracker.phase,'post_segment')
    [segment_data,segment_record,left_state,right_state]=evolve_segment( ...
        tracker.state_left,tracker.state_right,1,tracker.global_time,library,cfg);
    segment_file=fullfile(outdir,'segments','segment_01_post_rebracketing.mat');
    save(segment_file,'segment_data','-v7.3');segment_record.raw_file=segment_file;
    segments(end+1,1)=segment_record;tracker.global_time=segment_record.global_end_time;
    tracker.state_left=left_state;tracker.state_right=right_state;
    tracker.phase='complete';
    if segment_record.threshold_reached
        tracker.status='one_rebracketing_cycle_complete';
    else
        tracker.status='post_segment_did_not_reach_threshold';
    end
    write_checkpoint();
end

diagnostics=make_diagnostics(segments,cfg);
summary=make_summary(tracker,case_records,segments,diagnostics,cfg);
writetable(struct2table(case_records,'AsArray',true), ...
    fullfile(outdir,'dynamic_rebisection_cases.csv'));
writetable(struct2table(segments,'AsArray',true), ...
    fullfile(outdir,'dynamic_edge_segments.csv'));
writetable(struct2table(diagnostics,'AsArray',true), ...
    fullfile(outdir,'edge_diagnostics.csv'));
writetable(struct2table(summary,'AsArray',true), ...
    fullfile(outdir,'dynamic_edge_summary.csv'));
make_figures(outdir,segments,case_records,cfg);
results=struct('configuration',cfg,'bracket',bracket,'tracker',tracker, ...
    'case_records',case_records,'raw_files',{raw_files},'segments',segments, ...
    'diagnostics',diagnostics,'summary',summary);
save(final_file,'results','-v7.3');
brusselator_write_metadata(outdir,cfg, ...
    brusselator_run_metadata(cfg,mfilename,toc(timer)));
disp(struct2table(case_records,'AsArray',true));
disp(struct2table(diagnostics,'AsArray',true));disp(struct2table(summary,'AsArray',true));
fprintf('Dynamic %s edge tracking saved in: %s\n', ...
    cfg.dynamic_pair.bracket_id,outdir);

    function write_checkpoint()
        save(checkpoint,'cfg','tracker','case_records','raw_files','segments','-v7.3');
        if ~isempty(case_records)
            writetable(struct2table(case_records,'AsArray',true), ...
                fullfile(outdir,'progress_rebisection.csv'));
        end
    end
end

function [left,right,bracket,library]=load_sources(cfg)
if isfield(cfg.dynamic_pair,'left_source_file')
    left_file=cfg.dynamic_pair.left_source_file;
    right_file=cfg.dynamic_pair.right_source_file;
    left=load_case(left_file);right=load_case(right_file);
    if ~isfield(left,'lambda') && isfield(cfg.dynamic_pair,'left_coordinate')
        left.lambda=cfg.dynamic_pair.left_coordinate;
        right.lambda=cfg.dynamic_pair.right_coordinate;
    end
    bracket=table(string(cfg.dynamic_pair.bracket_id),left.lambda, ...
        string(cfg.dynamic_pair.left_outcome),right.lambda, ...
        string(cfg.dynamic_pair.right_outcome),string(left_file), ...
        string(right_file),'VariableNames',{'bracket_id','left_lambda', ...
        'left_outcome','right_lambda','right_outcome','left_source_file', ...
        'right_source_file'});
else
    table_data=readtable(cfg.dynamic_pair.bracket_file,'TextType','string', ...
        'Delimiter',',','VariableNamingRule','preserve');
    index=find(table_data.bracket_id==string(cfg.dynamic_pair.bracket_id),1);
    if isempty(index),error('Configured bracket id is absent.');end
    bracket=table_data(index,:);
    left=load_case(char(bracket.left_source_file));
    right=load_case(char(bracket.right_source_file));
end
[library,~]=brusselator_build_ml_outcome_library(cfg);
end

function data=load_case(filename)
loaded=load(filename);
if isfield(loaded,'case_data')
    data=loaded.case_data;
elseif isfield(loaded,'edge_case')
    data=loaded.edge_case;
elseif isfield(loaded,'case_result')
    source=loaded.case_result;
    data=source.frozen;
    data.initial_state=source.landing_state;
    data.lambda=source.protocol.Tup;
else
    error('Source file %s has no supported case structure.',filename);
end
end

function validate_sources(left,right,bracket,cfg)
if ~strcmp(left.classification.outcome,cfg.dynamic_pair.left_outcome)|| ...
        ~strcmp(right.classification.outcome,cfg.dynamic_pair.right_outcome)
    error('Source outcomes do not match the configured ordered bracket.');
end
left_column='left_lambda';right_column='right_lambda';
if isfield(cfg.dynamic_pair,'left_coordinate_column')
    left_column=cfg.dynamic_pair.left_coordinate_column;
    right_column=cfg.dynamic_pair.right_coordinate_column;
end
if abs(left.lambda-bracket.(left_column))>1e-12|| ...
        abs(right.lambda-bracket.(right_column))>1e-12
    error('Source coordinates do not match the bracket table.');
end
direct=pair_distance(left.initial_state,right.initial_state,cfg.grid.Lx);
reflected=pair_reflected_distance(left.initial_state,right.initial_state,cfg.grid.Lx);
if direct>=reflected,error('Source pair is not in its maintained common orientation.');end
end

function [data,record,left_final,right_final]=evolve_segment( ...
        left_initial,right_initial,index,global_start,library,cfg)
par=brusselator_make_parameters(cfg,@(t)cfg.edge_tracking.frozen_b);
fprintf('Evolving %s segment %d from global edge time %.3f.\n', ...
    cfg.dynamic_pair.bracket_id,index,global_start);
initial_separation=pair_distance(left_initial,right_initial,cfg.grid.Lx);
if initial_separation>=cfg.dynamic_pair.separation_threshold
    % Adjacent points in a forcing parameter can already be far apart in
    % state space after a long, sensitive nonautonomous evolution. Record
    % this faithfully and proceed directly to guarded frozen rebracketing.
    S=0;UL=left_initial(:,2)';VL=left_initial(:,3)';
    UR=right_initial(:,2)';VR=right_initial(:,3)';
    separation=initial_separation;reflected=pair_reflected_distance( ...
        left_initial,right_initial,cfg.grid.Lx);
    shadow_U=0.5*(UL+UR);shadow_V=0.5*(VL+VR);
    shadow_classification=struct('outcome','not_evaluated_zero_duration');
    shadow_history=struct('sample_times',[],'distances',[]);
    data=struct('segment_index',index,'S',S,'U_left',UL,'V_left',VL, ...
        'U_right',UR,'V_right',VR,'shadow_U',shadow_U,'shadow_V',shadow_V, ...
        'separation',separation,'reflected_separation',reflected, ...
        'shadow_classification_diagnostic',shadow_classification, ...
        'shadow_distance_history',shadow_history, ...
        'warning',['source states already exceed the separation threshold; ', ...
            'zero pre-rebracketing shadow duration is recorded']);
    record=empty_segment_record();record.segment_index=index;
    record.global_start_time=global_start;record.global_end_time=global_start;
    record.duration=0;record.initial_separation=initial_separation;
    record.final_separation=initial_separation;
    record.maximum_separation=initial_separation;record.threshold_reached=true;
    left_final=left_initial;right_final=right_initial;
    return;
end
[~,SL,VL,UL]=solve_brusselator_1d_forced(left_initial,par, ...
    cfg.dynamic_pair.segment_maximum_duration,0);
[~,SR,VR,UR]=solve_brusselator_1d_forced(right_initial,par, ...
    cfg.dynamic_pair.segment_maximum_duration,0);
if numel(SL)~=numel(SR)||max(abs(SL-SR))>1e-12
    error('Paired segment time grids differ.');
end
[separation,reflected]=pair_distance_history(UL,VL,UR,VR,cfg.grid.Lx);
stop=find(separation>=cfg.dynamic_pair.separation_threshold,1);
threshold=~isempty(stop);if ~threshold,stop=numel(SL);end
S=SL(1:stop);UL=UL(1:stop,:);VL=VL(1:stop,:);
UR=UR(1:stop,:);VR=VR(1:stop,:);separation=separation(1:stop);
reflected=reflected(1:stop);shadow_U=0.5*(UL+UR);shadow_V=0.5*(VL+VR);
[shadow_classification,shadow_history]=brusselator_classify_five_way_trajectory( ...
    S,shadow_U,shadow_V,library,'late_window',min(cfg.edge_tracking.late_window,S(end)), ...
    'sample_interval',cfg.edge_tracking.distance_sample_interval, ...
    'median_threshold',cfg.edge_tracking.median_distance_threshold, ...
    'maximum_threshold',cfg.edge_tracking.maximum_distance_threshold, ...
    'Lx',cfg.grid.Lx);
data=struct('segment_index',index,'S',S,'U_left',UL,'V_left',VL, ...
    'U_right',UR,'V_right',VR,'shadow_U',shadow_U,'shadow_V',shadow_V, ...
    'separation',separation,'reflected_separation',reflected, ...
    'shadow_classification_diagnostic',shadow_classification, ...
    'shadow_distance_history',shadow_history, ...
    'warning','midpoint shadow is not an exact PDE trajectory');
record=empty_segment_record();record.segment_index=index;
record.global_start_time=global_start;record.global_end_time=global_start+S(end);
record.duration=S(end);record.initial_separation=separation(1);
record.final_separation=separation(end);record.maximum_separation=max(separation);
record.threshold_reached=threshold;
x=left_initial(:,1);left_final=[x,UL(end,:)',VL(end,:)'];
right_final=[x,UR(end,:)',VR(end,:)'];
end

function diagnostics=make_diagnostics(segments,cfg)
template=struct('segment_index',NaN,'side','','recurrence_available',false, ...
    'period',NaN,'period_cv',NaN,'crossing_count',NaN, ...
    'minimum_exact_edge_distance',NaN,'exact_edge_shadow_duration',NaN, ...
    'period_relative_error_to_exact_edge',NaN,'error','');
diagnostics=template([]);
for s=1:numel(segments)
    loaded=load(segments(s).raw_file,'segment_data');segment=loaded.segment_data;
    for side={'left','right'}
        item=template;item.segment_index=segments(s).segment_index;item.side=side{1};
        try
            if segment.S(end)-cfg.dynamic_pair.recurrence_sample_delay>10
                analysis=brusselator_periodic_recurrence(segment.S, ...
                    segment.(['U_',side{1}]),segment.(['V_',side{1}]), ...
                    cfg.dynamic_pair.recurrence_sample_delay, ...
                    cfg.dynamic_pair.recurrence_threshold);
                item.recurrence_available=true;item.period=analysis.section.period;
                item.period_cv=analysis.section.period_cv;
                item.crossing_count=analysis.section.crossing_count;
            end
        catch exception
            item.error=exception.message;
        end
        if cfg.dynamic_pair.exact_edge_comparison
            exact=load(cfg.dynamic_pair.exact_edge_file,'results');edge=exact.results;
            [minimum,shadow]=exact_edge_distance(segment.S, ...
                segment.(['U_',side{1}]),segment.(['V_',side{1}]),edge,cfg);
            item.minimum_exact_edge_distance=minimum;
            item.exact_edge_shadow_duration=shadow;
            if item.recurrence_available
                item.period_relative_error_to_exact_edge= ...
                    abs(item.period-edge.period)/edge.period;
            end
        end
        diagnostics(end+1,1)=item; %#ok<AGROW>
    end
end
end

function [minimum,shadow]=exact_edge_distance(S,U,V,edge,cfg)
times=(S(1):0.5:S(end))';sample_U=interp1(S,U,times,'pchip');
sample_V=interp1(S,V,times,'pchip');distance=zeros(size(times));
for k=1:numel(times)
    value=brusselator_physical_orbit_distance(sample_U(k,:),sample_V(k,:), ...
        edge.orbit_U,edge.orbit_V,cfg.grid.Lx);
    distance(k)=value.distance;
end
minimum=min(distance);mask=distance<cfg.dynamic_pair.exact_edge_neighborhood_threshold;
change=diff([false;mask;false]);starts=find(change==1);ends=find(change==-1)-1;
if isempty(starts),shadow=0;else,shadow=max(times(ends)-times(starts));end
end

function summary=make_summary(tracker,cases,segments,diagnostics,cfg)
post=segments([segments.segment_index]==1);if isempty(post),post_duration=NaN;
else,post_duration=post.duration;end
exact_rows=diagnostics(isfinite([diagnostics.minimum_exact_edge_distance]));
if isempty(exact_rows),minimum_exact=NaN;maximum_shadow=NaN;exact_recovered=false;
else
    minimum_exact=min([exact_rows.minimum_exact_edge_distance]);
    maximum_shadow=max([exact_rows.exact_edge_shadow_duration]);
    exact_recovered=minimum_exact<cfg.dynamic_pair.exact_edge_neighborhood_threshold;
end
summary=struct('bracket_id',cfg.dynamic_pair.bracket_id, ...
    'left_outcome',cfg.dynamic_pair.left_outcome, ...
    'right_outcome',cfg.dynamic_pair.right_outcome,'status',tracker.status, ...
    'source_shadow_duration',segments(1).duration, ...
    'rebisection_steps',numel(cases), ...
    'renewed_pair_separation',tracker.restart_separation, ...
    'post_rebracketing_shadow_duration',post_duration, ...
    'minimum_exact_edge_distance',minimum_exact, ...
    'maximum_exact_edge_shadow_duration',maximum_shadow, ...
    'known_exact_edge_recovered',exact_recovered, ...
    'interpretation',cfg.dynamic_pair.claim_scope);
end

function record=empty_case_record(varargin)
record=struct('step',NaN,'alpha',NaN,'outcome','', ...
    'tested_duration',NaN,'decision_time',NaN,'runtime_seconds',NaN, ...
    'alpha_left_before',NaN,'alpha_right_before',NaN, ...
    'alpha_left_after',NaN,'alpha_right_after',NaN, ...
    'state_separation_after',NaN,'update','','raw_file','');
if nargin>0,record=record([]);end
end

function record=summarize_case(data,step,alpha,bL,bR,aL,aR,separation,update,raw)
record=empty_case_record();record.step=step;record.alpha=alpha;
record.outcome=data.classification.outcome;record.tested_duration=data.duration;
record.decision_time=data.decision_time;record.runtime_seconds=data.runtime_seconds;
record.alpha_left_before=bL;record.alpha_right_before=bR;
record.alpha_left_after=aL;record.alpha_right_after=aR;
record.state_separation_after=separation;record.update=update;record.raw_file=raw;
end

function record=empty_segment_record(varargin)
record=struct('segment_index',NaN,'global_start_time',NaN, ...
    'global_end_time',NaN,'duration',NaN,'initial_separation',NaN, ...
    'final_separation',NaN,'maximum_separation',NaN, ...
    'threshold_reached',false,'raw_file','');
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
[distance,~]=pair_distance_history(left(:,2)',left(:,3)', ...
    right(:,2)',right(:,3)',Lx);
end

function distance=pair_reflected_distance(left,right,Lx)
[~,distance]=pair_distance_history(left(:,2)',left(:,3)', ...
    right(:,2)',right(:,3)',Lx);
end

function make_figures(outdir,segments,cases,cfg)
fig=figure('Visible','off','Color','w');ax=axes(fig);hold(ax,'on');
for k=1:numel(segments)
    loaded=load(segments(k).raw_file,'segment_data');d=loaded.segment_data;
    plot(ax,segments(k).global_start_time+d.S,d.separation,'LineWidth',1.5, ...
        'DisplayName',sprintf('segment %d',segments(k).segment_index));
end
yline(ax,cfg.dynamic_pair.separation_threshold,'--k','threshold');
set(ax,'YScale','log');xlabel(ax,'assembled edge time');ylabel(ax,'pair separation');
title(ax,sprintf('%s dynamic edge tracking',cfg.dynamic_pair.bracket_id));
grid(ax,'on');legend(ax,'Location','best');
exportgraphics(fig,fullfile(outdir,'dynamic_edge_separation.png'),'Resolution',300);close(fig);
fig=figure('Visible','off','Color','w');ax=axes(fig);
if isempty(cases),plot(ax,NaN,NaN);else
    semilogy(ax,[cases.step],[cases.state_separation_after],'o-','LineWidth',1.4);
end
yline(ax,cfg.dynamic_pair.restart_separation_target,'--k','restart target');
xlabel(ax,'rebisection step');ylabel(ax,'renewed separation');grid(ax,'on');
title(ax,sprintf('%s five-way rebisection',cfg.dynamic_pair.bracket_id));
exportgraphics(fig,fullfile(outdir,'dynamic_rebisection.png'),'Resolution',300);close(fig);
end
