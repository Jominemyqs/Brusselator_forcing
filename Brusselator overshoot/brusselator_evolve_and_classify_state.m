function data = brusselator_evolve_and_classify_state( ...
        case_id,initial_state,library,cfg,partial_file)
%BRUSSELATOR_EVOLVE_AND_CLASSIFY_STATE Adaptively classify one frozen state.
%   Uses the existing physical five-way classifier without changing its
%   thresholds. Unresolved/ambiguous trajectories are extended in configured
%   chunks up to the dynamic-edge maximum duration.

if nargin<5,partial_file='';end
par=brusselator_make_parameters(cfg,@(t)cfg.edge_tracking.frozen_b);
maximum=cfg.dynamic_edge.classification_maximum_duration;
chunk=cfg.dynamic_edge.classification_chunk_duration;
known={library.label};
if ~isempty(partial_file)&&isfile(partial_file)
    loaded=load(partial_file,'partial_case');partial=loaded.partial_case;
    if ~strcmp(partial.case_id,case_id)|| ...
            max(abs(partial.initial_state(:)-initial_state(:)))>1e-12
        error('brusselator_evolve_and_classify_state:CheckpointMismatch', ...
            'Partial checkpoint does not match case %s.',case_id);
    end
    state=partial.state;S=partial.S;U=partial.U;V=partial.V;
    tested=partial.tested;classification=partial.classification;
    history=partial.history;accumulated_runtime=partial.accumulated_runtime;
    fprintf('  %s: resuming adaptive classification from t=%.0f.\n',case_id,S(end));
else
    state=initial_state;S=[];U=[];V=[];tested=[];classification=[];history=[];
    accumulated_runtime=0;
end
fprintf('  %s: adaptive five-way classification to at most t=%.0f ...\n', ...
    case_id,maximum);
while true
    previous=0;if ~isempty(S),previous=S(end);end
    if previous==0
        target=cfg.dynamic_edge.classification_initial_duration;
    else
        target=min(maximum,previous+chunk);
    end
    if (~isempty(classification)&&ismember(classification.outcome,known))|| ...
            previous>=maximum
        break;
    end
    duration=target-previous;chunk_timer=tic;
    [state,S2,V2,U2]=solve_brusselator_1d_forced(state,par,duration,0);
    [chunk_classification,chunk_history]= ...
        brusselator_classify_five_way_trajectory( ...
        S2,U2,V2,library,'late_window',cfg.edge_tracking.late_window, ...
        'sample_interval',cfg.edge_tracking.distance_sample_interval, ...
        'median_threshold',cfg.edge_tracking.median_distance_threshold, ...
        'maximum_threshold',cfg.edge_tracking.maximum_distance_threshold, ...
        'Lx',cfg.grid.Lx);
    if isempty(S)
        S=S2;U=U2;V=V2;
        history=chunk_history;
    else
        S=[S;previous+S2(2:end)]; %#ok<AGROW>
        U=[U;U2(2:end,:)];V=[V;V2(2:end,:)]; %#ok<AGROW>
        history.sample_times=[history.sample_times; ...
            previous+chunk_history.sample_times(2:end)];
        history.distances=[history.distances;chunk_history.distances(2:end,:)];
        history.uses_reflection=[history.uses_reflection; ...
            chunk_history.uses_reflection(2:end,:)];
        history.reference_phase_index=[history.reference_phase_index; ...
            chunk_history.reference_phase_index(2:end,:)];
    end
    classification=chunk_classification;
    classification.trajectory_end_time=target;
    tested(end+1,1)=target; %#ok<AGROW>
    accumulated_runtime=accumulated_runtime+toc(chunk_timer);
    history.sample_indices=(1:numel(history.sample_times))';
    fprintf('    tested t=%g: %s.\n',target,classification.outcome);
    if ~isempty(partial_file)
        partial_case=struct('case_id',case_id,'initial_state',initial_state, ...
            'state',state,'S',S,'U',U,'V',V,'tested',tested, ...
            'classification',classification,'history',history, ...
            'accumulated_runtime',accumulated_runtime,'complete', ...
            ismember(classification.outcome,known)||target>=maximum);
        save(partial_file,'partial_case','-v7.3');
    end
end
runtime=accumulated_runtime;
decision_time=first_commitment_time(history,classification.outcome,cfg);
recurrence=recurrence_check(S,U,V,classification.outcome,cfg);
fprintf('    outcome=%s, tested=%g, decision=%g, late d=[%s] (%.1f s).\n', ...
    classification.outcome,S(end),decision_time, ...
    sprintf(' %.3e',classification.late_median_distances),runtime);
data=struct('case_id',case_id,'initial_state',initial_state, ...
    'final_state',state,'S',S,'U',U,'V',V,'tested_durations',tested, ...
    'duration',S(end),'runtime_seconds',runtime,'decision_time',decision_time, ...
    'classification',classification,'history',history,'recurrence',recurrence);
end

function time=first_commitment_time(history,outcome,cfg)
time=NaN;index=find(strcmp(history.labels,outcome),1);if isempty(index),return;end
count=max(1,ceil(cfg.edge_tracking.decision_sustain_duration/ ...
    cfg.edge_tracking.distance_sample_interval));
for k=1:(numel(history.sample_times)-count+1)
    values=history.distances(k:k+count-1,index);
    if median(values)<cfg.edge_tracking.median_distance_threshold && ...
            max(values)<cfg.edge_tracking.maximum_distance_threshold
        time=history.sample_times(k);return;
    end
end
end

function recurrence=recurrence_check(S,U,V,outcome,cfg)
recurrence=struct('available',false,'period',NaN,'distance',NaN, ...
    'analysis',[],'error','');
if ~strcmp(outcome,'unresolved'),return;end
start=max(S(1),S(end)-cfg.edge_tracking.recurrence_analysis_window);
try
    analysis=brusselator_periodic_recurrence(S,U,V,start, ...
        cfg.edge_tracking.recurrence_threshold);
    recurrence.available=true;recurrence.analysis=analysis;
    best=analysis.recurrence.best;
    if isfinite(best.return_time)
        recurrence.period=best.return_time;
        recurrence.distance=best.best_median_rel_diff_full_state;
    else
        recurrence.period=analysis.section.period;
    end
catch exception
    recurrence.error=exception.message;
end
end
