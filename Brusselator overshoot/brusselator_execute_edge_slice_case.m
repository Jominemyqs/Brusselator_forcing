function data = brusselator_execute_edge_slice_case( ...
        case_id,lambda,duration,endpoint_B,endpoint_A,library,cfg)
%BRUSSELATOR_EXECUTE_EDGE_SLICE_CASE Evolve and classify one slice state.
%   Unresolved cases receive the configured adaptive extension. The returned
%   structure contains complete fields, physical five-way distances, decision
%   time, and an exploratory recurrence diagnostic when still unresolved.

x=endpoint_B(:,1);
u0=(1-lambda)*endpoint_B(:,2)+lambda*endpoint_A(:,2);
v0=(1-lambda)*endpoint_B(:,3)+lambda*endpoint_A(:,3);
initial_state=[x,u0,v0];
par=brusselator_make_parameters(cfg,@(t)cfg.edge_tracking.frozen_b);
fprintf('  %s: lambda=%.10f, frozen duration %.0f ...\n',case_id,lambda,duration);
timer=tic;
[final_state,S,V,U]=solve_brusselator_1d_forced(initial_state,par,duration,0);
[classification,history]=classify(S,U,V,library,cfg);
was_extended=false;
target=cfg.edge_tracking.unresolved_extension_duration;
if strcmp(classification.outcome,'unresolved') && duration<target
    fprintf('    unresolved at %.0f; extending to %.0f.\n',duration,target);
    [final_state,S2,V2,U2]=solve_brusselator_1d_forced( ...
        final_state,par,target-duration,0);
    S=[S;duration+S2(2:end)];U=[U;U2(2:end,:)];V=[V;V2(2:end,:)];
    duration=target;was_extended=true;
    [classification,history]=classify(S,U,V,library,cfg);
end
runtime=toc(timer);
decision_time=first_commitment_time(history,classification.outcome,cfg);
recurrence=recurrence_check(S,U,V,classification.outcome,cfg);
fprintf('    outcome=%s, decision=%g, late d=[%s] (%.1f s).\n', ...
    classification.outcome,decision_time, ...
    sprintf(' %.3e',classification.late_median_distances),runtime);
data=struct('case_id',case_id,'lambda',lambda,'duration',duration, ...
    'was_extended',was_extended,'runtime_seconds',runtime, ...
    'decision_time',decision_time,'initial_state',initial_state, ...
    'final_state',final_state,'S',S,'U',U,'V',V, ...
    'classification',classification,'history',history,'recurrence',recurrence);
end

function [classification,history]=classify(S,U,V,library,cfg)
[classification,history]=brusselator_classify_five_way_trajectory( ...
    S,U,V,library,'late_window',cfg.edge_tracking.late_window, ...
    'sample_interval',cfg.edge_tracking.distance_sample_interval, ...
    'median_threshold',cfg.edge_tracking.median_distance_threshold, ...
    'maximum_threshold',cfg.edge_tracking.maximum_distance_threshold, ...
    'Lx',cfg.grid.Lx);
end

function time=first_commitment_time(history,outcome,cfg)
time=NaN;class_index=find(strcmp(history.labels,outcome),1);
if isempty(class_index),return;end
count=max(1,ceil(cfg.edge_tracking.decision_sustain_duration/ ...
    cfg.edge_tracking.distance_sample_interval));
for k=1:(numel(history.sample_times)-count+1)
    values=history.distances(k:k+count-1,class_index);
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
