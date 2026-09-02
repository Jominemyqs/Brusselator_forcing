function results = run_BA_edge_tracking()
%RUN_BA_EDGE_TRACKING Full-state initial-condition edge tracking for B--A.
%   Uses the exact forcing-return states from the symmetric-ramp 40/45
%   bracket, aligns their reflection representatives, and bisects their full
%   (u,v) interpolation under the frozen b=10 PDE. A midpoint updates the
%   bracket only after unique A or B classification against the validated
%   A/B/C/R5/R6 orbit library. Other outcomes and persistent unresolved
%   dynamics stop the binary search and are retained as scientific results.

cfg=brusselator_edge_tracking_config();
study_id=cfg.experiment_name;
outdir=fullfile(cfg.output.root,study_id);
checkpoint_file=fullfile(outdir,'progress_edge_tracking.mat');
final_file=fullfile(outdir,'edge_tracking_results.mat');
if ~isfile(cfg.edge_tracking.stage2_file)
    error('run_BA_edge_tracking:MissingStage2', ...
        'Expected Stage 2 result at %s.',cfg.edge_tracking.stage2_file);
end
if isfile(final_file)
    error('run_BA_edge_tracking:OutputExists', ...
        'Refusing to overwrite completed edge-tracking output: %s',final_file);
end

source=load(cfg.edge_tracking.stage2_file,'results');
source=source.results;
validate_source_configuration(source,cfg);
library=build_library(source.outcomes);
[endpoint_B,endpoint_A,endpoint_info]=construct_endpoints(source,cfg);
frozen_par=brusselator_make_parameters(cfg,@(t)cfg.edge_tracking.frozen_b);

if isfile(checkpoint_file)
    saved=load(checkpoint_file,'records','raw_files','bracket','next_step', ...
        'status','candidate_case_id');
    records=saved.records;raw_files=saved.raw_files;bracket=saved.bracket;
    next_step=saved.next_step;status=saved.status;
    candidate_case_id=saved.candidate_case_id;
    fprintf('Resuming %s at bisection step %d.\n',study_id,next_step);
else
    if isfolder(outdir)
        error('run_BA_edge_tracking:IncompleteOutput', ...
            'Existing output lacks a valid checkpoint: %s',outdir);
    end
    mkdir(outdir);mkdir(fullfile(outdir,'raw_cases'));
    record_template=empty_record();records=record_template([]);raw_files={};
    bracket=struct('lambda_B',0,'lambda_A',1);
    next_step=1;status='initializing';candidate_case_id='';
    brusselator_write_metadata(outdir,cfg, ...
        brusselator_run_metadata(cfg,mfilename,NaN));
end
run_timer=tic;

[records,raw_files]=ensure_case(records,raw_files,'endpoint_B','endpoint',0, ...
    NaN,cfg.edge_tracking.classification_duration,endpoint_B,endpoint_A, ...
    library,cfg,frozen_par,outdir);
write_checkpoint(checkpoint_file,outdir,cfg,endpoint_info,records,raw_files, ...
    bracket,next_step,status,candidate_case_id);
[records,raw_files]=ensure_case(records,raw_files,'endpoint_A','endpoint',1, ...
    NaN,cfg.edge_tracking.classification_duration,endpoint_B,endpoint_A, ...
    library,cfg,frozen_par,outdir);
record_B=records(strcmp({records.case_id},'endpoint_B'));
record_A=records(strcmp({records.case_id},'endpoint_A'));
if ~strcmp(record_B.outcome,cfg.edge_tracking.left_expected_outcome) || ...
        ~strcmp(record_A.outcome,cfg.edge_tracking.right_expected_outcome)
    error('run_BA_edge_tracking:EndpointGate', ...
        'Fresh frozen endpoint classifications are %s and %s, expected B and A.', ...
        record_B.outcome,record_A.outcome);
end
if strcmp(status,'initializing')
    status='running_binary_bisection';
    write_checkpoint(checkpoint_file,outdir,cfg,endpoint_info,records,raw_files, ...
        bracket,next_step,status,candidate_case_id);
end

for step=next_step:cfg.edge_tracking.maximum_bisection_steps
    if ~strcmp(status,'running_binary_bisection')
        break;
    end
    lambda=0.5*(bracket.lambda_B+bracket.lambda_A);
    case_id=sprintf('bisection_%02d',step);
    fprintf('B--A edge bisection %d/%d: lambda=%.10f in [%.10f, %.10f].\n', ...
        step,cfg.edge_tracking.maximum_bisection_steps,lambda, ...
        bracket.lambda_B,bracket.lambda_A);
    before=bracket;
    [records,raw_files]=ensure_case(records,raw_files,case_id,'bisection', ...
        lambda,step,cfg.edge_tracking.classification_duration,endpoint_B, ...
        endpoint_A,library,cfg,frozen_par,outdir);
    index=find(strcmp({records.case_id},case_id),1);
    outcome=records(index).outcome;
    if strcmp(outcome,'B')
        bracket.lambda_B=lambda;
    elseif strcmp(outcome,'A')
        bracket.lambda_A=lambda;
    elseif ismember(outcome,{'C','R5','R6'})
        status=['third_outcome_',outcome,'_obstructs_binary_bracket'];
        candidate_case_id=case_id;
    elseif strcmp(outcome,'ambiguous_multiple_outcomes')
        status='ambiguous_midpoint_obstructs_binary_bracket';
        candidate_case_id=case_id;
    else
        status='persistent_unresolved_edge_candidate';
        candidate_case_id=case_id;
    end
    records(index)=set_bracket_fields(records(index),before,bracket);
    next_step=step+1;
    write_checkpoint(checkpoint_file,outdir,cfg,endpoint_info,records,raw_files, ...
        bracket,next_step,status,candidate_case_id);
end

if strcmp(status,'running_binary_bisection')
    status='finite_time_binary_bracket_refined';
end
if strcmp(status,'finite_time_binary_bracket_refined') && isempty(candidate_case_id)
    lambda=0.5*(bracket.lambda_B+bracket.lambda_A);
    candidate_case_id='final_midpoint_candidate';
    fprintf('Long candidate integration at lambda=%.10f for %.0f time units.\n', ...
        lambda,cfg.edge_tracking.candidate_duration);
    [records,raw_files]=ensure_case(records,raw_files,candidate_case_id, ...
        'long_candidate',lambda,cfg.edge_tracking.maximum_bisection_steps+1, ...
        cfg.edge_tracking.candidate_duration,endpoint_B,endpoint_A,library, ...
        cfg,frozen_par,outdir);
    index=find(strcmp({records.case_id},candidate_case_id),1);
    records(index)=set_bracket_fields(records(index),bracket,bracket);
    status=['long_candidate_',records(index).outcome];
    write_checkpoint(checkpoint_file,outdir,cfg,endpoint_info,records,raw_files, ...
        bracket,next_step,status,candidate_case_id);
end

summary=summarize(records,bracket,status,candidate_case_id,cfg);
writetable(struct2table(records,'AsArray',true), ...
    fullfile(outdir,'edge_tracking_summary.csv'));
writetable(struct2table(summary,'AsArray',true), ...
    fullfile(outdir,'edge_tracking_bracket.csv'));
make_figures(outdir,records,raw_files,candidate_case_id,library,cfg);
results=struct('configuration',cfg,'endpoint_info',endpoint_info, ...
    'endpoint_B',endpoint_B,'endpoint_A',endpoint_A,'library',library, ...
    'records',records,'raw_files',{raw_files},'bracket',bracket, ...
    'status',status,'candidate_case_id',candidate_case_id,'summary',summary);
save(final_file,'results','-v7.3');
brusselator_write_metadata(outdir,cfg, ...
    brusselator_run_metadata(cfg,mfilename,toc(run_timer)));
disp(struct2table(records,'AsArray',true));
disp(struct2table(summary,'AsArray',true));
fprintf('B--A edge-tracking pilot saved in: %s\n',outdir);
end

function validate_source_configuration(source,cfg)
if source.configuration.grid.N~=cfg.grid.N || ...
        abs(source.configuration.grid.Lx-cfg.grid.Lx)>1e-12 || ...
        abs(source.configuration.forcing.b0-cfg.edge_tracking.frozen_b)>1e-12
    error('run_BA_edge_tracking:SourceConfiguration', ...
        'Stage 2 grid, domain, or frozen parameter disagrees with edge config.');
end
end

function library=build_library(outcomes)
template=struct('label','','type','','period',NaN,'U',[],'V',[]);
library=repmat(template,numel(outcomes),1);
for k=1:numel(outcomes)
    library(k).label=outcomes(k).label;
    library(k).type=outcomes(k).type;
    library(k).period=outcomes(k).period;
    if strcmp(outcomes(k).type,'stationary')
        library(k).U=outcomes(k).orbit_U;
        library(k).V=outcomes(k).orbit_V;
    else
        library(k).U=outcomes(k).dense_template.U;
        library(k).V=outcomes(k).dense_template.V;
    end
end
if ~isequal({library.label},{'A','B','C','R5','R6'})
    error('run_BA_edge_tracking:OutcomeLibrary', ...
        'The validated Stage 2 outcome library is not A/B/C/R5/R6.');
end
end

function [B,A,info]=construct_endpoints(source,cfg)
records=source.landing_records;
iB=find(strcmp({records.protocol_id},cfg.edge_tracking.left_protocol),1);
iA=find(strcmp({records.protocol_id},cfg.edge_tracking.right_protocol),1);
if isempty(iB)||isempty(iA)
    error('run_BA_edge_tracking:MissingEndpoints', ...
        'The requested Stage 2 landing protocols are unavailable.');
end
uB=source.landing_states.raw_U(records(iB).state_index,:);
vB=source.landing_states.raw_V(records(iB).state_index,:);
uA=source.landing_states.raw_U(records(iA).state_index,:);
vA=source.landing_states.raw_V(records(iA).state_index,:);
direct=direct_distance(uA,vA,uB,vB,cfg.grid.Lx);
reflected=direct_distance(fliplr(uA),fliplr(vA),uB,vB,cfg.grid.Lx);
A_reflected=reflected<direct;
if A_reflected,uA=fliplr(uA);vA=fliplr(vA);end
x=linspace(0,cfg.grid.Lx,cfg.grid.N)';
B=[x,uB',vB'];A=[x,uA',vA'];
info=struct('B_protocol',records(iB).protocol_id, ...
    'A_protocol',records(iA).protocol_id,'B_state_index',records(iB).state_index, ...
    'A_state_index',records(iA).state_index,'A_reflected_for_alignment',A_reflected, ...
    'unaligned_endpoint_distance',direct,'reflected_endpoint_distance',reflected, ...
    'selected_endpoint_distance',min(direct,reflected), ...
    'interpolation',cfg.edge_tracking.interpolation);
end

function record=empty_record()
record=struct('case_id','','stage','','bisection_step',NaN,'lambda',NaN, ...
    'duration',NaN,'was_extended',false,'outcome','unclassified', ...
    'decision_time',NaN,'runtime_seconds',NaN, ...
    'initial_distance_A',NaN,'initial_distance_B',NaN, ...
    'initial_distance_C',NaN,'initial_distance_R5',NaN,'initial_distance_R6',NaN, ...
    'late_median_distance_A',NaN,'late_median_distance_B',NaN, ...
    'late_median_distance_C',NaN,'late_median_distance_R5',NaN, ...
    'late_median_distance_R6',NaN,'late_maximum_distance_A',NaN, ...
    'late_maximum_distance_B',NaN,'late_maximum_distance_C',NaN, ...
    'late_maximum_distance_R5',NaN,'late_maximum_distance_R6',NaN, ...
    'final_distance_A',NaN,'final_distance_B',NaN,'final_distance_C',NaN, ...
    'final_distance_R5',NaN,'final_distance_R6',NaN, ...
    'recurrence_available',false,'recurrence_period',NaN, ...
    'recurrence_distance',NaN,'bracket_B_before',NaN,'bracket_A_before',NaN, ...
    'bracket_B_after',NaN,'bracket_A_after',NaN,'bracket_width_after',NaN, ...
    'raw_file','');
end

function [records,raw_files]=ensure_case(records,raw_files,case_id,stage,lambda, ...
        step,duration,B,A,library,cfg,frozen_par,outdir)
existing=find(strcmp({records.case_id},case_id),1);
if ~isempty(existing)
    fprintf('Skipping saved edge case %s.\n',case_id);
    return;
end
raw_file=fullfile(outdir,'raw_cases',[case_id,'.mat']);
if isfile(raw_file)
    saved=load(raw_file,'case_result');
    record=saved.case_result.record;record.raw_file=raw_file;
    records(end+1,1)=record;raw_files{end+1,1}=raw_file;
    fprintf('Recovered completed edge case %s from its raw file.\n',case_id);
    return;
end
[record,case_result]=execute_case(case_id,stage,lambda,step,duration, ...
    B,A,library,cfg,frozen_par);
record.raw_file=raw_file;
case_result.record=record;
save(raw_file,'case_result','-v7.3');
records(end+1,1)=record;
raw_files{end+1,1}=raw_file;
end

function [record,data]=execute_case(case_id,stage,lambda,step,duration, ...
        B,A,library,cfg,frozen_par)
x=B(:,1);u0=(1-lambda)*B(:,2)+lambda*A(:,2);
v0=(1-lambda)*B(:,3)+lambda*A(:,3);
initial_state=[x,u0,v0];
fprintf('  %s: lambda=%.10f, frozen duration %.0f ...\n',case_id,lambda,duration);
timer=tic;
[final_state,S,V,U]=solve_brusselator_1d_forced( ...
    initial_state,frozen_par,duration,0);
[classification,history]=classify(S,U,V,library,cfg);
was_extended=false;
if strcmp(classification.outcome,'unresolved') && ...
        duration<cfg.edge_tracking.unresolved_extension_duration
    target=cfg.edge_tracking.unresolved_extension_duration;
    fprintf('    unresolved at %.0f; extending to %.0f.\n',duration,target);
    [final_state,S2,V2,U2]=solve_brusselator_1d_forced(final_state,frozen_par, ...
        target-duration,0);
    S=[S;duration+S2(2:end)];U=[U;U2(2:end,:)];V=[V;V2(2:end,:)];
    duration=target;was_extended=true;
    [classification,history]=classify(S,U,V,library,cfg);
end
runtime=toc(timer);
decision_time=first_commitment_time(history,classification.outcome,cfg);
recurrence=recurrence_check(S,U,V,classification.outcome,cfg);
record=empty_record();record.case_id=case_id;record.stage=stage;
record.bisection_step=step;record.lambda=lambda;record.duration=duration;
record.was_extended=was_extended;record.outcome=classification.outcome;
record.decision_time=decision_time;record.runtime_seconds=runtime;
record=set_distance_fields(record,history.distances(1,:), ...
    classification.late_median_distances,classification.late_maximum_distances, ...
    classification.final_distances);
record.recurrence_available=recurrence.available;
record.recurrence_period=recurrence.period;
record.recurrence_distance=recurrence.distance;
fprintf('    outcome=%s, decision=%g, late d=[%s] (%.1f s).\n', ...
    record.outcome,record.decision_time,sprintf(' %.3e', ...
    classification.late_median_distances),runtime);
data=struct('case_id',case_id,'stage',stage,'lambda',lambda, ...
    'initial_state',initial_state,'final_state',final_state,'S',S,'U',U,'V',V, ...
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
time=NaN;labels=history.labels;
class_index=find(strcmp(labels,outcome),1);
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
recurrence=struct('available',false,'period',NaN,'distance',NaN,'analysis',[], ...
    'error','');
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

function record=set_distance_fields(record,initial,medians,maxima,final)
names={'A','B','C','R5','R6'};
for k=1:numel(names)
    record.(['initial_distance_',names{k}])=initial(k);
    record.(['late_median_distance_',names{k}])=medians(k);
    record.(['late_maximum_distance_',names{k}])=maxima(k);
    record.(['final_distance_',names{k}])=final(k);
end
end

function record=set_bracket_fields(record,before,after)
record.bracket_B_before=before.lambda_B;record.bracket_A_before=before.lambda_A;
record.bracket_B_after=after.lambda_B;record.bracket_A_after=after.lambda_A;
record.bracket_width_after=after.lambda_A-after.lambda_B;
end

function write_checkpoint(file,outdir,cfg,endpoint_info,records,raw_files, ...
        bracket,next_step,status,candidate_case_id)
writetable(struct2table(records,'AsArray',true), ...
    fullfile(outdir,'progress_summary.csv'));
save(file,'cfg','endpoint_info','records','raw_files','bracket','next_step', ...
    'status','candidate_case_id','-v7.3');
end

function summary=summarize(records,bracket,status,candidate_case_id,cfg)
times=[records.decision_time];times(~isfinite(times))=[records(~isfinite(times)).duration];
[longest,index]=max(times);
summary=struct('number_cases',numel(records), ...
    'number_A',sum(strcmp({records.outcome},'A')), ...
    'number_B',sum(strcmp({records.outcome},'B')), ...
    'number_third_outcome',sum(ismember({records.outcome},{'C','R5','R6'})), ...
    'number_unresolved',sum(strcmp({records.outcome},'unresolved')), ...
    'lambda_B',bracket.lambda_B,'lambda_A',bracket.lambda_A, ...
    'bracket_width',bracket.lambda_A-bracket.lambda_B, ...
    'bracket_midpoint',0.5*(bracket.lambda_A+bracket.lambda_B), ...
    'status',status,'candidate_case_id',candidate_case_id, ...
    'longest_observed_lifetime',longest, ...
    'longest_lifetime_case_id',records(index).case_id, ...
    'frozen_b',cfg.edge_tracking.frozen_b, ...
    'interpretation',cfg.edge_tracking.interpretation);
end

function make_figures(outdir,records,raw_files,candidate_case_id,library,cfg)
outcomes={records.outcome};colors=zeros(numel(records),3);
for k=1:numel(records),colors(k,:)=outcome_color(outcomes{k});end
lifetimes=[records.decision_time];
lifetimes(~isfinite(lifetimes))=[records(~isfinite(lifetimes)).duration];
fig=figure('Color','w','Position',[100 100 900 520]);ax=axes(fig);hold(ax,'on');
scatter(ax,[records.lambda],lifetimes,65,colors,'filled','MarkerEdgeColor','k');
for k=1:numel(records)
    text(ax,records(k).lambda,lifetimes(k),[' ',records(k).case_id], ...
        'FontSize',7,'Interpreter','none');
end
xlabel(ax,'lambda: 0 = B landing, 1 = A landing');ylabel(ax,'decision time or tested duration');
title(ax,'Finite-time outcome lifetime along the B--A landing-state segment');
style_axes(ax);exportgraphics(fig,fullfile(outdir,'edge_lifetime_by_lambda.png'), ...
    'Resolution',300);close(fig);

candidate_index=find(strcmp({records.case_id},candidate_case_id),1);
if isempty(candidate_index)
    [~,candidate_index]=max(lifetimes);
end
loaded=load(raw_files{candidate_index},'case_result');data=loaded.case_result;
fig=figure('Color','w','Position',[100 100 950 580]);ax=axes(fig);hold(ax,'on');
for k=1:numel(library)
    semilogy(ax,data.history.sample_times,max(data.history.distances(:,k),eps), ...
        'LineWidth',1.3,'DisplayName',library(k).label);
end
yline(ax,cfg.edge_tracking.median_distance_threshold,'--k', ...
    'HandleVisibility','off');xlabel(ax,'frozen b=10 time');
ylabel(ax,'physical phase/reflection-aware distance');
title(ax,sprintf('Candidate %s at lambda=%.10f', ...
    records(candidate_index).case_id,records(candidate_index).lambda), ...
    'Interpreter','none');style_axes(ax);style_legend(legend(ax,'Location','best'));
exportgraphics(fig,fullfile(outdir,'candidate_distance_histories.png'),'Resolution',300);
close(fig);

fig=figure('Color','w','Position',[100 100 950 580]);ax=axes(fig);
imagesc(ax,linspace(0,cfg.grid.Lx,cfg.grid.N),data.S,data.V);axis(ax,'xy');
colorbar(ax);xlabel(ax,'x');ylabel(ax,'frozen b=10 time');
title(ax,sprintf('Candidate v(x,t): %s',records(candidate_index).case_id), ...
    'Interpreter','none');style_axes(ax);
exportgraphics(fig,fullfile(outdir,'candidate_v_spacetime.png'),'Resolution',300);
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
set(ax,'Color','w','XColor','k','YColor','k','ZColor','k', ...
    'GridColor',[0.75 0.75 0.75]);grid(ax,'on');box(ax,'on');
ax.Title.Color='k';ax.XLabel.Color='k';ax.YLabel.Color='k';
end

function style_legend(leg)
set(leg,'Color','w','TextColor','k','EdgeColor',[0.25 0.25 0.25], ...
    'Interpreter','none');
end

function distance=direct_distance(u,v,ref_u,ref_v,Lx)
N=numel(u);dx=Lx/(N-1);weights=dx*ones(1,N);weights([1,end])=0.5*dx;
distance=sqrt(sum(((u-ref_u).^2+(v-ref_v).^2).*weights)/ ...
    max(sum((ref_u.^2+ref_v.^2).*weights),eps));
end
