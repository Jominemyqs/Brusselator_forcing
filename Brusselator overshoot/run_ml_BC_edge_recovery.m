function results = run_ml_BC_edge_recovery(cfg)
%RUN_ML_BC_EDGE_RECOVERY Test whether frozen ML brackets recover exact E_BC.

if nargin<1
    cfg=brusselator_ml_edge_recovery_config();
end
outdir=fullfile(cfg.output.root,cfg.experiment_name);
final_file=fullfile(outdir,'ml_BC_edge_recovery.mat');
if isfile(final_file)
    error('run_ml_BC_edge_recovery:OutputExists', ...
        'Refusing to overwrite completed recovery result: %s',final_file);
end
if isfield(cfg.ml_edge,'frozen_manifest_file')
    frozen_manifest=cfg.ml_edge.frozen_manifest_file;
elseif isfield(cfg.ml_edge,'combined_pair_file')
    frozen_manifest=fullfile(cfg.ml_edge.model_directory,'frozen_pair_manifest.json');
else
    frozen_manifest=fullfile(cfg.ml_edge.model_directory,'frozen_model_manifest.json');
end
required={cfg.ml_labeling.geometry_file,cfg.ml_edge.endpoint_labels_file, ...
    cfg.ml_edge.exact_edge_file,frozen_manifest};
if ~all(cellfun(@isfile,required))
    error('run_ml_BC_edge_recovery:MissingSource', ...
        'Frozen models, endpoint labels, geometry, or exact E_BC are missing.');
end
timer=tic;
if ~isfolder(outdir)
    mkdir(outdir);mkdir(fullfile(outdir,'raw_cases'));
    brusselator_write_metadata(outdir,cfg, ...
        brusselator_run_metadata(cfg,mfilename,NaN));
elseif ~isfile(fullfile(outdir,'configuration_and_metadata.mat'))
    error('run_ml_BC_edge_recovery:UnsafeResume', ...
        'Existing output directory has no recognized metadata: %s',outdir);
end

G=load(cfg.ml_labeling.geometry_file,'results');geometry=G.results;
labels=readtable(cfg.ml_edge.endpoint_labels_file,'TextType','string');
[library,library_provenance]=brusselator_build_ml_outcome_library(cfg);
E=load(cfg.ml_edge.exact_edge_file,'results');edge=E.results;
if ~strcmp(edge.status,'newton_converged') || abs(edge.period-4.0029149135788)>1e-8
    error('run_ml_BC_edge_recovery:ExactEdgeGate', ...
        'The configured exact E_BC result is not the verified N=400 orbit.');
end
template=struct('U',edge.orbit_U,'V',edge.orbit_V, ...
    'period',edge.period,'phase_fraction',edge.orbit_times/edge.period);
model_results=repmat(failed_model_result(''),numel(cfg.ml_edge.model_names),1);
all_records=empty_bisection_record([]);
for model_index=1:numel(cfg.ml_edge.model_names)
    name=cfg.ml_edge.model_names{model_index};
    if isfield(cfg.ml_edge,'combined_pair_file')
        proposals=combined_pair_proposal(cfg.ml_edge.combined_pair_file,model_index);
    else
        proposal_file=fullfile(cfg.ml_edge.model_directory, ...
            [name,'_edge_bracket_proposals.csv']);
        proposals=readtable(proposal_file,'TextType','string');
    end
    selected=select_verified_bracket(proposals,labels,name);
    if isempty(selected)
        fprintf('\nML edge recovery for %s model: no verified B/C pair among four frozen proposals.\n',name);
        model_result=failed_model_result(name);
        model_results(model_index,1)=model_result;
        continue;
    end
    [B,C]=endpoint_states(selected,labels,geometry);
    fprintf('\nML edge recovery for %s model: proposal %d, states %d/%d.\n', ...
        name,selected.proposal_rank,selected.first_state_index, ...
        selected.second_state_index);
    [model_result,records]=bisect_and_diagnose(name,B,C,selected, ...
        library,template,cfg,outdir);
    model_results(model_index,1)=model_result;
    all_records=[all_records;records(:)]; %#ok<AGROW>
    writetable(struct2table(all_records,'AsArray',true), ...
        fullfile(outdir,'bisection_progress.csv'));
end

function proposals=combined_pair_proposal(filename,pair_index)
batch=readtable(filename,'TextType','string');
if pair_index>height(batch)
    error('run_ml_BC_edge_recovery:PairIndex', ...
        'Requested pair %d but the frozen batch has only %d rows.', ...
        pair_index,height(batch));
end
row=batch(pair_index,:);
ids=[row.first_state_index,row.second_state_index];
if all(ismember({'first_probability_C','second_probability_C'}, ...
        row.Properties.VariableNames))
    probability=[row.first_probability_C,row.second_probability_C];
elseif all(ismember({'first_latent_mean','second_latent_mean'}, ...
        row.Properties.VariableNames))
    latent=[row.first_latent_mean,row.second_latent_mean];
    probability=0.5*(1+erf(latent/sqrt(2)));
else
    error('run_ml_BC_edge_recovery:PairProbability', ...
        'Combined pair table lacks endpoint probabilities or latent means.');
end
[~,order]=sort(probability,'ascend');
if ismember('proposal_rank',row.Properties.VariableNames)
    rank=row.proposal_rank;
else
    rank=row.batch_step;
end
proposals=table(rank,ids(order(1)),ids(order(2)), ...
    probability(order(1)),probability(order(2)), ...
    'VariableNames',{'proposal_rank','B_state_index','C_state_index', ...
    'B_probability_C','C_probability_C'});
end
summary=make_summary(model_results,cfg);
writetable(struct2table(all_records,'AsArray',true), ...
    fullfile(outdir,'bisection_cases.csv'));
writetable(struct2table(summary,'AsArray',true), ...
    fullfile(outdir,'edge_recovery_summary.csv'));
make_figure(outdir,model_results,cfg);
results=struct('configuration',cfg,'library_provenance',library_provenance, ...
    'exact_edge_period',edge.period,'models',model_results, ...
    'bisection_records',all_records,'summary',summary);
save(final_file,'results','-v7.3');
brusselator_write_metadata(outdir,cfg, ...
    brusselator_run_metadata(cfg,mfilename,toc(timer)));
disp(struct2table(summary,'AsArray',true));
fprintf('ML edge-recovery study saved in: %s\n',outdir);
end

function selected=select_verified_bracket(proposals,labels,name)
for k=1:height(proposals)
    ids=[proposals.B_state_index(k),proposals.C_state_index(k)];
    first=labels.final_label(labels.state_index==ids(1));
    second=labels.final_label(labels.state_index==ids(2));
    if isscalar(first) && isscalar(second) && ...
            isequal(sort([first,second]),sort(["B","C"]))
        selected=struct('model',name,'proposal_rank',proposals.proposal_rank(k), ...
            'first_state_index',ids(1),'second_state_index',ids(2), ...
            'first_actual_label',char(first),'second_actual_label',char(second), ...
            'first_predicted_probability_C',proposals.B_probability_C(k), ...
            'second_predicted_probability_C',proposals.C_probability_C(k));
        return;
    end
end
selected=[];
end

function [B,C]=endpoint_states(selected,labels,geometry)
ids=[selected.first_state_index,selected.second_state_index];
states=cell(2,1);
for k=1:2
    states{k}=[geometry.center.x(:),geometry.states.raw_U(ids(k),:)', ...
        geometry.states.raw_V(ids(k),:)'];
end
first=labels.final_label(labels.state_index==ids(1));
if first=="B",B=states{1};C=states{2};else,B=states{2};C=states{1};end
end

function [result,records]=bisect_and_diagnose(name,B,C,selected,library,template,cfg,outdir)
records=empty_bisection_record([]);alpha_B=0;alpha_C=1;
state_B=B;state_C=C;status='binary_bisection_complete';
for step=1:cfg.ml_edge.maximum_bisection_steps
    alpha=0.5*(alpha_B+alpha_C);midpoint=interpolate(state_B,state_C,0.5);
    case_id=sprintf('%s_bisection_%02d',name,step);
    raw_file=fullfile(outdir,'raw_cases',[case_id,'.mat']);
    if isfile(raw_file)
        loaded=load(raw_file,'case_data');case_data=loaded.case_data;
    else
        partial_file=[raw_file(1:end-4),'_partial.mat'];
        case_data=brusselator_evolve_and_classify_state( ...
            case_id,midpoint,library,cfg,partial_file);
        save(raw_file,'case_data','-v7.3');
        if isfile(partial_file),delete(partial_file);end
    end
    outcome=case_data.classification.outcome;
    before_B=alpha_B;before_C=alpha_C;
    if strcmp(outcome,'B')
        state_B=midpoint;alpha_B=alpha;
    elseif strcmp(outcome,'C')
        state_C=midpoint;alpha_C=alpha;
    else
        status=['obstructed_by_',outcome];
    end
    record=empty_bisection_record();record.model=name;record.step=step;
    record.alpha=alpha;record.outcome=outcome;record.alpha_B_before=before_B;
    record.alpha_C_before=before_C;record.alpha_B_after=alpha_B;
    record.alpha_C_after=alpha_C;record.bracket_width=alpha_C-alpha_B;
    record.tested_duration=case_data.duration;record.decision_time=case_data.decision_time;
    record.runtime_seconds=case_data.runtime_seconds;record.raw_file=raw_file;
    records(end+1,1)=record; %#ok<AGROW>
    fprintf('  %s step %d: %s, width %.4e.\n',name,step,outcome,alpha_C-alpha_B);
    if ~strcmp(status,'binary_bisection_complete'),break;end
end
candidate=interpolate(state_B,state_C,0.5);
candidate_file=fullfile(outdir,'raw_cases',[name,'_candidate.mat']);
if isfile(candidate_file)
    loaded=load(candidate_file,'candidate_data');candidate_data=loaded.candidate_data;
else
    par=brusselator_make_parameters(cfg,@(t)cfg.edge_tracking.frozen_b);
    solve_timer=tic;
    [final_state,S,V,U]=solve_brusselator_1d_forced( ...
        candidate,par,cfg.ml_edge.candidate_duration,0);
    [classification,history]=brusselator_classify_five_way_trajectory( ...
        S,U,V,library,'late_window',cfg.edge_tracking.late_window, ...
        'sample_interval',cfg.edge_tracking.distance_sample_interval, ...
        'median_threshold',cfg.edge_tracking.median_distance_threshold, ...
        'maximum_threshold',cfg.edge_tracking.maximum_distance_threshold, ...
        'Lx',cfg.grid.Lx);
    [distance_times,distances]=edge_distances(S,U,V,template,cfg);
    [shadow_start,shadow_end,shadow_duration]=longest_shadow( ...
        distance_times,distances,cfg.ml_edge.edge_neighborhood_threshold);
    period=period_diagnostic(S,U,V,shadow_start,shadow_end,template.period);
    candidate_data=struct('initial_state',candidate,'final_state',final_state, ...
        'S',S,'U',U,'V',V,'classification',classification,'history',history, ...
        'distance_times',distance_times,'edge_orbit_distances',distances, ...
        'shadow_start',shadow_start,'shadow_end',shadow_end, ...
        'shadow_duration',shadow_duration,'period_diagnostic',period, ...
        'runtime_seconds',toc(solve_timer));
    save(candidate_file,'candidate_data','-v7.3');
end
[minimum_distance,index]=min(candidate_data.edge_orbit_distances);
period_value=candidate_data.period_diagnostic.best_period;
period_error=abs(period_value-template.period)/template.period;
result=struct('model',name,'selected_proposal_rank',selected.proposal_rank, ...
    'selected_first_state_index',selected.first_state_index, ...
    'selected_second_state_index',selected.second_state_index, ...
    'verified_BC_bracket_found',true, ...
    'bisection_status',status,'bisection_steps',numel(records), ...
    'final_bracket_width',alpha_C-alpha_B, ...
    'candidate_final_outcome',candidate_data.classification.outcome, ...
    'minimum_edge_orbit_distance',minimum_distance, ...
    'minimum_distance_time',candidate_data.distance_times(index), ...
    'edge_shadow_duration',candidate_data.shadow_duration, ...
    'direct_period_estimate',period_value,'period_relative_error',period_error, ...
    'period_diagnostic_available',candidate_data.period_diagnostic.available, ...
    'recovered_exact_edge_neighborhood', ...
        minimum_distance<cfg.ml_edge.edge_neighborhood_threshold, ...
    'period_matches_when_available',~candidate_data.period_diagnostic.available || ...
        period_error<cfg.ml_edge.period_relative_tolerance, ...
    'candidate_raw_file',candidate_file,'distance_times',candidate_data.distance_times, ...
    'edge_orbit_distances',candidate_data.edge_orbit_distances);
end

function record=empty_bisection_record(varargin)
record=struct('model','','step',NaN,'alpha',NaN,'outcome','', ...
    'alpha_B_before',NaN,'alpha_C_before',NaN,'alpha_B_after',NaN, ...
    'alpha_C_after',NaN,'bracket_width',NaN,'tested_duration',NaN, ...
    'decision_time',NaN,'runtime_seconds',NaN,'raw_file','');
if nargin>0,record=record([]);end
end

function state=interpolate(B,C,alpha)
state=B;state(:,2:3)=(1-alpha)*B(:,2:3)+alpha*C(:,2:3);
end

function [times,distances]=edge_distances(S,U,V,template,cfg)
times=(S(1):cfg.ml_edge.orbit_distance_sample_interval:S(end))';
sample_U=interp1(S,U,times,'pchip');sample_V=interp1(S,V,times,'pchip');
distances=zeros(size(times));
for k=1:numel(times)
    item=brusselator_physical_orbit_distance(sample_U(k,:),sample_V(k,:), ...
        template.U,template.V,cfg.grid.Lx);
    distances(k)=item.distance;
end
end

function [start_time,end_time,duration]=longest_shadow(times,distances,threshold)
mask=distances<threshold;change=diff([false;mask;false]);
starts=find(change==1);ends=find(change==-1)-1;
if isempty(starts),start_time=NaN;end_time=NaN;duration=0;return;end
[duration,index]=max(times(ends)-times(starts));
start_time=times(starts(index));end_time=times(ends(index));
end

function period=period_diagnostic(S,U,V,start_time,end_time,guess)
period=struct('available',false,'best_period',NaN,'best_distance',NaN, ...
    'best_uses_reflection',false,'error','');
if ~isfinite(start_time)||end_time-start_time<50,return;end
selected=S>=start_time&S<=end_time;
S2=S(selected)-S(find(selected,1));U2=U(selected,:);V2=V(selected,:);
try
    item=brusselator_direct_period_recurrence(S2,U2,V2,S2(1),guess);
    period.available=true;period.best_period=item.best_period;
    period.best_distance=item.best_distance;
    period.best_uses_reflection=item.best_uses_reflection;
catch exception
    period.error=exception.message;
end
end

function summary=make_summary(models,cfg)
summary=repmat(struct('model','','verified_BC_bracket_found',false, ...
    'recovered_exact_edge_neighborhood',false,'minimum_edge_orbit_distance',NaN, ...
    'edge_shadow_duration',NaN,'direct_period_estimate',NaN, ...
    'exact_edge_period',4.0029149135788,'period_relative_error',NaN, ...
    'candidate_final_outcome','','interpretation',''),numel(models),1);
for k=1:numel(models)
    summary(k).model=models(k).model;
    summary(k).verified_BC_bracket_found=models(k).verified_BC_bracket_found;
    summary(k).recovered_exact_edge_neighborhood=models(k).recovered_exact_edge_neighborhood;
    summary(k).minimum_edge_orbit_distance=models(k).minimum_edge_orbit_distance;
    summary(k).edge_shadow_duration=models(k).edge_shadow_duration;
    summary(k).direct_period_estimate=models(k).direct_period_estimate;
    summary(k).period_relative_error=models(k).period_relative_error;
    summary(k).candidate_final_outcome=models(k).candidate_final_outcome;
    summary(k).interpretation=cfg.ml_edge.claim_scope;
end
end

function result=failed_model_result(name)
result=struct('model',name,'selected_proposal_rank',NaN, ...
    'selected_first_state_index',NaN,'selected_second_state_index',NaN, ...
    'verified_BC_bracket_found',false, ...
    'bisection_status','no_verified_BC_bracket_in_four_label_blind_proposals', ...
    'bisection_steps',0,'final_bracket_width',NaN, ...
    'candidate_final_outcome','not_run','minimum_edge_orbit_distance',NaN, ...
    'minimum_distance_time',NaN,'edge_shadow_duration',NaN, ...
    'direct_period_estimate',NaN,'period_relative_error',NaN, ...
    'period_diagnostic_available',false,'recovered_exact_edge_neighborhood',false, ...
    'period_matches_when_available',false,'candidate_raw_file','', ...
    'distance_times',[],'edge_orbit_distances',[]);
end

function make_figure(outdir,models,cfg)
fig=figure('Visible','off','Color','w');ax=axes(fig);hold(ax,'on');
for k=1:numel(models)
    if isempty(models(k).distance_times),continue;end
    semilogy(ax,models(k).distance_times,max(models(k).edge_orbit_distances,eps), ...
        'LineWidth',1.4,'DisplayName',models(k).model);
end
yline(ax,cfg.ml_edge.edge_neighborhood_threshold,'--k','edge neighborhood');
xlabel(ax,'frozen time');ylabel(ax,'distance to exact E_{BC} orbit');
title(ax,'Dynamical usefulness of label-blind ML brackets');grid(ax,'on');
legend(ax,'Location','best');
exportgraphics(fig,fullfile(outdir,'ML_bracket_edge_recovery.png'),'Resolution',300);
close(fig);
end
