function results = build_stage2_landing_geometry()
%BUILD_STAGE2_LANDING_GEOMETRY Construct weighted, phase-aware POD geometry.
%   Builds an auditable five-outcome orbit library and deduplicated forcing
%   landing-state manifest. Raw, reflected, feature-canonical, and
%   family-continuity-oriented states are preserved. The primary POD uses
%   explicit equal aggregate outcome weights and separate u/v scales; an
%   unscaled-component POD is computed as a sensitivity analysis.

study_id = 'stage2_landing_geometry_v2';
run_timer = tic;
outdir = fullfile('experiment_outputs',study_id);
if isfolder(outdir)
    error('build_stage2_landing_geometry:OutputExists', ...
        'Refusing to overwrite existing Stage 2 output: %s',outdir);
end
cfg = brusselator_stage2_geometry_config();
validate_source_gates(cfg);
[outcomes,phase_records] = load_outcome_library(cfg);
[landing_records,landing_states] = load_landing_states(cfg);
[landing_records,landing_states,family_paths,family_membership] = ...
    orient_landing_states(landing_records,landing_states,cfg);
[landing_records,physical_distances] = add_physical_distances( ...
    landing_records,landing_states,outcomes,cfg);

[training_U,training_V,training_weights,training_labels] = ...
    build_weighted_training_set(outcomes,cfg);
primary = brusselator_weighted_pod(training_U,training_V,training_weights, ...
    cfg.grid.Lx,true,cfg.stage2.pod.maximum_components);
unscaled = brusselator_weighted_pod(training_U,training_V,training_weights, ...
    cfg.grid.Lx,false,cfg.stage2.pod.maximum_components);
[projections,reconstruction_records] = project_all_states( ...
    outcomes,landing_states,family_paths,training_U,training_V, ...
    training_labels,primary,unscaled,cfg);
scaling_sensitivity = compare_scalings(primary,unscaled, ...
    projections.landing.feature_canonical.primary, ...
    projections.landing.feature_canonical.unscaled,cfg);
[landing_aware,landing_aware_diagnostic] = build_landing_aware_diagnostic( ...
    outcomes,landing_states,training_U,training_V,training_weights, ...
    projections,primary,cfg);

mkdir(outdir);
brusselator_write_metadata(outdir,cfg, ...
    brusselator_run_metadata(cfg,mfilename,NaN));
writetable(struct2table(landing_records,'AsArray',true), ...
    fullfile(outdir,'landing_state_manifest.csv'));
writetable(struct2table(phase_records,'AsArray',true), ...
    fullfile(outdir,'phase_section_diagnostics.csv'));
writetable(struct2table(family_membership,'AsArray',true), ...
    fullfile(outdir,'family_orientation_manifest.csv'));
writetable(struct2table(reconstruction_records,'AsArray',true), ...
    fullfile(outdir,'pod_reconstruction_summary.csv'));
writetable(struct2table(scaling_sensitivity,'AsArray',true), ...
    fullfile(outdir,'pod_scaling_sensitivity.csv'));
write_landing_aware_tables(outdir,landing_aware,landing_aware_diagnostic,cfg);
write_pod_tables(outdir,primary,unscaled,training_labels,training_weights);
make_figures(outdir,outcomes,landing_records,physical_distances, ...
    family_paths,projections,primary,unscaled,reconstruction_records, ...
    landing_aware,landing_aware_diagnostic);

results = struct('configuration',cfg,'outcomes',outcomes, ...
    'phase_records',phase_records,'landing_records',landing_records, ...
    'landing_states',landing_states,'family_paths',{family_paths}, ...
    'family_membership',family_membership, ...
    'physical_distances',physical_distances, ...
    'training',struct('U',training_U,'V',training_V, ...
        'weights',training_weights,'labels',{training_labels}), ...
    'primary_pod',primary,'unscaled_pod',unscaled, ...
    'projections',projections, ...
    'reconstruction_records',reconstruction_records, ...
    'scaling_sensitivity',scaling_sensitivity, ...
    'landing_aware_pod',landing_aware, ...
    'landing_aware_diagnostic',landing_aware_diagnostic);
save(fullfile(outdir,'stage2_landing_geometry.mat'),'results','-v7.3');
brusselator_write_metadata(outdir,cfg, ...
    brusselator_run_metadata(cfg,mfilename,toc(run_timer)));
disp(struct2table(phase_records,'AsArray',true));
disp(struct2table(scaling_sensitivity,'AsArray',true));
fprintf('Stage 2 landing geometry saved in: %s\n',outdir);
end

function validate_source_gates(cfg)
fields = fieldnames(cfg.stage2.sources);
for k = 1:numel(fields)
    value = cfg.stage2.sources.(fields{k});
    if ischar(value) && ~strcmp(fields{k},'ABC_reference_loader') && ~isfile(value)
        error('build_stage2_landing_geometry:MissingSource', ...
            'Missing Stage 2 source %s: %s',fields{k},value);
    end
end
b = load(cfg.stage2.sources.B_local_gate,'results');
if ~all(contains(string({b.results.records.classification}),'returns'))
    error('build_stage2_landing_geometry:BGateFailed', ...
        'The saved B local-return gate is not fully supporting.');
end
c = load(cfg.stage2.sources.C_numerical_gate,'results');
if ~c.results.summary.proceed_to_L40_forcing_scan
    error('build_stage2_landing_geometry:CGateFailed', ...
        'The saved C fixed-Lx numerical gate failed.');
end
c_local = load(cfg.stage2.sources.C_local_gate,'results');
if ~all(contains(string({c_local.results.records.classification}),'returns'))
    error('build_stage2_landing_geometry:CLocalGateFailed', ...
        'The saved C local-return gate is not fully supporting.');
end
r = load(cfg.stage2.sources.recurrent_local_gate,'results');
if ~all([r.results.gates.finite_time_distinct_attractor_supported])
    error('build_stage2_landing_geometry:RecurrentGateFailed', ...
        'R5/R6 have not passed their saved finite-time evidence gates.');
end
end

function [outcomes,records] = load_outcome_library(cfg)
[refs,~] = brusselator_three_way_references();
slice = load(refs.info.slice_file,'results');
slice = slice.results;
b_index = find(abs([slice.records.lambda]-1)<1e-12,1);
if isempty(b_index)
    error('build_stage2_landing_geometry:MissingBTrajectory', ...
        'No lambda=1 B trajectory is saved.');
end
c_source = load(refs.info.candidate_file,'results');
c_source = c_source.results;
validation = load(cfg.stage2.sources.recurrent_validation,'results');
validation = validation.results;

template = struct('label','','type','','period',NaN,'S',[],'U',[],'V',[], ...
    'dense_template',struct(),'orbit_U',[],'orbit_V',[], ...
    'phase_fixed',struct(),'orbit_uses_reflection',false, ...
    'phase_diagnostic',struct());
outcomes = repmat(template,5,1);
outcomes(1).label = 'A'; outcomes(1).type = 'stationary';
orientation = brusselator_reflection_canonicalize(refs.A.U,refs.A.V,cfg.grid.Lx);
outcomes(1).orbit_U = orientation.canonical_U;
outcomes(1).orbit_V = orientation.canonical_V;
outcomes(1).phase_fixed = struct('times',NaN,'raw_U',refs.A.U, ...
    'raw_V',refs.A.V,'canonical_U',orientation.canonical_U, ...
    'canonical_V',orientation.canonical_V, ...
    'uses_reflection',orientation.uses_reflection,'branch_id',0, ...
    'distance_to_template_branch',0);

b_traj = slice.trajectory_data{b_index};
outcomes(2) = periodic_outcome(template,'B',refs.B.phase_times(end), ...
    b_traj.S,b_traj.U,b_traj.V,cfg);
outcomes(3) = periodic_outcome(template,'C',refs.C.phase_times(end), ...
    c_source.trajectory.S,c_source.trajectory.U,c_source.trajectory.V,cfg);
classes = {'period_5p0585','period_6p5837'};
labels = {'R5','R6'};
for k = 1:2
    index = find(strcmp({validation.validation_records.class_id},classes{k}) & ...
        [validation.validation_records.N]==400 & ...
        strcmp({validation.validation_records.factor},'tolerance'),1);
    loaded = load(validation.validation_raw_files{index},'case_result');
    item = loaded.case_result;
    outcomes(k+3) = periodic_outcome(template,labels{k}, ...
        item.section.section.period,item.S,item.U,item.V,cfg);
end

record_template = struct('label','','type','','period',NaN, ...
    'number_crossings_per_period',NaN,'minimum_relative_transversality',NaN, ...
    'transversality_pass',false,'temporal_resampling_counts','', ...
    'maximum_temporal_phase_error',NaN,'temporal_resampling_pass',false, ...
    'spatial_resampling_count',NaN,'spatial_phase_error',NaN, ...
    'spatial_resampling_pass',false, ...
    'maximum_branch_distance_modulo_reflection',NaN, ...
    'branches_reflection_equivalent',false,'branch_policy','', ...
    'number_late_source_crossings',NaN,'expected_late_source_crossings',NaN, ...
    'source_crossing_count_pass',false,'maximum_source_branch_distance',NaN, ...
    'branch_reproducibility_pass',false,'phase_section_pass',false);
records = repmat(record_template,5,1);
records(1) = record_template; records(1).label='A'; records(1).type='stationary';
records(1).phase_section_pass = true; records(1).branch_policy='not_applicable';
for k = 2:5
    d = outcomes(k).phase_diagnostic;
    records(k) = record_template;
    records(k).label = outcomes(k).label; records(k).type='periodic';
    records(k).period = outcomes(k).period;
    records(k).number_crossings_per_period = d.number_crossings_per_period;
    records(k).minimum_relative_transversality = d.minimum_relative_transversality;
    records(k).transversality_pass = d.transversality_pass;
    records(k).temporal_resampling_counts = mat2str(d.temporal_resampling_counts);
    records(k).maximum_temporal_phase_error = max(d.temporal_phase_errors);
    records(k).temporal_resampling_pass = d.temporal_resampling_pass;
    records(k).spatial_resampling_count = d.spatial_resampling_count;
    records(k).spatial_phase_error = d.spatial_phase_error;
    records(k).spatial_resampling_pass = d.spatial_resampling_pass;
    records(k).maximum_branch_distance_modulo_reflection = ...
        d.maximum_branch_distance_modulo_reflection;
    records(k).branches_reflection_equivalent = d.branches_reflection_equivalent;
    records(k).branch_policy = d.branch_policy;
    records(k).number_late_source_crossings = d.number_late_source_crossings;
    records(k).expected_late_source_crossings = d.expected_late_source_crossings;
    records(k).source_crossing_count_pass = d.source_crossing_count_pass;
    records(k).maximum_source_branch_distance = d.maximum_source_branch_distance;
    records(k).branch_reproducibility_pass = d.branch_reproducibility_pass;
    records(k).phase_section_pass = d.phase_section_pass;
end
if ~all([records.phase_section_pass])
    failed = strjoin({records(~[records.phase_section_pass]).label},', ');
    error('build_stage2_landing_geometry:PhaseSectionFailed', ...
        'Phase-section validation failed for: %s',failed);
end
end

function outcome = periodic_outcome(template,label,period,S,U,V,cfg)
outcome = template; outcome.label=label; outcome.type='periodic';
outcome.period=period; outcome.S=S(:); outcome.U=U; outcome.V=V;
outcome.dense_template = dense_template(S,U,V,period, ...
    cfg.stage2.phase.dense_template_samples);
[diagnostic,phase_fixed] = brusselator_phase_section_diagnostics( ...
    outcome.dense_template,S,U,V,cfg.grid.Lx,cfg.stage2.phase);
outcome.phase_diagnostic = diagnostic;
outcome.phase_fixed = phase_fixed;
anchor = brusselator_reflection_canonicalize( ...
    phase_fixed.raw_U(1,:),phase_fixed.raw_V(1,:),cfg.grid.Lx);
outcome.orbit_uses_reflection = anchor.uses_reflection;
n = cfg.stage2.phase.full_orbit_plot_samples;
indices = floor((0:n-1)*size(outcome.dense_template.U,1)/n)+1;
outcome.orbit_U = outcome.dense_template.U(indices,:);
outcome.orbit_V = outcome.dense_template.V(indices,:);
if outcome.orbit_uses_reflection
    outcome.orbit_U = fliplr(outcome.orbit_U);
    outcome.orbit_V = fliplr(outcome.orbit_V);
end
end

function template = dense_template(S,U,V,period,count)
S=S(:); start=S(end)-period; phase=(0:count-1)'/count;
query=start+phase*period;
template=struct('phase_fraction',phase,'phase_times',phase*period, ...
    'period',period,'start_time',start, ...
    'U',interp1(S,U,query,'pchip'),'V',interp1(S,V,query,'pchip'));
end

function [records,states] = load_landing_states(cfg)
pilot=load(cfg.stage2.sources.forcing_pilot,'results'); pilot=pilot.results;
refinement=load(cfg.stage2.sources.ramp_refinement,'results'); refinement=refinement.results;
validation=load(cfg.stage2.sources.recurrent_validation,'results'); validation=validation.results;
sources={pilot,refinement}; source_names={'forcing_pilot','ramp_refinement'};
template=struct('state_index',NaN,'source_study','','source_raw_file','', ...
    'protocol_id','','family','','bmax',NaN,'Tup',NaN,'Thold',NaN,'Tdown',NaN, ...
    'forcing_end',NaN,'landing_time',NaN,'grid_N',NaN,'grid_Lx',NaN, ...
    'original_outcome','','final_label','', ...
    'label_status','','bracket_role','','feature_uses_reflection',false, ...
    'feature_index',NaN,'feature_value',NaN,'reflection_invariance_error',NaN, ...
    'distance_A',NaN,'distance_B',NaN, ...
    'distance_C',NaN,'distance_R5',NaN,'distance_R6',NaN);
records=template([]); raw_U=zeros(0,cfg.grid.N); raw_V=raw_U; keys={};
for s=1:numel(sources)
    source=sources{s};
    for k=1:numel(source.records)
        loaded=load(source.raw_files{k},'case_result'); item=loaded.case_result;
        protocol=item.protocol; discovery=item.discovery;
        validate_landing_source(item,cfg,source.raw_files{k});
        key=sprintf('%.10g|%.10g|%.10g|%.10g',protocol.bmax,protocol.Tup, ...
            protocol.Thold,protocol.Tdown);
        duplicate=find(strcmp(keys,key),1);
        u=discovery.landing_state(:,2)'; v=discovery.landing_state(:,3)';
        if ~isempty(duplicate)
            mismatch=direct_state_distance(u,v,raw_U(duplicate,:),raw_V(duplicate,:),cfg.grid.Lx);
            if mismatch>1e-10
                error('build_stage2_landing_geometry:DuplicateMismatch', ...
                    'Duplicate forcing tuple %s disagrees by %.3e.',key,mismatch);
            end
            continue;
        end
        keys{end+1,1}=key; %#ok<AGROW>
        raw_U(end+1,:)=u; raw_V(end+1,:)=v; %#ok<AGROW>
        record=template; record.state_index=size(raw_U,1);
        record.source_study=source_names{s}; record.source_raw_file=source.raw_files{k};
        record.protocol_id=protocol.id; record.family=normalize_family(protocol.family);
        record.bmax=protocol.bmax; record.Tup=protocol.Tup;
        record.Thold=protocol.Thold; record.Tdown=protocol.Tdown;
        record.forcing_end=discovery.forcing_end;
        record.landing_time=discovery.landing_time;
        record.grid_N=discovery.configuration.grid.N;
        record.grid_Lx=discovery.configuration.grid.Lx;
        record.original_outcome=discovery.classification.outcome;
        [record.final_label,record.label_status]=landing_label( ...
            protocol.id,record.original_outcome,validation.cross_preparation_records);
        record.bracket_role=bracket_role(record);
        records(end+1,1)=record; %#ok<AGROW>
    end
end
states=struct('raw_U',raw_U,'raw_V',raw_V, ...
    'reflected_U',fliplr(raw_U),'reflected_V',fliplr(raw_V), ...
    'feature_canonical_U',[],'feature_canonical_V',[]);
end

function validate_landing_source(item,cfg,source_file)
d=item.discovery;
if d.configuration.grid.N~=cfg.grid.N || ...
        abs(d.configuration.grid.Lx-cfg.grid.Lx)>1e-12
    error('build_stage2_landing_geometry:GridMismatch', ...
        'Landing source is not N=%d, Lx=%g: %s',cfg.grid.N,cfg.grid.Lx,source_file);
end
if ~isequal(size(d.landing_state),[cfg.grid.N,3]) || ...
        any(~isfinite(d.landing_state),'all')
    error('build_stage2_landing_geometry:InvalidLandingState', ...
        'Landing source lacks a finite full (x,u,v) state: %s',source_file);
end
expected_x=linspace(0,cfg.grid.Lx,cfg.grid.N)';
if max(abs(d.landing_state(:,1)-expected_x))>1e-12
    error('build_stage2_landing_geometry:GridCoordinates', ...
        'Landing-state coordinates disagree with the declared grid: %s',source_file);
end
if abs(d.landing_time-d.forcing_end)>1e-9
    error('build_stage2_landing_geometry:LandingTime', ...
        'Landing state was not saved at forcing return: %s',source_file);
end
[time_error,index]=min(abs(d.S(:)-d.landing_time));
if time_error>1e-9
    error('build_stage2_landing_geometry:MissingLandingSnapshot', ...
        'The exact landing time is absent from the saved trajectory: %s',source_file);
end
snapshot_error=direct_state_distance(d.landing_state(:,2)', ...
    d.landing_state(:,3)',d.U(index,:),d.V(index,:),cfg.grid.Lx);
if snapshot_error>1e-12
    error('build_stage2_landing_geometry:LandingSnapshotMismatch', ...
        'Landing-state fields disagree with the saved forcing-end row: %s',source_file);
end
end

function family = normalize_family(value)
if startsWith(value,'symmetric_ramp')
    family='symmetric_ramp';
else
    family=value;
end
end

function [label,status] = landing_label(protocol_id,outcome,cross_records)
index=find(strcmp({cross_records.protocol_id},protocol_id),1);
if ~isempty(index)
    item=cross_records(index);
    if strcmp(item.class_id,'period_5p0585') && ...
            strcmp(item.classification,'converges_to_representative_orbit_set')
        label='R5'; status='validated_outcome_preparation'; return;
    elseif strcmp(item.class_id,'period_6p5837') && ...
            strcmp(item.classification,'converges_to_representative_orbit_set')
        label='R6'; status='validated_outcome_preparation'; return;
    elseif strcmp(item.class_id,'period_6p5837')
        label='R6_candidate_same_period'; status='provisional'; return;
    end
end
switch outcome
    case 'A_stationary_neighborhood'
        label='A'; status='validated_reference_neighborhood';
    case 'B_periodic_reflection_neighborhood'
        label='B'; status='validated_reference_neighborhood';
    case 'C_periodic_direct_neighborhood'
        label='C'; status='validated_reference_neighborhood';
    otherwise
        label='unresolved'; status='unresolved';
end
end

function role = bracket_role(record)
role='';
if abs(record.bmax-11.14)>1e-10 || record.Thold~=80 || record.Tup~=record.Tdown
    return;
end
switch record.Tup
    case 16, role='R5_B_left_ramp16';
    case 17, role='R5_B_right_ramp17';
    case 40, role='B_A_left_ramp40';
    case 45, role='B_A_right_ramp45';
end
end

function [records,states,paths,membership] = orient_landing_states(records,states,cfg)
orientation=brusselator_reflection_canonicalize(states.raw_U,states.raw_V,cfg.grid.Lx);
reflected_orientation=brusselator_reflection_canonicalize( ...
    states.reflected_U,states.reflected_V,cfg.grid.Lx);
states.feature_canonical_U=orientation.canonical_U;
states.feature_canonical_V=orientation.canonical_V;
for k=1:numel(records)
    records(k).feature_uses_reflection=orientation.uses_reflection(k);
    records(k).feature_index=orientation.selected_feature_index(k);
    records(k).feature_value=orientation.selected_feature_value(k);
    records(k).reflection_invariance_error=direct_state_distance( ...
        orientation.canonical_U(k,:),orientation.canonical_V(k,:), ...
        reflected_orientation.canonical_U(k,:), ...
        reflected_orientation.canonical_V(k,:),cfg.grid.Lx);
end
if max([records.reflection_invariance_error]) > ...
        cfg.stage2.reflection.invariance_test_tolerance
    error('build_stage2_landing_geometry:ReflectionInvariance', ...
        'Feature canonicalization failed the raw/reflected invariance gate.');
end
path_specs=family_path_specs(records);
paths=cell(numel(path_specs),1);
member_template=struct('path_name','','state_index',NaN,'protocol_id','', ...
    'order_value',NaN,'order_rank',NaN,'continuity_uses_reflection',false, ...
    'raw_step_distance',NaN,'feature_step_distance',NaN, ...
    'continuity_step_distance',NaN);
membership=member_template([]);
for p=1:numel(path_specs)
    indices=path_specs(p).indices; values=path_specs(p).values;
    [values,order]=sort(values); indices=indices(order);
    continuity_U=zeros(numel(indices),cfg.grid.N); continuity_V=continuity_U;
    flags=false(numel(indices),1); raw_steps=nan(numel(indices),1);
    feature_steps=raw_steps; continuity_steps=raw_steps;
    first=indices(1); continuity_U(1,:)=states.feature_canonical_U(first,:);
    continuity_V(1,:)=states.feature_canonical_V(first,:);
    flags(1)=records(first).feature_uses_reflection;
    for k=2:numel(indices)
        current=indices(k); previous=indices(k-1);
        raw_steps(k)=direct_state_distance(states.raw_U(current,:),states.raw_V(current,:), ...
            states.raw_U(previous,:),states.raw_V(previous,:),cfg.grid.Lx);
        feature_steps(k)=direct_state_distance(states.feature_canonical_U(current,:), ...
            states.feature_canonical_V(current,:),states.feature_canonical_U(previous,:), ...
            states.feature_canonical_V(previous,:),cfg.grid.Lx);
        direct=direct_state_distance(states.raw_U(current,:),states.raw_V(current,:), ...
            continuity_U(k-1,:),continuity_V(k-1,:),cfg.grid.Lx);
        reflected=direct_state_distance(states.reflected_U(current,:), ...
            states.reflected_V(current,:),continuity_U(k-1,:),continuity_V(k-1,:),cfg.grid.Lx);
        if reflected<direct
            continuity_U(k,:)=states.reflected_U(current,:);
            continuity_V(k,:)=states.reflected_V(current,:); flags(k)=true;
            continuity_steps(k)=reflected;
        else
            continuity_U(k,:)=states.raw_U(current,:);
            continuity_V(k,:)=states.raw_V(current,:); flags(k)=false;
            continuity_steps(k)=direct;
        end
    end
    paths{p}=struct('name',path_specs(p).name,'state_indices',indices, ...
        'order_values',values,'raw_U',states.raw_U(indices,:), ...
        'raw_V',states.raw_V(indices,:),'feature_U',states.feature_canonical_U(indices,:), ...
        'feature_V',states.feature_canonical_V(indices,:), ...
        'continuity_U',continuity_U,'continuity_V',continuity_V, ...
        'continuity_uses_reflection',flags,'raw_step_distance',raw_steps, ...
        'feature_step_distance',feature_steps,'continuity_step_distance',continuity_steps);
    for k=1:numel(indices)
        item=member_template; item.path_name=path_specs(p).name;
        item.state_index=indices(k); item.protocol_id=records(indices(k)).protocol_id;
        item.order_value=values(k); item.order_rank=k;
        item.continuity_uses_reflection=flags(k); item.raw_step_distance=raw_steps(k);
        item.feature_step_distance=feature_steps(k);
        item.continuity_step_distance=continuity_steps(k);
        membership(end+1,1)=item; %#ok<AGROW>
    end
end
end

function specs=family_path_specs(records)
bmax=[records.bmax]; up=[records.Tup]; hold=[records.Thold]; down=[records.Tdown];
specs(1)=struct('name','amplitude','indices',find(up==40 & hold==80 & down==40), ...
    'values',bmax(up==40 & hold==80 & down==40));
specs(2)=struct('name','hold','indices',find(abs(bmax-11.14)<1e-10 & up==40 & down==40), ...
    'values',hold(abs(bmax-11.14)<1e-10 & up==40 & down==40));
mask=abs(bmax-11.14)<1e-10 & hold==80 & up==down;
specs(3)=struct('name','symmetric_ramp','indices',find(mask),'values',up(mask));
mask=strcmp({records.family},'asymmetric_ramp');
specs(4)=struct('name','asymmetric_ramp','indices',find(mask),'values',up(mask));
end

function [records,distances] = add_physical_distances(records,states,outcomes,cfg)
distances=zeros(numel(records),numel(outcomes));
for i=1:numel(records)
    for j=1:numel(outcomes)
        if strcmp(outcomes(j).type,'stationary')
            ref_U=outcomes(j).orbit_U; ref_V=outcomes(j).orbit_V;
        else
            ref_U=outcomes(j).dense_template.U;
            ref_V=outcomes(j).dense_template.V;
        end
        item=brusselator_physical_orbit_distance(states.raw_U(i,:),states.raw_V(i,:), ...
            ref_U,ref_V,cfg.grid.Lx);
        distances(i,j)=item.distance;
    end
    records(i).distance_A=distances(i,1); records(i).distance_B=distances(i,2);
    records(i).distance_C=distances(i,3); records(i).distance_R5=distances(i,4);
    records(i).distance_R6=distances(i,5);
end
end

function [U,V,weights,labels]=build_weighted_training_set(outcomes,cfg)
U=outcomes(1).orbit_U; V=outcomes(1).orbit_V; weights=0.2; labels={'A'};
for k=2:5
    count=size(outcomes(k).orbit_U,1);
    U=[U;outcomes(k).orbit_U]; %#ok<AGROW>
    V=[V;outcomes(k).orbit_V]; %#ok<AGROW>
    weights=[weights;repmat(0.2/count,count,1)]; %#ok<AGROW>
    labels=[labels;repmat({outcomes(k).label},count,1)]; %#ok<AGROW>
end
if abs(sum(weights)-1)>1e-12
    error('build_stage2_landing_geometry:WeightSum','POD weights do not sum to one.');
end
expected=cfg.stage2.pod.class_aggregate_weights;
for label={'A','B','C','R5','R6'}
    aggregate=sum(weights(strcmp(labels,label{1})));
    if abs(aggregate-expected.(label{1}))>1e-12
        error('build_stage2_landing_geometry:ClassWeight', ...
            'Aggregate weight for %s is incorrect.',label{1});
    end
end
end

function [projections,reconstruction] = project_all_states(outcomes,states,paths, ...
        training_U,training_V,training_labels,primary,unscaled,cfg)
counts=cfg.stage2.pod.diagnostic_component_counts;
[~,train_primary_errors]=brusselator_project_pod(training_U,training_V,primary,counts);
[~,train_unscaled_errors]=brusselator_project_pod(training_U,training_V,unscaled,counts);
variants={'raw','reflected','feature_canonical'};
projections=struct();
for k=1:numel(variants)
    switch variants{k}
        case 'raw', U=states.raw_U; V=states.raw_V;
        case 'reflected', U=states.reflected_U; V=states.reflected_V;
        otherwise, U=states.feature_canonical_U; V=states.feature_canonical_V;
    end
    [p,e1]=brusselator_project_pod(U,V,primary,counts);
    [q,e2]=brusselator_project_pod(U,V,unscaled,counts);
    projections.landing.(variants{k})=struct('primary',p,'unscaled',q, ...
        'primary_errors',e1,'unscaled_errors',e2);
end
projections.outcomes=cell(numel(outcomes),1);
for k=1:numel(outcomes)
    [p,e1]=brusselator_project_pod(outcomes(k).orbit_U,outcomes(k).orbit_V,primary,counts);
    [q,e2]=brusselator_project_pod(outcomes(k).orbit_U,outcomes(k).orbit_V,unscaled,counts);
    [pf,~]=brusselator_project_pod(outcomes(k).phase_fixed.canonical_U, ...
        outcomes(k).phase_fixed.canonical_V,primary,counts);
    projections.outcomes{k}=struct('label',outcomes(k).label, ...
        'primary',p,'unscaled',q,'phase_fixed_primary',pf, ...
        'primary_errors',e1,'unscaled_errors',e2);
end
projections.family_paths=cell(numel(paths),1);
for k=1:numel(paths)
    [raw,~]=brusselator_project_pod(paths{k}.raw_U,paths{k}.raw_V,primary,counts);
    [feature,~]=brusselator_project_pod(paths{k}.feature_U,paths{k}.feature_V,primary,counts);
    [continuity,~]=brusselator_project_pod(paths{k}.continuity_U,paths{k}.continuity_V,primary,counts);
    projections.family_paths{k}=struct('name',paths{k}.name,'raw',raw, ...
        'feature',feature,'continuity',continuity);
end
reconstruction=reconstruction_summary(training_labels,train_primary_errors, ...
    train_unscaled_errors,repmat({'training'},numel(training_labels),1),counts);
landing_labels=repmat({'landing'},size(states.raw_U,1),1);
landing_groups=repmat({'all_landings'},size(states.raw_U,1),1);
landing_summary=reconstruction_summary(landing_groups, ...
    projections.landing.feature_canonical.primary_errors, ...
    projections.landing.feature_canonical.unscaled_errors,landing_labels,counts);
reconstruction=[reconstruction;landing_summary];
end

function records=reconstruction_summary(groups,primary,unscaled,dataset,counts)
unique_groups=unique(groups,'stable');
template=struct('dataset','','group','','component_count',NaN, ...
    'primary_median_relative_error',NaN,'primary_maximum_relative_error',NaN, ...
    'unscaled_median_relative_error',NaN,'unscaled_maximum_relative_error',NaN);
records=template([]);
for g=1:numel(unique_groups)
    selected=strcmp(groups,unique_groups{g});
    for k=1:numel(counts)
        item=template; item.dataset=dataset{find(selected,1)}; item.group=unique_groups{g};
        item.component_count=counts(k);
        item.primary_median_relative_error=median(primary(selected,k));
        item.primary_maximum_relative_error=max(primary(selected,k));
        item.unscaled_median_relative_error=median(unscaled(selected,k));
        item.unscaled_maximum_relative_error=max(unscaled(selected,k));
        records(end+1,1)=item; %#ok<AGROW>
    end
end
end

function records=compare_scalings(primary,unscaled,primary_scores,unscaled_scores,cfg)
counts=[2,3,5];
template=struct('component_count',NaN,'maximum_principal_angle_degrees',NaN, ...
    'median_principal_angle_degrees',NaN,'landing_pairwise_distance_correlation',NaN);
records=repmat(template,numel(counts),1);
N=cfg.grid.N;
for k=1:numel(counts)
    count=counts(k);
    P=primary.modes(:,1:count);
    P(1:N,:)=primary.scale_U*P(1:N,:);
    P(N+1:end,:)=primary.scale_V*P(N+1:end,:);
    [P,~]=qr(P,0); Q=unscaled.modes(:,1:count);
    singular=svd(P'*Q); angles=acosd(min(max(singular,-1),1));
    d1=pairwise_vector(primary_scores(:,1:count));
    d2=pairwise_vector(unscaled_scores(:,1:count));
    correlation=corrcoef(d1,d2);
    records(k)=struct('component_count',count, ...
        'maximum_principal_angle_degrees',max(angles), ...
        'median_principal_angle_degrees',median(angles), ...
        'landing_pairwise_distance_correlation',correlation(1,2));
end
end

function vector=pairwise_vector(scores)
n=size(scores,1); vector=zeros(n*(n-1)/2,1); count=0;
for i=1:n
    for j=i+1:n
        count=count+1; vector(count)=norm(scores(i,:)-scores(j,:));
    end
end
end

function [landing_aware,diagnostic]=build_landing_aware_diagnostic( ...
        outcomes,states,training_U,training_V,training_weights,projections,primary,cfg)
counts=cfg.stage2.pod.diagnostic_component_counts;
[~,training_errors]=brusselator_project_pod(training_U,training_V,primary,counts);
landing_errors=projections.landing.feature_canonical.primary_errors;
component_count=counts(end);
training_max=max(training_errors(:,end));
landing_max=max(landing_errors(:,end));
error_ratio=landing_max/max(training_max,eps);
triggered=error_ratio>=cfg.stage2.pod.landing_aware_trigger_maximum_error_ratio;
diagnostic=struct('triggered',triggered,'trigger_component_count',component_count, ...
    'orbit_library_maximum_relative_error',training_max, ...
    'landing_maximum_relative_error',landing_max, ...
    'maximum_error_ratio',error_ratio, ...
    'trigger_ratio',cfg.stage2.pod.landing_aware_trigger_maximum_error_ratio, ...
    'interpretation',['diagnostic basis only; the balanced outcome-only POD ', ...
        'remains the primary coordinate system']);
landing_aware=struct('triggered',triggered,'model',[], ...
    'sample_weights',[],'sample_labels',{{}},'landing_scores',[], ...
    'landing_errors',[],'outcome_scores',{{}},'outcome_errors',{{}});
if ~triggered
    return;
end
n_landing=size(states.feature_canonical_U,1);
outcome_weight=cfg.stage2.pod.landing_aware_outcome_aggregate_weight;
landing_weight=cfg.stage2.pod.landing_aware_landing_aggregate_weight;
weights=[outcome_weight*training_weights; ...
    repmat(landing_weight/n_landing,n_landing,1)];
labels=[repmat({'balanced_outcome_library'},numel(training_weights),1); ...
    repmat({'forcing_landing'},n_landing,1)];
combined_U=[training_U;states.feature_canonical_U];
combined_V=[training_V;states.feature_canonical_V];
model=brusselator_weighted_pod(combined_U,combined_V,weights, ...
    cfg.grid.Lx,true,cfg.stage2.pod.maximum_components);
[landing_scores,landing_aware_errors]=brusselator_project_pod( ...
    states.feature_canonical_U,states.feature_canonical_V,model,counts);
outcome_scores=cell(numel(outcomes),1); outcome_errors=outcome_scores;
for k=1:numel(outcomes)
    [outcome_scores{k},outcome_errors{k}]=brusselator_project_pod( ...
        outcomes(k).orbit_U,outcomes(k).orbit_V,model,counts);
end
landing_aware=struct('triggered',true,'model',model, ...
    'sample_weights',weights,'sample_labels',{labels}, ...
    'landing_scores',landing_scores,'landing_errors',landing_aware_errors, ...
    'outcome_scores',{outcome_scores},'outcome_errors',{outcome_errors});
end

function write_landing_aware_tables(outdir,landing_aware,diagnostic,cfg)
writetable(struct2table(diagnostic,'AsArray',true), ...
    fullfile(outdir,'landing_aware_pod_diagnostic.csv'));
if ~landing_aware.triggered
    return;
end
weights=table((1:numel(landing_aware.sample_weights))', ...
    landing_aware.sample_labels,landing_aware.sample_weights, ...
    'VariableNames',{'sample_index','sample_group','sample_weight'});
writetable(weights,fullfile(outdir,'landing_aware_pod_training_weights.csv'));
model=landing_aware.model;
writetable(table((1:numel(model.explained_percent))',model.explained_percent, ...
    'VariableNames',{'component','explained_percent'}), ...
    fullfile(outdir,'landing_aware_pod_explained_variance.csv'));
counts=cfg.stage2.pod.diagnostic_component_counts;
writetable(table(counts(:),median(landing_aware.landing_errors,1)', ...
    max(landing_aware.landing_errors,[],1)', ...
    'VariableNames',{'component_count','median_landing_relative_error', ...
        'maximum_landing_relative_error'}), ...
    fullfile(outdir,'landing_aware_pod_reconstruction_summary.csv'));
end

function write_pod_tables(outdir,primary,unscaled,labels,weights)
weight_table=table((1:numel(weights))',labels,weights, ...
    'VariableNames',{'sample_index','outcome_label','sample_weight'});
writetable(weight_table,fullfile(outdir,'pod_training_weights.csv'));
n=max(numel(primary.explained_percent),numel(unscaled.explained_percent));
component=(1:n)'; p=nan(n,1); u=nan(n,1);
p(1:numel(primary.explained_percent))=primary.explained_percent;
u(1:numel(unscaled.explained_percent))=unscaled.explained_percent;
writetable(table(component,p,u,'VariableNames', ...
    {'component','primary_explained_percent','unscaled_explained_percent'}), ...
    fullfile(outdir,'pod_explained_variance.csv'));
end

function make_figures(outdir,outcomes,records,distances,paths,proj,primary,unscaled,reconstruction,landing_aware,landing_aware_diagnostic)
labels={records.final_label}; colors=label_colors(labels);
fig=figure('Color','w','Position',[100 100 1000 650]); ax=axes(fig); hold(ax,'on');
for k=1:numel(outcomes)
    s=proj.outcomes{k}.phase_fixed_primary;
    scatter(ax,s(:,1),s(:,2),90,outcome_color(outcomes(k).label),'filled', ...
        'DisplayName',[outcomes(k).label,' phase-fixed']);
end
plot_landing_groups(ax,records,proj.landing.feature_canonical.primary,45);
mark_brackets(ax,records,proj.landing.feature_canonical.primary);
xlabel(ax,'PC1'); ylabel(ax,'PC2'); title(ax,'Phase-fixed outcomes and canonical landing states');
style_axes(ax); style_legend(legend(ax,'Location','eastoutside'));
exportgraphics(fig,fullfile(outdir,'phase_fixed_outcomes_and_landings_PC1_PC2.png'),'Resolution',300); close(fig);

fig=figure('Color','w','Position',[100 100 1000 650]); ax=axes(fig); hold(ax,'on');
for k=1:numel(outcomes)
    s=proj.outcomes{k}.primary;
    if size(s,1)==1
        scatter(ax,s(:,1),s(:,2),100,outcome_color(outcomes(k).label),'filled', ...
            'DisplayName',outcomes(k).label);
    else
        plot(ax,[s(:,1);s(1,1)],[s(:,2);s(1,2)],'LineWidth',2, ...
            'Color',outcome_color(outcomes(k).label),'DisplayName',outcomes(k).label);
    end
end
plot_landing_groups(ax,records,proj.landing.feature_canonical.primary,35);
mark_brackets(ax,records,proj.landing.feature_canonical.primary);
xlabel(ax,'PC1'); ylabel(ax,'PC2'); title(ax,'Full periodic orbits and forcing landing states');
style_axes(ax); style_legend(legend(ax,'Location','eastoutside'));
exportgraphics(fig,fullfile(outdir,'full_orbit_loops_and_landings_PC1_PC2.png'),'Resolution',300); close(fig);

fig=figure('Color','w','Position',[100 100 1000 700]); ax=axes(fig); hold(ax,'on');
for k=1:numel(outcomes)
    s=proj.outcomes{k}.primary;
    if size(s,1)==1
        scatter3(ax,s(:,1),s(:,2),s(:,3),100,outcome_color(outcomes(k).label),'filled');
    else
        plot3(ax,[s(:,1);s(1,1)],[s(:,2);s(1,2)],[s(:,3);s(1,3)], ...
            'LineWidth',2,'Color',outcome_color(outcomes(k).label));
    end
end
scatter3(ax,proj.landing.feature_canonical.primary(:,1), ...
    proj.landing.feature_canonical.primary(:,2), ...
    proj.landing.feature_canonical.primary(:,3),35,colors,'filled','MarkerEdgeColor','k');
xlabel(ax,'PC1');ylabel(ax,'PC2');zlabel(ax,'PC3');title(ax,'Full-orbit Stage 2 geometry');
style_axes(ax); view(ax,35,25);
exportgraphics(fig,fullfile(outdir,'full_orbit_loops_PC1_PC2_PC3.png'),'Resolution',300); close(fig);

fig=figure('Color','w','Position',[100 100 1200 420]); tiledlayout(fig,1,3,'TileSpacing','compact');
variants={'raw','feature','continuity'}; titles={'Raw orientation','Feature canonical','Continuity oriented'};
for q=1:3
    ax=nexttile; hold(ax,'on');
    for p=1:numel(paths)
        s=proj.family_paths{p}.(variants{q});
        plot(ax,s(:,1),s(:,2),'-o','LineWidth',1.2,'MarkerSize',4, ...
            'DisplayName',paths{p}.name);
    end
    title(ax,titles{q});xlabel(ax,'PC1');ylabel(ax,'PC2');style_axes(ax);
    if q==3,style_legend(legend(ax,'Location','eastoutside'));end
end
exportgraphics(fig,fullfile(outdir,'reflection_continuity_diagnostic.png'),'Resolution',300); close(fig);

fig=figure('Color','w','Position',[100 100 1100 480]); tiledlayout(fig,1,2,'TileSpacing','compact');
ax=nexttile; scatter(ax,proj.landing.feature_canonical.primary(:,1), ...
    proj.landing.feature_canonical.primary(:,2),45,colors,'filled','MarkerEdgeColor','k');
title(ax,'Primary POD: separate u/v RMS scales');xlabel(ax,'PC1');ylabel(ax,'PC2');style_axes(ax);
ax=nexttile; scatter(ax,proj.landing.feature_canonical.unscaled(:,1), ...
    proj.landing.feature_canonical.unscaled(:,2),45,colors,'filled','MarkerEdgeColor','k');
title(ax,'Sensitivity POD: original u/v units');xlabel(ax,'PC1');ylabel(ax,'PC2');style_axes(ax);
exportgraphics(fig,fullfile(outdir,'pod_scaling_sensitivity.png'),'Resolution',300);close(fig);

fig=figure('Color','w','Position',[100 100 1000 450]); tiledlayout(fig,1,2,'TileSpacing','compact');
ax=nexttile; plot(ax,cumsum(primary.explained_percent),'-o','LineWidth',1.4);hold(ax,'on');
plot(ax,cumsum(unscaled.explained_percent),'-s','LineWidth',1.4);xlabel(ax,'Component');
ylabel(ax,'Cumulative explained variance (%)');style_axes(ax);
style_legend(legend(ax,{'primary','unscaled'},'Location','southeast'));
ax=nexttile; rows=strcmp({reconstruction.dataset},'landing');
items=reconstruction(rows); counts=[items.component_count];
plot(ax,counts,[items.primary_median_relative_error],'-o','LineWidth',1.4);hold(ax,'on');
plot(ax,counts,[items.unscaled_median_relative_error],'-s','LineWidth',1.4);
xlabel(ax,'Components');ylabel(ax,'Median landing reconstruction error');style_axes(ax);
style_legend(legend(ax,{'primary','unscaled'}));
exportgraphics(fig,fullfile(outdir,'POD_explained_variance_and_reconstruction_error.png'),'Resolution',300);close(fig);

fig=figure('Color','w','Position',[100 100 900 800]);ax=axes(fig);imagesc(ax,distances);colorbar(ax);
xticks(ax,1:5);xticklabels(ax,{'A','B','C','R5','R6'});yticks(ax,1:numel(records));
yticklabels(ax,{records.protocol_id});set(ax,'FontSize',7,'TickLabelInterpreter','none');
title(ax,'Physical phase/reflection-aware landing distances');
xlabel(ax,'Validated outcome');ylabel(ax,'Forcing protocol');style_axes(ax);
exportgraphics(fig,fullfile(outdir,'landing_distance_to_each_validated_orbit.png'),'Resolution',300);close(fig);

if landing_aware.triggered
    fig=figure('Color','w','Position',[100 100 1050 650]);ax=axes(fig);hold(ax,'on');
    for k=1:numel(outcomes)
        s=landing_aware.outcome_scores{k};
        if size(s,1)==1
            scatter(ax,s(:,1),s(:,2),100,outcome_color(outcomes(k).label),'filled', ...
                'DisplayName',['reference ',outcomes(k).label]);
        else
            plot(ax,[s(:,1);s(1,1)],[s(:,2);s(1,2)],'LineWidth',2, ...
                'Color',outcome_color(outcomes(k).label), ...
                'DisplayName',['reference ',outcomes(k).label]);
        end
    end
    plot_landing_groups(ax,records,landing_aware.landing_scores,35);
    xlabel(ax,'diagnostic PC1');ylabel(ax,'diagnostic PC2');
    title(ax,sprintf('Landing-aware POD (diagnostic only; trigger ratio %.2f)', ...
        landing_aware_diagnostic.maximum_error_ratio));
    style_axes(ax);style_legend(legend(ax,'Location','eastoutside'));
    exportgraphics(fig,fullfile(outdir,'landing_aware_POD_diagnostic.png'),'Resolution',300);
    close(fig);
end
end

function colors=label_colors(labels)
colors=zeros(numel(labels),3);
for k=1:numel(labels),colors(k,:)=outcome_color(labels{k});end
end

function color=outcome_color(label)
switch label
    case 'A',color=[0.10 0.35 0.80];
    case 'B',color=[0.90 0.35 0.10];
    case 'C',color=[0.85 0.65 0.05];
    case 'R5',color=[0.50 0.20 0.75];
    case 'R6',color=[0.10 0.60 0.30];
    case 'R6_candidate_same_period',color=[0.50 0.80 0.55];
    otherwise,color=[0.45 0.45 0.45];
end
end

function mark_brackets(ax,records,scores)
indices=find(~cellfun(@isempty,{records.bracket_role}));
scatter(ax,scores(indices,1),scores(indices,2),120,'s','MarkerEdgeColor','k', ...
    'MarkerFaceColor','none','LineWidth',2,'DisplayName','transition brackets');
for k=indices(:)'
    text(ax,scores(k,1),scores(k,2),sprintf(' T_r=%.0f',records(k).Tup), ...
        'FontSize',7,'Color','k','Interpreter','none','VerticalAlignment','bottom');
end
end

function plot_landing_groups(ax,records,scores,marker_size)
labels=unique({records.final_label},'stable');
for k=1:numel(labels)
    selected=strcmp({records.final_label},labels{k});
    marker='o';
    if any(strcmp({records(selected).label_status},'provisional'))
        marker='^';
    end
    scatter(ax,scores(selected,1),scores(selected,2),marker_size, ...
        outcome_color(labels{k}),'filled','Marker',marker,'MarkerEdgeColor','k', ...
        'DisplayName',['landing ',labels{k}]);
end
end

function style_axes(ax)
set(ax,'Color','w','XColor','k','YColor','k','ZColor','k', ...
    'GridColor',[0.75 0.75 0.75]);grid(ax,'on');box(ax,'on');
ax.Title.Color='k';ax.XLabel.Color='k';ax.YLabel.Color='k';
if ~isempty(ax.ZLabel),ax.ZLabel.Color='k';end
end

function style_legend(leg)
set(leg,'Color','w','TextColor','k','EdgeColor',[0.25 0.25 0.25], ...
    'Interpreter','none');
end

function distance=direct_state_distance(u,v,reference_u,reference_v,Lx)
N=numel(u);dx=Lx/(N-1);weights=dx*ones(1,N);weights([1,end])=0.5*dx;
numerator=sum(((u-reference_u).^2+(v-reference_v).^2).*weights);
denominator=sum((reference_u.^2+reference_v.^2).*weights);
distance=sqrt(numerator/max(denominator,eps));
end
