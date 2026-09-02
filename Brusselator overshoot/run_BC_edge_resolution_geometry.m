function results = run_BC_edge_resolution_geometry()
%RUN_BC_EDGE_RESOLUTION_GEOMETRY Continue exact B/C and retest edge branches.
%   Newton-converges B and C at N=400, continues both exact periodic orbits
%   to N=200 and N=600, verifies their leading Floquet spectra, then follows
%   both signs of the consistently oriented unstable edge direction at every
%   grid. Classification uses exact grid-matched B/C orbit libraries and
%   retains unresolved outcomes rather than forcing a binary label.

cfg = brusselator_BC_edge_resolution_config();
outdir = fullfile(cfg.output.root,cfg.experiment_name);
checkpoint = fullfile(outdir,'BC_edge_resolution_checkpoint.mat');
final_file = fullfile(outdir,'BC_edge_resolution_geometry.mat');
if isfile(final_file)
    error('run_BC_edge_resolution_geometry:OutputExists', ...
        'Refusing to overwrite completed result: %s',final_file);
end
if ~isfile(cfg.BC_resolution.stage2_file) || ...
        ~isfile(cfg.BC_resolution.edge_resolution_file)
    error('run_BC_edge_resolution_geometry:MissingSource', ...
        'The Stage-2 library or exact edge-resolution study is missing.');
end
stage2_loaded = load(cfg.BC_resolution.stage2_file,'results');
stage2 = stage2_loaded.results;
edge_loaded = load(cfg.BC_resolution.edge_resolution_file,'results');
edge_resolution = edge_loaded.results;

orbit_specs = make_orbit_specs(cfg);
branch_specs = make_branch_specs(cfg);
orbit_template = empty_orbit_record();
branch_template = empty_branch_record();
if isfile(checkpoint)
    saved = load(checkpoint,'orbit_records','orbit_files','next_orbit_case', ...
        'branch_records','branch_files','next_branch_case');
    orbit_records=saved.orbit_records; orbit_files=saved.orbit_files;
    next_orbit_case=saved.next_orbit_case;
    branch_records=saved.branch_records; branch_files=saved.branch_files;
    next_branch_case=saved.next_branch_case;
    fprintf('Resuming matched B--edge--C study.\n');
else
    if isfolder(outdir)
        error('run_BC_edge_resolution_geometry:IncompleteOutput', ...
            'Existing output lacks a checkpoint: %s',outdir);
    end
    mkdir(outdir); mkdir(fullfile(outdir,'seeds')); mkdir(fullfile(outdir,'branches'));
    orbit_records=orbit_template([]); orbit_files=struct( ...
        'label',{},'N',{},'orbit_file',{},'floquet_file',{});
    next_orbit_case=1; branch_records=branch_template([]); branch_files={};
    next_branch_case=1;
    brusselator_write_metadata(outdir,cfg, ...
        brusselator_run_metadata(cfg,mfilename,NaN));
    save_checkpoint();
end
timer=tic;

for case_index=next_orbit_case:numel(orbit_specs)
    spec=orbit_specs(case_index);
    fprintf('\nExact attractor continuation: %s, N=%d.\n',spec.label,spec.N);
    if spec.N==cfg.BC_resolution.source_resolution
        seed_file=make_stage2_seed(stage2,spec.label,spec.N,cfg,outdir);
        source_kind='stage2_phase_fixed_seed';
    else
        source_file=attractor_orbit_file(outdir,spec.label, ...
            cfg.BC_resolution.source_resolution);
        if ~isfile(source_file)
            error('run_BC_edge_resolution_geometry:MissingN400Attractor', ...
                'The exact N=400 %s orbit must be computed first.',spec.label);
        end
        source_loaded=load(source_file,'results');
        seed_file=make_regridded_seed(source_loaded.results,spec.label,spec.N,cfg,outdir);
        source_kind='regridded_exact_N400_seed';
    end

    newton_cfg=brusselator_periodic_edge_newton_config();
    newton_cfg.grid.N=spec.N; newton_cfg.output.root=outdir;
    newton_cfg.experiment_name=sprintf('%s_N%04d_newton',spec.label,spec.N);
    newton_cfg.newton.seed_file=seed_file;
    orbit_file=fullfile(outdir,newton_cfg.experiment_name, ...
        'periodic_edge_orbit_newton.mat');
    if isfile(orbit_file)
        loaded=load(orbit_file,'results'); orbit=loaded.results;
    else
        orbit=run_periodic_edge_orbit_newton(newton_cfg);
    end
    if ~strcmp(orbit.status,'newton_converged')
        error('run_BC_edge_resolution_geometry:AttractorNewtonFailure', ...
            '%s at N=%d did not Newton converge.',spec.label,spec.N);
    end

    floquet_cfg=brusselator_periodic_edge_floquet_config();
    floquet_cfg.grid.N=spec.N; floquet_cfg.output.root=outdir;
    floquet_cfg.experiment_name=sprintf('%s_N%04d_floquet',spec.label,spec.N);
    floquet_cfg.floquet.orbit_file=orbit_file;
    floquet_cfg.floquet.predicted_transverse_multiplier=NaN;
    floquet_file=fullfile(outdir,floquet_cfg.experiment_name, ...
        'periodic_edge_orbit_floquet.mat');
    if isfile(floquet_file)
        loaded=load(floquet_file,'results'); floquet=loaded.results;
    else
        floquet=run_periodic_edge_orbit_floquet(floquet_cfg);
    end
    record=summarize_attractor(spec,source_kind,orbit,floquet, ...
        orbit_file,floquet_file,cfg);
    orbit_records(end+1,1)=record; %#ok<AGROW>
    orbit_files(end+1,1)=struct('label',spec.label,'N',spec.N, ...
        'orbit_file',orbit_file,'floquet_file',floquet_file); %#ok<AGROW>
    next_orbit_case=case_index+1;
    writetable(struct2table(orbit_records,'AsArray',true), ...
        fullfile(outdir,'attractor_progress.csv'));
    save_checkpoint();
end

orientation_reference=edge_direction_reference(edge_resolution,cfg);
for case_index=next_branch_case:numel(branch_specs)
    spec=branch_specs(case_index);
    fprintf('\nMatched edge branch: N=%d, sign=%+d.\n',spec.N,spec.sign);
    edge=load_edge_case(edge_resolution,spec.N);
    B=load_attractor(outdir,'B',spec.N);
    C=load_attractor(outdir,'C',spec.N);
    library=exact_BC_library(B,C,cfg.BC_resolution.exact_library_phase_samples);
    [direction_y,orientation_factor,orientation_inner_product]= ...
        oriented_edge_direction(edge,orientation_reference,cfg);
    state=edge.state; N=spec.N; w=spatial_weights(state(:,1));
    state_norm=norm(sqrt([w;w]).*[state(:,2);state(:,3)]);
    delta=spec.sign*cfg.BC_resolution.branch_relative_amplitude* ...
        state_norm*direction_y;
    initial_state=state;
    initial_state(:,2)=initial_state(:,2)+delta(1:N);
    initial_state(:,3)=initial_state(:,3)+delta(N+1:end);
    case_id=sprintf('N%04d_sign_%s',N,sign_token(spec.sign));
    raw_file=fullfile(outdir,'branches',[case_id,'.mat']);
    partial_file=fullfile(outdir,'branches',[case_id,'_partial.mat']);
    branch_cfg=cfg; branch_cfg.grid.N=N;
    if isfile(raw_file)
        loaded=load(raw_file,'edge_case'); edge_case=loaded.edge_case;
    else
        edge_case=brusselator_evolve_and_classify_state( ...
            case_id,initial_state,library,branch_cfg,partial_file);
        save(raw_file,'edge_case','library','-v7.3');
    end
    record=branch_template;
    record.N=N; record.sign=spec.sign;
    record.relative_amplitude=cfg.BC_resolution.branch_relative_amplitude;
    record.orientation_factor=orientation_factor;
    record.orientation_inner_product=orientation_inner_product;
    record.initial_relative_distance=state_distance(initial_state,state,cfg.grid.Lx);
    record.outcome=edge_case.classification.outcome;
    record.tested_duration=edge_case.duration;
    record.decision_time=edge_case.decision_time;
    record.late_median_distance_B=edge_case.classification.late_median_distances(1);
    record.late_median_distance_C=edge_case.classification.late_median_distances(2);
    record.late_maximum_distance_B=edge_case.classification.late_maximum_distances(1);
    record.late_maximum_distance_C=edge_case.classification.late_maximum_distances(2);
    record.raw_file=raw_file;
    branch_records(end+1,1)=record; branch_files{end+1,1}=raw_file; %#ok<AGROW>
    next_branch_case=case_index+1;
    writetable(struct2table(branch_records,'AsArray',true), ...
        fullfile(outdir,'branch_progress.csv'));
    save_checkpoint();
end

[orbit_records,orbit_order]=sort_orbit_records(orbit_records);
orbit_files=orbit_files(orbit_order);
[branch_records,branch_order]=sort_branch_records(branch_records);
branch_files=branch_files(branch_order);
comparison_records=attractor_comparisons(orbit_files,cfg);
resolution_records=make_resolution_records(orbit_records,comparison_records,cfg);
branch_resolution_records=summarize_branches(branch_records,cfg);
summary=make_summary(orbit_records,branch_resolution_records,cfg);
writetable(struct2table(orbit_records,'AsArray',true), ...
    fullfile(outdir,'exact_attractor_continuation.csv'));
writetable(struct2table(comparison_records,'AsArray',true), ...
    fullfile(outdir,'attractor_orbit_comparisons.csv'));
writetable(struct2table(resolution_records,'AsArray',true), ...
    fullfile(outdir,'attractor_resolution_diagnostics.csv'));
writetable(struct2table(branch_records,'AsArray',true), ...
    fullfile(outdir,'matched_edge_branch_outcomes.csv'));
writetable(struct2table(branch_resolution_records,'AsArray',true), ...
    fullfile(outdir,'matched_edge_branch_summary.csv'));
writetable(struct2table(summary,'AsArray',true), ...
    fullfile(outdir,'BC_edge_geometry_summary.csv'));
make_figures(outdir,orbit_records,comparison_records,branch_records,cfg);
results=struct('configuration',cfg,'orbit_records',orbit_records, ...
    'orbit_files',orbit_files,'comparison_records',comparison_records, ...
    'resolution_records',resolution_records,'branch_records',branch_records, ...
    'branch_files',{branch_files}, ...
    'branch_resolution_records',branch_resolution_records,'summary',summary);
save(final_file,'results','-v7.3');
brusselator_write_metadata(outdir,cfg, ...
    brusselator_run_metadata(cfg,mfilename,toc(timer)));
disp(struct2table(orbit_records,'AsArray',true));
disp(struct2table(comparison_records,'AsArray',true));
disp(struct2table(branch_records,'AsArray',true));
disp(struct2table(summary,'AsArray',true));
fprintf('Matched B--edge--C grid study saved in: %s\n',outdir);

    function save_checkpoint()
        save(checkpoint,'cfg','orbit_records','orbit_files','next_orbit_case', ...
            'branch_records','branch_files','next_branch_case','-v7.3');
    end
end

function specs=make_orbit_specs(cfg)
template=struct('label','','N',NaN); specs=repmat(template,6,1); index=0;
for label=cfg.BC_resolution.outcomes
    order=[cfg.BC_resolution.source_resolution, ...
        cfg.BC_resolution.resolutions( ...
        cfg.BC_resolution.resolutions~=cfg.BC_resolution.source_resolution)];
    for N=order
        index=index+1; specs(index)=struct('label',label{1},'N',N);
    end
end
end

function specs=make_branch_specs(cfg)
template=struct('N',NaN,'sign',NaN);
specs=repmat(template,numel(cfg.BC_resolution.resolutions)*2,1); index=0;
for N=cfg.BC_resolution.resolutions
    for sign_value=cfg.BC_resolution.branch_signs
        index=index+1; specs(index)=struct('N',N,'sign',sign_value);
    end
end
end

function file=make_stage2_seed(stage2,label,N,cfg,outdir)
index=find(strcmp({stage2.outcomes.label},label),1);
outcome=stage2.outcomes(index); x=linspace(0,cfg.grid.Lx,N)';
u=outcome.phase_fixed.canonical_U(1,:)';
v=outcome.phase_fixed.canonical_V(1,:)';
state=[x,u,v]; target_cfg=brusselator_periodic_edge_newton_config();
target_cfg.grid.N=N;
seed=make_seed_struct(state,outcome.period,label,N,target_cfg, ...
    'stage2_phase_fixed_candidate');
file=fullfile(outdir,'seeds',sprintf('%s_N%04d_seed.mat',label,N));
save_if_consistent(file,seed,label,N);
end

function file=make_regridded_seed(source,label,N,cfg,outdir)
x=linspace(0,cfg.grid.Lx,N)';
state=[x,interp1(source.state(:,1),source.state(:,2),x,'pchip'), ...
    interp1(source.state(:,1),source.state(:,3),x,'pchip')];
target_cfg=brusselator_periodic_edge_newton_config(); target_cfg.grid.N=N;
seed=make_seed_struct(state,source.period,label,N,target_cfg, ...
    'regridded_exact_N400_attractor');
file=fullfile(outdir,'seeds',sprintf('%s_N%04d_seed.mat',label,N));
save_if_consistent(file,seed,label,N);
end

function seed=make_seed_struct(state,period,label,N,target_cfg,status)
tangent=brusselator_frozen_rhs(state,target_cfg,target_cfg.edge_tracking.frozen_b);
seed=struct('state',state,'period',period,'event_index',NaN, ...
    'reference_state',state,'reference_tangent',tangent, ...
    'phase_condition','<X-X_ref,F_target_grid(X_ref)>_w=0', ...
    'flow_residual',NaN,'phase_condition_residual',0,'status',status, ...
    'outcome_label',label,'target_resolution',N);
end

function save_if_consistent(file,seed,label,N)
if isfile(file)
    loaded=load(file,'seed');
    if loaded.seed.target_resolution~=N || ...
            ~strcmp(loaded.seed.outcome_label,label)
        error('run_BC_edge_resolution_geometry:SeedMismatch', ...
            'Existing seed does not match %s N=%d.',label,N);
    end
else
    save(file,'seed','-v7.3');
end
end

function file=attractor_orbit_file(outdir,label,N)
file=fullfile(outdir,sprintf('%s_N%04d_newton',label,N), ...
    'periodic_edge_orbit_newton.mat');
end

function record=summarize_attractor(spec,source_kind,orbit,floquet, ...
        orbit_file,floquet_file,cfg)
nontrivial=true(size(floquet.multipliers)); nontrivial(floquet.phase_index)=false;
largest=max(abs(floquet.multipliers(nontrivial)));
record=empty_orbit_record(); record.label=spec.label; record.N=spec.N;
record.dx=cfg.grid.Lx/(spec.N-1); record.source_kind=source_kind;
record.period=orbit.period; record.newton_residual=orbit.summary.flow_residual;
record.ode45_verification_residual=orbit.summary.ode45_verification_flow_residual;
record.state_correction_from_seed=orbit.summary.state_correction_from_seed;
record.nontrivial_unstable_count=numel(floquet.unstable_indices);
record.largest_nontrivial_multiplier_modulus=largest;
record.phase_multiplier_error=abs(floquet.multipliers(floquet.phase_index)-1);
record.maximum_eigenpair_residual=max(floquet.eigen_residuals);
record.attracting_supported=numel(floquet.unstable_indices)==0 && ...
    largest<1-cfg.BC_resolution.attracting_modulus_margin;
record.orbit_file=orbit_file; record.floquet_file=floquet_file;
end

function record=empty_orbit_record()
record=struct('label','','N',NaN,'dx',NaN,'source_kind','', ...
    'period',NaN,'newton_residual',NaN,'ode45_verification_residual',NaN, ...
    'state_correction_from_seed',NaN,'nontrivial_unstable_count',NaN, ...
    'largest_nontrivial_multiplier_modulus',NaN,'phase_multiplier_error',NaN, ...
    'maximum_eigenpair_residual',NaN,'attracting_supported',false, ...
    'orbit_file','','floquet_file','');
end

function record=empty_branch_record()
record=struct('N',NaN,'sign',NaN,'relative_amplitude',NaN, ...
    'orientation_factor',NaN,'orientation_inner_product',NaN, ...
    'initial_relative_distance',NaN,'outcome','','tested_duration',NaN, ...
    'decision_time',NaN,'late_median_distance_B',NaN, ...
    'late_median_distance_C',NaN,'late_maximum_distance_B',NaN, ...
    'late_maximum_distance_C',NaN,'raw_file','');
end

function reference=edge_direction_reference(edge_resolution,cfg)
index=find([edge_resolution.case_files.N]==cfg.BC_resolution.source_resolution,1);
loaded=load(edge_resolution.case_files(index).floquet_file,'results'); edge=loaded.results;
unstable=edge.unstable_indices;
q=real(edge.weighted_eigenvectors(:,unstable));
w=spatial_weights(edge.state(:,1)); y=q./sqrt([w;w]);
reference=struct('x',edge.state(:,1),'direction_y',y/norm(sqrt([w;w]).*y));
end

function edge=load_edge_case(edge_resolution,N)
index=find([edge_resolution.case_files.N]==N,1);
loaded=load(edge_resolution.case_files(index).floquet_file,'results'); edge=loaded.results;
if numel(edge.unstable_indices)~=1
    error('run_BC_edge_resolution_geometry:EdgeUnstableDimension', ...
        'N=%d edge orbit does not have exactly one unstable mode.',N);
end
end

function attractor=load_attractor(outdir,label,N)
loaded=load(attractor_orbit_file(outdir,label,N),'results'); attractor=loaded.results;
end

function library=exact_BC_library(B,C,count)
template=struct('label','','type','periodic','U',[],'V',[],'period',NaN);
library=repmat(template,2,1); labels={'B','C'}; sources={B,C};
phase=(0:count-1)'/count;
for k=1:2
    source=sources{k}; theta=source.orbit_times(:)/source.period;
    [theta,indices]=unique(theta,'stable');
    library(k).label=labels{k}; library(k).period=source.period;
    library(k).U=interp1(theta,source.orbit_U(indices,:),phase,'pchip');
    library(k).V=interp1(theta,source.orbit_V(indices,:),phase,'pchip');
end
end

function [direction_y,factor,inner_product]=oriented_edge_direction(edge,reference,cfg)
N=size(edge.state,1); w=spatial_weights(edge.state(:,1));
q=real(edge.weighted_eigenvectors(:,edge.unstable_indices)); q=q/norm(q);
direction_y=q./sqrt([w;w]);
x_common=reference.x; x=edge.state(:,1); half=N;
candidate=[interp1(x,direction_y(1:half),x_common,'pchip'); ...
    interp1(x,direction_y(half+1:end),x_common,'pchip')];
w_common=spatial_weights(x_common);
inner_product=sum(w_common.*(candidate(1:numel(x_common)).* ...
    reference.direction_y(1:numel(x_common))+ ...
    candidate(numel(x_common)+1:end).* ...
    reference.direction_y(numel(x_common)+1:end)));
factor=sign(inner_product);
if factor==0
    error('run_BC_edge_resolution_geometry:DirectionOrientation', ...
        'Unable to orient the N=%d unstable direction.',N);
end
direction_y=factor*direction_y;
inner_product=factor*inner_product;
state_norm=norm(sqrt([w;w]).*direction_y);
direction_y=direction_y/state_norm;
if abs(cfg.grid.Lx-(x(end)-x(1)))>1e-10
    error('run_BC_edge_resolution_geometry:DomainMismatch','Edge domain mismatch.');
end
end

function token=sign_token(value)
if value<0,token='minus';else,token='plus';end
end

function [records,order]=sort_orbit_records(records)
labels=string({records.label}); key=1000*(labels=="C")+[records.N];
[~,order]=sort(key); records=records(order);
end

function [records,order]=sort_branch_records(records)
key=10*[records.N]+([records.sign]>0); [~,order]=sort(key); records=records(order);
end

function comparisons=attractor_comparisons(files,cfg)
template=struct('label','','N_coarse',NaN,'N_fine',NaN, ...
    'period_difference',NaN,'full_orbit_distance',NaN, ...
    'temporal_shift_fraction',NaN,'uses_reflection',false);
comparisons=repmat(template,4,1); index=0;
for label=cfg.BC_resolution.outcomes
    selected=files(strcmp({files.label},label{1}));
    [~,order]=sort([selected.N]); selected=selected(order);
    orbits=cell(3,1);
    for k=1:3
        loaded=load(selected(k).orbit_file,'results'); orbits{k}=loaded.results;
    end
    for k=1:2
        index=index+1;
        [distance,shift,reflection]=orbit_distance(orbits{k},orbits{k+1},cfg);
        comparisons(index).label=label{1}; comparisons(index).N_coarse=selected(k).N;
        comparisons(index).N_fine=selected(k+1).N;
        comparisons(index).period_difference=abs(orbits{k}.period-orbits{k+1}.period);
        comparisons(index).full_orbit_distance=distance;
        comparisons(index).temporal_shift_fraction=shift;
        comparisons(index).uses_reflection=reflection;
    end
end
end

function [distance,best_shift,uses_reflection]=orbit_distance(A,B,cfg)
x=linspace(0,cfg.grid.Lx,cfg.BC_resolution.comparison_spatial_points);
M=cfg.BC_resolution.comparison_phase_points; phase=(0:M-1)'/M;
OA=common_orbit(A,x,phase); OB=common_orbit(B,x,phase); w=trap_weights(x);
energy_A=mean(sum((OA.U.^2+OA.V.^2).*w,2));
energy_B=mean(sum((OB.U.^2+OB.V.^2).*w,2)); denominator=0.5*(energy_A+energy_B);
distance=Inf; best_shift=NaN; uses_reflection=false;
for shift=0:M-1
    U=circshift(OB.U,shift,1); V=circshift(OB.V,shift,1);
    direct=sqrt(mean(sum(((OA.U-U).^2+(OA.V-V).^2).*w,2))/denominator);
    reflected=sqrt(mean(sum(((OA.U-fliplr(U)).^2+ ...
        (OA.V-fliplr(V)).^2).*w,2))/denominator);
    if direct<distance,distance=direct;best_shift=shift/M;uses_reflection=false;end
    if reflected<distance,distance=reflected;best_shift=shift/M;uses_reflection=true;end
end
end

function O=common_orbit(source,x,phase)
theta=source.orbit_times(:)/source.period; [theta,indices]=unique(theta,'stable');
U=interp1(theta,source.orbit_U(indices,:),phase,'pchip');
V=interp1(theta,source.orbit_V(indices,:),phase,'pchip');
source_x=source.state(:,1)';
O=struct('U',interp1(source_x,U',x,'pchip')', ...
    'V',interp1(source_x,V',x,'pchip')');
end

function records=make_resolution_records(orbits,comparisons,cfg)
template=struct('label','','period_extrapolate',NaN, ...
    'leading_stable_modulus_extrapolate',NaN, ...
    'coarse_middle_orbit_distance',NaN,'middle_fine_orbit_distance',NaN, ...
    'orbit_distance_ratio',NaN);
records=repmat(template,2,1);
for k=1:2
    label=cfg.BC_resolution.outcomes{k}; selected=orbits(strcmp({orbits.label},label));
    [~,order]=sort([selected.N]); selected=selected(order);
    h=[selected.dx]'; design=[ones(3,1),h.^2];
    period_fit=design\[selected.period]';
    multiplier_fit=design\[selected.largest_nontrivial_multiplier_modulus]';
    pairs=comparisons(strcmp({comparisons.label},label));
    records(k).label=label; records(k).period_extrapolate=period_fit(1);
    records(k).leading_stable_modulus_extrapolate=multiplier_fit(1);
    records(k).coarse_middle_orbit_distance=pairs(1).full_orbit_distance;
    records(k).middle_fine_orbit_distance=pairs(2).full_orbit_distance;
    records(k).orbit_distance_ratio=pairs(1).full_orbit_distance/ ...
        max(pairs(2).full_orbit_distance,eps);
end
end

function records=summarize_branches(branches,cfg)
template=struct('N',NaN,'minus_outcome','','plus_outcome','', ...
    'opposite_BC_supported',false,'maximum_decision_time',NaN);
records=repmat(template,numel(cfg.BC_resolution.resolutions),1);
for k=1:numel(records)
    N=cfg.BC_resolution.resolutions(k); selected=branches([branches.N]==N);
    minus=selected([selected.sign]<0); plus=selected([selected.sign]>0);
    records(k).N=N; records(k).minus_outcome=minus.outcome;
    records(k).plus_outcome=plus.outcome;
    records(k).opposite_BC_supported=isequal(sort({minus.outcome,plus.outcome}),{'B','C'});
    records(k).maximum_decision_time=max([selected.decision_time]);
end
end

function summary=make_summary(orbits,branches,cfg)
pass_orbits=all([orbits.attracting_supported]) && ...
    max([orbits.newton_residual])<cfg.BC_resolution.newton_residual_gate && ...
    max([orbits.ode45_verification_residual])< ...
        cfg.BC_resolution.verification_residual_gate && ...
    max([orbits.phase_multiplier_error])< ...
        cfg.BC_resolution.phase_multiplier_error_gate && ...
    max([orbits.maximum_eigenpair_residual])<cfg.BC_resolution.eigenpair_residual_gate;
pass_branches=all([branches.opposite_BC_supported]);
if pass_orbits&&pass_branches
    status='matched_B_edge_C_geometry_persists_across_tested_resolutions';
else
    status='matched_B_edge_C_resolution_gate_not_passed';
end
summary=struct('status',status,'exact_attractor_case_count',numel(orbits), ...
    'all_BC_orbits_attracting_supported',all([orbits.attracting_supported]), ...
    'all_edge_branch_pairs_reach_opposite_BC',pass_branches, ...
    'maximum_attractor_newton_residual',max([orbits.newton_residual]), ...
    'maximum_attractor_verification_residual', ...
        max([orbits.ode45_verification_residual]), ...
    'maximum_phase_multiplier_error',max([orbits.phase_multiplier_error]), ...
    'maximum_attractor_eigenpair_residual',max([orbits.maximum_eigenpair_residual]), ...
    'claim_scope',cfg.BC_resolution.claim_scope);
end

function make_figures(outdir,orbits,comparisons,branches,cfg)
fig=figure('Color','w','Position',[100 100 1100 800]); layout=tiledlayout(fig,2,2);
colors=lines(2);
for j=1:2
    label=cfg.BC_resolution.outcomes{j}; selected=orbits(strcmp({orbits.label},label));
    [~,order]=sort([selected.N]); selected=selected(order); h2=[selected.dx].^2;
    ax=nexttile(layout); plot(ax,h2,[selected.period],'o-', ...
        'Color',colors(j,:),'MarkerFaceColor',colors(j,:),'LineWidth',1.3);
    xlabel(ax,'dx^2');ylabel(ax,'period');title(ax,[label,' period']);style_axes(ax);
    ax=nexttile(layout); plot(ax,h2,[selected.largest_nontrivial_multiplier_modulus], ...
        'o-','Color',colors(j,:),'MarkerFaceColor',colors(j,:),'LineWidth',1.3);
    yline(ax,1,'k--');xlabel(ax,'dx^2');ylabel(ax,'largest nontrivial |mu|');
    title(ax,[label,' Floquet stability']);style_axes(ax);
end
exportgraphics(fig,fullfile(outdir,'BC_exact_orbit_convergence.png'),'Resolution',300);
close(fig);

fig=figure('Color','w','Position',[100 100 950 550]); ax=axes(fig); hold(ax,'on');
for j=1:numel(comparisons)
    marker='o'; if strcmp(comparisons(j).label,'C'),marker='s';end
    semilogy(ax,comparisons(j).N_fine,comparisons(j).full_orbit_distance,marker, ...
        'MarkerSize',8,'LineWidth',1.5,'DisplayName',sprintf('%s %d--%d', ...
        comparisons(j).label,comparisons(j).N_coarse,comparisons(j).N_fine));
end
xlabel(ax,'fine-grid N');ylabel(ax,'aligned orbit distance');
title(ax,'Exact B/C orbit grid differences');legend(ax);style_axes(ax);
exportgraphics(fig,fullfile(outdir,'BC_orbit_distances.png'),'Resolution',300);close(fig);

fig=figure('Color','w','Position',[100 100 900 520]);ax=axes(fig);hold(ax,'on');
for k=1:numel(branches)
    y=double(strcmp(branches(k).outcome,'C'));
    if branches(k).sign<0
        color=[.75 .45 .25];
    else
        color=[.25 .45 .75];
    end
    scatter(ax,branches(k).N+8*branches(k).sign,y,85, ...
        color,'filled');
end
yticks(ax,[0 1]);yticklabels(ax,{'B','C'});xlabel(ax,'N (sign offset)');
title(ax,'Grid-matched unstable edge branches');style_axes(ax);
exportgraphics(fig,fullfile(outdir,'matched_edge_branch_outcomes.png'),'Resolution',300);
close(fig);
end

function weights=spatial_weights(x)
dx=x(2)-x(1);weights=dx*ones(numel(x),1);weights([1,end])=0.5*dx;
end

function w=trap_weights(x)
dx=x(2)-x(1);w=dx*ones(size(x));w([1,end])=0.5*dx;
end

function distance=state_distance(A,B,Lx)
N=size(A,1);dx=Lx/(N-1);w=dx*ones(N,1);w([1,end])=0.5*dx;
num=sum(((A(:,2)-B(:,2)).^2+(A(:,3)-B(:,3)).^2).*w);
den=0.5*(sum((A(:,2).^2+A(:,3).^2).*w)+sum((B(:,2).^2+B(:,3).^2).*w));
distance=sqrt(num/max(den,eps));
end

function style_axes(ax)
set(ax,'Color','w','XColor','k','YColor','k','GridColor',[.75 .75 .75]);
grid(ax,'on');box(ax,'on');ax.Title.Color='k';ax.XLabel.Color='k';ax.YLabel.Color='k';
end
