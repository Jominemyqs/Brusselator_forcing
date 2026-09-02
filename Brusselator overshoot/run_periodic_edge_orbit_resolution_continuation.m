function results = run_periodic_edge_orbit_resolution_continuation()
%RUN_PERIODIC_EDGE_ORBIT_RESOLUTION_CONTINUATION Continue exact edge orbit in N.
%   Regrids the Newton-converged N=400 orbit to N=200 and N=600, reconverges
%   the shooting equation, computes the leading Floquet spectrum, and
%   compares phase-fixed states and complete periodic orbits in a common
%   physical norm. The N=400 exact results are reused without modification.

cfg = brusselator_periodic_edge_resolution_config();
outdir = fullfile(cfg.output.root,cfg.experiment_name);
checkpoint_file = fullfile(outdir,'resolution_continuation_checkpoint.mat');
final_file = fullfile(outdir,'periodic_edge_orbit_resolution_continuation.mat');
if isfile(final_file)
    error('run_periodic_edge_orbit_resolution_continuation:OutputExists', ...
        'Refusing to overwrite completed result: %s',final_file);
end
required = {cfg.resolution_continuation.source_orbit_file, ...
    cfg.resolution_continuation.source_floquet_file};
if any(~cellfun(@isfile,required))
    error('run_periodic_edge_orbit_resolution_continuation:MissingSource', ...
        'The N=400 Newton or Floquet source is missing.');
end
source_orbit_loaded = load(required{1},'results');
source_orbit = source_orbit_loaded.results;
source_floquet_loaded = load(required{2},'results');
source_floquet = source_floquet_loaded.results;
if ~strcmp(source_orbit.status,'newton_converged') || ...
        numel(source_floquet.unstable_indices)~=1
    error('run_periodic_edge_orbit_resolution_continuation:InvalidSource', ...
        'Source orbit must be Newton converged with one unstable multiplier.');
end

if isfile(checkpoint_file)
    saved = load(checkpoint_file,'case_records','case_files','next_case');
    case_records = saved.case_records;
    case_files = saved.case_files;
    next_case = saved.next_case;
    fprintf('Resuming exact-orbit resolution continuation at case %d.\n',next_case);
else
    if isfolder(outdir)
        error('run_periodic_edge_orbit_resolution_continuation:IncompleteOutput', ...
            'Existing output lacks a checkpoint: %s',outdir);
    end
    mkdir(outdir); mkdir(fullfile(outdir,'seeds'));
    record_template = empty_case_record();
    case_records = record_template([]);
    case_files = struct('N',{},'orbit_file',{},'floquet_file',{});
    next_case = 1;
    brusselator_write_metadata(outdir,cfg, ...
        brusselator_run_metadata(cfg,mfilename,NaN));
    save_checkpoint();
end
timer = tic;
resolutions = cfg.resolution_continuation.resolutions;

for case_index = next_case:numel(resolutions)
    target_N = resolutions(case_index);
    fprintf('\nExact edge-orbit continuation: N=%d.\n',target_N);
    if target_N == cfg.resolution_continuation.source_resolution
        orbit = source_orbit;
        floquet = source_floquet;
        orbit_file = required{1};
        floquet_file = required{2};
        source_kind = 'reused_exact_source';
    else
        seed_file = make_target_seed( ...
            source_orbit,target_N,cfg,outdir);
        newton_cfg = brusselator_periodic_edge_newton_config();
        newton_cfg.grid.N = target_N;
        newton_cfg.output.root = outdir;
        newton_cfg.experiment_name = sprintf('N%04d_newton',target_N);
        newton_cfg.newton.seed_file = seed_file;
        orbit_file = fullfile(outdir,newton_cfg.experiment_name, ...
            'periodic_edge_orbit_newton.mat');
        if isfile(orbit_file)
            loaded = load(orbit_file,'results'); orbit = loaded.results;
        else
            orbit = run_periodic_edge_orbit_newton(newton_cfg);
        end
        if ~strcmp(orbit.status,'newton_converged')
            error('run_periodic_edge_orbit_resolution_continuation:NewtonFailure', ...
                'The N=%d shooting solve did not converge.',target_N);
        end
        floquet_cfg = brusselator_periodic_edge_floquet_config();
        floquet_cfg.grid.N = target_N;
        floquet_cfg.output.root = outdir;
        floquet_cfg.experiment_name = sprintf('N%04d_floquet',target_N);
        floquet_cfg.floquet.orbit_file = orbit_file;
        floquet_file = fullfile(outdir,floquet_cfg.experiment_name, ...
            'periodic_edge_orbit_floquet.mat');
        if isfile(floquet_file)
            loaded = load(floquet_file,'results'); floquet = loaded.results;
        else
            floquet = run_periodic_edge_orbit_floquet(floquet_cfg);
        end
        source_kind = 'regridded_and_newton_converged';
    end

    record = summarize_case(target_N,source_kind,orbit,floquet, ...
        orbit_file,floquet_file,cfg);
    case_records(end+1,1) = record; %#ok<AGROW>
    case_files(end+1,1) = struct('N',target_N,'orbit_file',orbit_file, ...
        'floquet_file',floquet_file); %#ok<AGROW>
    next_case = case_index+1;
    writetable(struct2table(case_records,'AsArray',true), ...
        fullfile(outdir,'resolution_progress.csv'));
    save_checkpoint();
end

[case_records,order] = sort_records(case_records);
case_files = case_files(order);
orbit_results = cell(numel(case_files),1);
floquet_results = cell(numel(case_files),1);
for k = 1:numel(case_files)
    loaded = load(case_files(k).orbit_file,'results'); orbit_results{k}=loaded.results;
    loaded = load(case_files(k).floquet_file,'results'); floquet_results{k}=loaded.results;
end
comparison_records = compare_resolutions(orbit_results,case_records,cfg);
convergence = convergence_summary(case_records,comparison_records,cfg);
status = determine_status(case_records,cfg);
summary = struct('status',status,'resolution_count',numel(case_records), ...
    'minimum_resolution',min([case_records.N]), ...
    'maximum_resolution',max([case_records.N]), ...
    'period_range',max([case_records.period])-min([case_records.period]), ...
    'unstable_multiplier_range',max([case_records.unstable_multiplier])- ...
        min([case_records.unstable_multiplier]), ...
    'maximum_newton_residual',max([case_records.newton_residual]), ...
    'maximum_phase_multiplier_error',max([case_records.phase_multiplier_error]), ...
    'maximum_eigenpair_residual',max([case_records.maximum_eigenpair_residual]), ...
    'expected_order_period_extrapolate',convergence.period_extrapolate, ...
    'expected_order_multiplier_extrapolate',convergence.multiplier_extrapolate, ...
    'interpretation',cfg.resolution_continuation.interpretation);
writetable(struct2table(case_records,'AsArray',true), ...
    fullfile(outdir,'resolution_continuation_summary.csv'));
writetable(struct2table(comparison_records,'AsArray',true), ...
    fullfile(outdir,'pairwise_orbit_comparisons.csv'));
writetable(struct2table(convergence,'AsArray',true), ...
    fullfile(outdir,'convergence_diagnostics.csv'));
writetable(struct2table(summary,'AsArray',true), ...
    fullfile(outdir,'study_summary.csv'));
make_figures(outdir,case_records,orbit_results,comparison_records,convergence,cfg);
results = struct('configuration',cfg,'case_records',case_records, ...
    'case_files',case_files,'comparison_records',comparison_records, ...
    'convergence',convergence,'summary',summary);
save(final_file,'results','-v7.3');
brusselator_write_metadata(outdir,cfg, ...
    brusselator_run_metadata(cfg,mfilename,toc(timer)));
disp(struct2table(case_records,'AsArray',true));
disp(struct2table(comparison_records,'AsArray',true));
disp(struct2table(convergence,'AsArray',true));
disp(struct2table(summary,'AsArray',true));
fprintf('Exact edge-orbit resolution continuation saved in: %s\n',outdir);

    function save_checkpoint()
        save(checkpoint_file,'cfg','case_records','case_files','next_case','-v7.3');
    end
end

function seed_file = make_target_seed(source,target_N,cfg,outdir)
target_cfg = brusselator_periodic_edge_newton_config();
target_cfg.grid.N = target_N;
target_state = regrid_state(source.state,target_N,cfg.grid.Lx);
target_tangent = brusselator_frozen_rhs(target_state,target_cfg, ...
    target_cfg.edge_tracking.frozen_b);
seed = struct('state',target_state,'period',source.period, ...
    'event_index',NaN,'reference_state',target_state, ...
    'reference_tangent',target_tangent, ...
    'phase_condition','<X-X_ref,F_target_grid(X_ref)>_w=0', ...
    'flow_residual',NaN,'phase_condition_residual',0, ...
    'status','regridded_exact_orbit_continuation_seed', ...
    'source_resolution',size(source.state,1),'target_resolution',target_N, ...
    'regridding',cfg.resolution_continuation.regridding);
seed_file = fullfile(outdir,'seeds',sprintf('N%04d_seed.mat',target_N));
if ~isfile(seed_file)
    save(seed_file,'seed','-v7.3');
else
    loaded = load(seed_file,'seed');
    if loaded.seed.target_resolution~=target_N
        error('run_periodic_edge_orbit_resolution_continuation:SeedMismatch', ...
            'Existing seed file has the wrong target resolution.');
    end
end
end

function state = regrid_state(source,target_N,Lx)
x = linspace(0,Lx,target_N)';
state = [x,interp1(source(:,1),source(:,2),x,'pchip'), ...
    interp1(source(:,1),source(:,3),x,'pchip')];
end

function record = summarize_case(N,source_kind,orbit,floquet,orbit_file, ...
        floquet_file,cfg)
unstable = floquet.unstable_indices;
if isscalar(unstable)
    multiplier = real(floquet.multipliers(unstable));
    exponent = log(abs(multiplier))/orbit.period;
else
    multiplier = NaN; exponent = NaN;
end
record = empty_case_record();
record.N = N; record.dx = cfg.grid.Lx/(N-1); record.source_kind = source_kind;
record.period = orbit.period;
record.newton_residual = orbit.summary.flow_residual;
record.ode45_verification_residual = orbit.summary.ode45_verification_flow_residual;
record.state_correction_from_seed = orbit.summary.state_correction_from_seed;
record.unstable_count = numel(unstable);
record.unstable_multiplier = multiplier;
record.unstable_exponent = exponent;
record.phase_multiplier = real(floquet.multipliers(floquet.phase_index));
record.phase_multiplier_error = abs(record.phase_multiplier-1);
record.maximum_eigenpair_residual = max(floquet.eigen_residuals);
record.tangent_finite_difference_error = ...
    floquet.finite_difference_check.relative_error;
record.orbit_file = orbit_file; record.floquet_file = floquet_file;
end

function record = empty_case_record()
record = struct('N',NaN,'dx',NaN,'source_kind','', ...
    'period',NaN,'newton_residual',NaN,'ode45_verification_residual',NaN, ...
    'state_correction_from_seed',NaN,'unstable_count',NaN, ...
    'unstable_multiplier',NaN,'unstable_exponent',NaN, ...
    'phase_multiplier',NaN,'phase_multiplier_error',NaN, ...
    'maximum_eigenpair_residual',NaN,'tangent_finite_difference_error',NaN, ...
    'orbit_file','','floquet_file','');
end

function [records,order] = sort_records(records)
[~,order] = sort([records.N]); records = records(order);
end

function comparisons = compare_resolutions(orbits,records,cfg)
template = struct('N_coarse',NaN,'N_fine',NaN,'dx_coarse',NaN,'dx_fine',NaN, ...
    'period_difference',NaN,'phase_fixed_state_distance',NaN, ...
    'phase_fixed_uses_reflection',false,'full_orbit_distance',NaN, ...
    'full_orbit_temporal_shift_fraction',NaN,'full_orbit_uses_reflection',false);
comparisons = repmat(template,numel(orbits)-1,1);
x_common = linspace(0,cfg.grid.Lx, ...
    cfg.resolution_continuation.comparison_spatial_points);
phase = (0:cfg.resolution_continuation.comparison_phase_points-1)' / ...
    cfg.resolution_continuation.comparison_phase_points;
common = cell(numel(orbits),1);
for k = 1:numel(orbits)
    common{k} = orbit_on_common_grid(orbits{k},x_common,phase);
end
for k = 1:numel(orbits)-1
    [state_distance,uses_reflection] = phase_fixed_distance( ...
        common{k}.U(1,:),common{k}.V(1,:),common{k+1}.U(1,:), ...
        common{k+1}.V(1,:),x_common);
    [orbit_distance,shift,orbit_reflection] = aligned_orbit_distance( ...
        common{k},common{k+1},x_common);
    comparisons(k).N_coarse=records(k).N;
    comparisons(k).N_fine=records(k+1).N;
    comparisons(k).dx_coarse=records(k).dx;
    comparisons(k).dx_fine=records(k+1).dx;
    comparisons(k).period_difference=abs(records(k).period-records(k+1).period);
    comparisons(k).phase_fixed_state_distance=state_distance;
    comparisons(k).phase_fixed_uses_reflection=uses_reflection;
    comparisons(k).full_orbit_distance=orbit_distance;
    comparisons(k).full_orbit_temporal_shift_fraction=shift;
    comparisons(k).full_orbit_uses_reflection=orbit_reflection;
end
end

function common = orbit_on_common_grid(orbit,x_common,phase)
theta = orbit.orbit_times(:)/orbit.period;
[theta,unique_indices] = unique(theta,'stable');
U = interp1(theta,orbit.orbit_U(unique_indices,:),phase,'pchip');
V = interp1(theta,orbit.orbit_V(unique_indices,:),phase,'pchip');
x = orbit.state(:,1)';
common = struct('U',interp1(x,U',x_common,'pchip')', ...
    'V',interp1(x,V',x_common,'pchip')');
end

function [distance,uses_reflection] = phase_fixed_distance(u1,v1,u2,v2,x)
w = trap_weights(x);
denominator = max(0.5*(sum((u1.^2+v1.^2).*w)+ ...
    sum((u2.^2+v2.^2).*w)),eps);
direct = sqrt(sum(((u1-u2).^2+(v1-v2).^2).*w)/denominator);
reflected = sqrt(sum(((u1-fliplr(u2)).^2+(v1-fliplr(v2)).^2).*w)/denominator);
uses_reflection = reflected<direct; distance=min(direct,reflected);
end

function [distance,best_shift,uses_reflection] = aligned_orbit_distance(A,B,x)
w = trap_weights(x); M = size(A.U,1);
energy_A = mean(sum((A.U.^2+A.V.^2).*w,2));
energy_B = mean(sum((B.U.^2+B.V.^2).*w,2));
denominator = max(0.5*(energy_A+energy_B),eps);
best = Inf; best_shift=NaN; uses_reflection=false;
for shift = 0:M-1
    U = circshift(B.U,shift,1); V = circshift(B.V,shift,1);
    direct = sqrt(mean(sum(((A.U-U).^2+(A.V-V).^2).*w,2))/denominator);
    reflected = sqrt(mean(sum(((A.U-fliplr(U)).^2+ ...
        (A.V-fliplr(V)).^2).*w,2))/denominator);
    if direct<best
        best=direct; best_shift=shift/M; uses_reflection=false;
    end
    if reflected<best
        best=reflected; best_shift=shift/M; uses_reflection=true;
    end
end
distance=best;
end

function w = trap_weights(x)
dx=x(2)-x(1); w=dx*ones(size(x)); w([1,end])=0.5*dx;
end

function convergence = convergence_summary(records,comparisons,cfg)
h = [records.dx]'; periods=[records.period]'; multipliers=[records.unstable_multiplier]';
order = cfg.resolution_continuation.expected_spatial_order;
design = [ones(numel(h),1),h.^order];
period_fit = design\periods; multiplier_fit = design\multipliers;
period_prediction=design*period_fit; multiplier_prediction=design*multiplier_fit;
convergence = struct('expected_spatial_order',order, ...
    'period_extrapolate',period_fit(1),'period_h2_coefficient',period_fit(2), ...
    'period_fit_maximum_residual',max(abs(periods-period_prediction)), ...
    'multiplier_extrapolate',multiplier_fit(1), ...
    'multiplier_h2_coefficient',multiplier_fit(2), ...
    'multiplier_fit_maximum_residual',max(abs(multipliers-multiplier_prediction)), ...
    'coarse_to_middle_orbit_distance',comparisons(1).full_orbit_distance, ...
    'middle_to_fine_orbit_distance',comparisons(2).full_orbit_distance, ...
    'orbit_distance_ratio',comparisons(1).full_orbit_distance / ...
        max(comparisons(2).full_orbit_distance,eps), ...
    'interpretation',['h^2 extrapolates are diagnostics consistent with the ', ...
        'second-order spatial stencil, not proof of asymptotic convergence']);
end

function status = determine_status(records,cfg)
pass = all(strcmp({records.source_kind},'reused_exact_source') | ...
    strcmp({records.source_kind},'regridded_and_newton_converged')) && ...
    all([records.newton_residual]<cfg.resolution_continuation.required_newton_residual) && ...
    all([records.phase_multiplier_error]< ...
        cfg.resolution_continuation.required_phase_multiplier_error) && ...
    all([records.maximum_eigenpair_residual]< ...
        cfg.resolution_continuation.required_eigenpair_residual) && ...
    all([records.unstable_count]==1);
if pass
    status='exact_periodic_edge_orbit_persists_across_tested_resolutions';
else
    status='resolution_continuation_gate_not_passed';
end
end

function make_figures(outdir,records,orbits,comparisons,convergence,cfg)
h2=[records.dx].^2; periods=[records.period]; multipliers=[records.unstable_multiplier];
h2line=linspace(0,max(h2)*1.05,100);
fig=figure('Color','w','Position',[100 100 1100 850]); layout=tiledlayout(fig,2,2);
ax=nexttile(layout); hold(ax,'on'); plot(ax,h2,periods,'o','MarkerFaceColor',[.2 .5 .8]);
plot(ax,h2line,convergence.period_extrapolate+ ...
    convergence.period_h2_coefficient*h2line,'-k');
xlabel(ax,'dx^2'); ylabel(ax,'period'); title(ax,'Period convergence'); style_axes(ax);
ax=nexttile(layout); hold(ax,'on'); plot(ax,h2,multipliers,'o','MarkerFaceColor',[.85 .4 .15]);
plot(ax,h2line,convergence.multiplier_extrapolate+ ...
    convergence.multiplier_h2_coefficient*h2line,'-k');
xlabel(ax,'dx^2'); ylabel(ax,'unstable multiplier'); title(ax,'Floquet convergence'); style_axes(ax);
ax=nexttile(layout); hold(ax,'on'); colors=lines(numel(orbits));
for k=1:numel(orbits)
    plot(ax,orbits{k}.state(:,1),orbits{k}.state(:,2),'Color',colors(k,:), ...
        'LineWidth',1.2,'DisplayName',sprintf('N=%d',records(k).N));
end
xlabel(ax,'x'); ylabel(ax,'u'); title(ax,'Phase-fixed orbit state'); legend(ax); style_axes(ax);
ax=nexttile(layout); semilogy(ax,[comparisons.N_fine], ...
    [comparisons.full_orbit_distance],'o-','LineWidth',1.3, ...
    'MarkerFaceColor',[.3 .65 .4]);
xlabel(ax,'fine-grid N'); ylabel(ax,'aligned orbit distance');
title(ax,'Pairwise full-orbit differences'); style_axes(ax);
exportgraphics(fig,fullfile(outdir,'resolution_continuation.png'),'Resolution',300);
close(fig);

fig=figure('Color','w','Position',[100 100 950 650]); layout=tiledlayout(fig,2,1);
component_labels={'u','v'};
for component=1:2
    ax=nexttile(layout); hold(ax,'on');
    for k=1:numel(orbits)
        values=orbits{k}.state(:,component+1);
        plot(ax,orbits{k}.state(:,1),values,'Color',colors(k,:), ...
            'LineWidth',1.1,'DisplayName',sprintf('N=%d',records(k).N));
    end
    ylabel(ax,component_labels{component}); style_axes(ax);
    if component==1, title(ax,sprintf('Exact edge orbit, L=%g',cfg.grid.Lx)); legend(ax); end
end
xlabel(nexttile(layout,2),'x');
exportgraphics(fig,fullfile(outdir,'phase_fixed_profiles.png'),'Resolution',300);
close(fig);
end

function style_axes(ax)
set(ax,'Color','w','XColor','k','YColor','k','GridColor',[.75 .75 .75]);
grid(ax,'on'); box(ax,'on'); ax.Title.Color='k';
ax.XLabel.Color='k'; ax.YLabel.Color='k';
end
