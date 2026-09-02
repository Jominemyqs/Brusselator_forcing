function results = analyze_dynamic_edge_BC_fixed_phase()
%ANALYZE_DYNAMIC_EDGE_BC_FIXED_PHASE Audit edge templates on one section.
%   Rephases every completed B--C edge segment on the single hyperplane
%
%       <X-X_ref,F(X_ref)>_w = 0,
%
%   using upward crossings. B- and C-side exact trajectories are phase fixed
%   independently; their matched section states are averaged only after the
%   crossing calculation. This removes the per-segment centering drift in
%   the exploratory mean-v Poincare sections and produces an auditable seed
%   for the periodic shooting equation. It does not solve that equation.

cfg = brusselator_dynamic_edge_multicycle_config();
cfg.experiment_name = 'dynamic_edge_BC_fixed_phase_convergence_v1';
cfg.fixed_phase = struct( ...
    'source_file', fullfile('experiment_outputs', ...
        'dynamic_edge_tracking_BC_multicycle_v1', ...
        'dynamic_edge_tracking_BC_multicycle.mat'), ...
    'parent_segment_file', cfg.multicycle.parent_post_segment_file, ...
    'reference_event', 4, ...
    'minimum_time', cfg.multicycle.phase_sample_delay, ...
    'minimum_relative_transversality', 1e-3, ...
    'maximum_period_relative_deviation', 0.2, ...
    'phase_condition', ...
        '<X-X_ref,F(X_ref)>_w=0 with positive crossing derivative', ...
    'template_rule', ['average the separately phase-fixed exact B- and C-side ', ...
        'states at their first admissible fixed-section crossings'], ...
    'seed_rule', 'smallest direct one-period residual among audited events');

outdir = fullfile(cfg.output.root,cfg.experiment_name);
final_file = fullfile(outdir,'fixed_phase_convergence.mat');
if isfile(final_file)
    error('analyze_dynamic_edge_BC_fixed_phase:OutputExists', ...
        'Refusing to overwrite completed result: %s',final_file);
end
if ~isfile(cfg.fixed_phase.source_file) || ...
        ~isfile(cfg.fixed_phase.parent_segment_file)
    error('analyze_dynamic_edge_BC_fixed_phase:MissingSource', ...
        'The completed multicycle result or its parent segment is missing.');
end
mkdir(outdir);
timer = tic;

loaded = load(cfg.fixed_phase.source_file,'results');
source = loaded.results;
segments = load_segments(source,cfg);
reference_state = source.template_states{cfg.fixed_phase.reference_event};
F_reference = brusselator_frozen_rhs(reference_state,cfg, ...
    cfg.edge_tracking.frozen_b);
weights = spatial_weights(reference_state(:,1));

n_events = numel(segments);
record_template = empty_record();
records = repmat(record_template,n_events,1);
templates = cell(n_events,1);
side_diagnostics = cell(n_events,2);

for k = 1:n_events
    segment = segments{k};
    expected_period = source.template_records(k).period_mean;
    B = fixed_phase_side(segment.S,segment.U_B,segment.V_B, ...
        reference_state,F_reference,weights,expected_period,cfg);
    C = fixed_phase_side(segment.S,segment.U_C,segment.V_C, ...
        reference_state,F_reference,weights,expected_period,cfg);
    side_diagnostics{k,1} = B;
    side_diagnostics{k,2} = C;

    template = B.phase_state;
    template(:,2:3) = 0.5*(B.phase_state(:,2:3)+C.phase_state(:,2:3));
    period = 0.5*(B.period+C.period);
    [flow_residual,reflected_flow_residual] = compute_flow_residual( ...
        template,period,cfg);
    phase_residual = phase_value(template,reference_state,F_reference,weights);
    phase_pair_distance = state_distance(B.phase_state,C.phase_state,cfg.grid.Lx);
    reference_distance = state_distance(template,reference_state,cfg.grid.Lx);

    record = record_template;
    record.event_index = k;
    record.period_B = B.period;
    record.period_C = C.period;
    record.period_mean = period;
    record.period_difference = abs(B.period-C.period);
    record.period_cv_B = B.period_cv;
    record.period_cv_C = C.period_cv;
    record.crossing_count_B = B.crossing_count;
    record.crossing_count_C = C.crossing_count;
    record.first_crossing_time_B = B.phase_time;
    record.first_crossing_time_C = C.phase_time;
    record.phase_time_difference = abs(B.phase_time-C.phase_time);
    record.minimum_relative_transversality_B = min(B.relative_transversality);
    record.minimum_relative_transversality_C = min(C.relative_transversality);
    record.maximum_period_relative_deviation_B = B.maximum_period_relative_deviation;
    record.maximum_period_relative_deviation_C = C.maximum_period_relative_deviation;
    record.unique_branch_B = B.unique_branch;
    record.unique_branch_C = C.unique_branch;
    record.phase_pair_distance = phase_pair_distance;
    record.phase_condition_residual = phase_residual;
    record.reference_distance = reference_distance;
    record.flow_residual = flow_residual;
    record.reflected_flow_residual = reflected_flow_residual;
    if k > 1
        [best,direct,reflected,uses_reflection] = template_distance( ...
            templates{k-1},template,cfg.grid.Lx);
        record.previous_template_distance = best;
        record.previous_template_direct_distance = direct;
        record.previous_template_reflected_distance = reflected;
        record.previous_template_uses_reflection = uses_reflection;
        record.period_change = abs(period-records(k-1).period_mean);
    end
    record.audit_pass = B.audit_pass && C.audit_pass;
    records(k) = record;
    templates{k} = template;
end

eligible = find([records.audit_pass]);
if isempty(eligible)
    error('analyze_dynamic_edge_BC_fixed_phase:NoAuditedSeed', ...
        'No event passed the fixed-section crossing audit.');
end
[~,local_best] = min([records(eligible).flow_residual]);
best_event = eligible(local_best);
seed = struct('state',templates{best_event}, ...
    'period',records(best_event).period_mean, ...
    'event_index',best_event,'reference_state',reference_state, ...
    'reference_tangent',F_reference, ...
    'phase_condition',cfg.fixed_phase.phase_condition, ...
    'flow_residual',records(best_event).flow_residual, ...
    'phase_condition_residual',records(best_event).phase_condition_residual, ...
    'status','fixed_phase_shooting_seed_not_newton_converged');

summary = make_summary(records,best_event,source);
writetable(struct2table(records,'AsArray',true), ...
    fullfile(outdir,'fixed_phase_template_convergence.csv'));
writetable(struct2table(summary,'AsArray',true), ...
    fullfile(outdir,'fixed_phase_summary.csv'));
save(fullfile(outdir,'fixed_phase_newton_seed.mat'),'seed','-v7.3');
make_figure(outdir,records);
results = struct('configuration',cfg,'source_file',cfg.fixed_phase.source_file, ...
    'reference_state',reference_state,'reference_tangent',F_reference, ...
    'records',records,'templates',{templates}, ...
    'side_diagnostics',{side_diagnostics},'seed',seed,'summary',summary);
save(final_file,'results','-v7.3');
brusselator_write_metadata(outdir,cfg, ...
    brusselator_run_metadata(cfg,mfilename,toc(timer)));
disp(struct2table(records,'AsArray',true));
disp(struct2table(summary,'AsArray',true));
fprintf('Fixed-phase edge analysis saved in: %s\n',outdir);
end

function segments = load_segments(source,cfg)
segments = cell(numel(source.template_records),1);
loaded = load(cfg.fixed_phase.parent_segment_file,'segment_data');
segments{1} = loaded.segment_data;
for k = 2:numel(segments)
    file = source.event_records(k-1).raw_file;
    if ~isfile(file)
        error('analyze_dynamic_edge_BC_fixed_phase:MissingSegment', ...
            'Missing event segment: %s',file);
    end
    loaded = load(file,'segment_data');
    segments{k} = loaded.segment_data;
end
end

function diagnostics = fixed_phase_side(S,U,V,reference,Fref,weights, ...
        expected_period,cfg)
S = S(:);
valid = S >= cfg.fixed_phase.minimum_time;
S = S(valid); U = U(valid,:); V = V(valid,:);
g = zeros(numel(S),1);
for j = 1:numel(S)
    state = [reference(:,1),U(j,:)',V(j,:)'];
    g(j) = phase_value(state,reference,Fref,weights);
end
indices = find(g(1:end-1) <= 0 & g(2:end) > 0);
if numel(indices) < 4
    error('analyze_dynamic_edge_BC_fixed_phase:FewCrossings', ...
        'A side trajectory has fewer than four upward fixed-section crossings.');
end

n = numel(indices);
times = zeros(n,1); transversality = zeros(n,1);
states = cell(n,1); reference_distances = zeros(n,1);
for j = 1:n
    bracket = S(indices(j):indices(j)+1);
    phase_function = @(t) interpolated_phase_value(t,S,U,V, ...
        reference,Fref,weights);
    times(j) = fzero(phase_function,bracket);
    u = interp1(S,U,times(j),'pchip');
    v = interp1(S,V,times(j),'pchip');
    state = [reference(:,1),u',v'];
    states{j} = state;
    F = brusselator_frozen_rhs(state,cfg,cfg.edge_tracking.frozen_b);
    transversality(j) = physical_inner(F,Fref,weights) / ...
        max(physical_inner(Fref,Fref,weights),eps);
    reference_distances(j) = state_distance(state,reference,cfg.grid.Lx);
end
periods = diff(times);
period = median(periods);
period_cv = std(periods)/max(mean(periods),eps);
maximum_period_relative_deviation = max(abs(periods-expected_period)) / ...
    expected_period;
positive_transverse = transversality >= ...
    cfg.fixed_phase.minimum_relative_transversality;
regular_spacing = maximum_period_relative_deviation <= ...
    cfg.fixed_phase.maximum_period_relative_deviation;
unique_branch = regular_spacing && all(positive_transverse);
diagnostics = struct('phase_time',times(1),'phase_state',states{1}, ...
    'times',times,'states',{states},'periods',periods,'period',period, ...
    'period_cv',period_cv,'crossing_count',n, ...
    'relative_transversality',transversality, ...
    'reference_distances',reference_distances, ...
    'maximum_period_relative_deviation',maximum_period_relative_deviation, ...
    'unique_branch',unique_branch, ...
    'audit_pass',unique_branch && period_cv < 0.05);
end

function value = interpolated_phase_value(t,S,U,V,reference,Fref,weights)
u = interp1(S,U,t,'pchip');
v = interp1(S,V,t,'pchip');
state = [reference(:,1),u',v'];
value = phase_value(state,reference,Fref,weights);
end

function value = phase_value(state,reference,Fref,weights)
value = physical_inner(state(:,2:3)-reference(:,2:3),Fref,weights);
end

function value = physical_inner(A,B,weights)
value = sum(sum(A.*B,2).*weights);
end

function weights = spatial_weights(x)
dx = x(2)-x(1);
weights = dx*ones(numel(x),1);
weights([1,end]) = 0.5*dx;
end

function [direct,reflected] = compute_flow_residual(state,period,cfg)
par = brusselator_make_parameters(cfg,@(t)cfg.edge_tracking.frozen_b);
final_state = solve_brusselator_1d_forced(state,par,period,0);
direct = state_distance(state,final_state,cfg.grid.Lx);
reflected = reflected_state_distance(state,final_state,cfg.grid.Lx);
end

function [best,direct,reflected,uses_reflection] = template_distance(A,B,Lx)
direct = state_distance(A,B,Lx);
reflected = reflected_state_distance(A,B,Lx);
best = min(direct,reflected);
uses_reflection = reflected < direct;
end

function distance = state_distance(A,B,Lx)
[distance,~] = distance_components(A,B,Lx);
end

function distance = reflected_state_distance(A,B,Lx)
[~,distance] = distance_components(A,B,Lx);
end

function [direct,reflected] = distance_components(A,B,Lx)
N = size(A,1); dx = Lx/(N-1); w = dx*ones(N,1); w([1,end]) = 0.5*dx;
energy_A = sum((A(:,2).^2+A(:,3).^2).*w);
energy_B = sum((B(:,2).^2+B(:,3).^2).*w);
denominator = max(0.5*(energy_A+energy_B),eps);
direct = sqrt(sum(((A(:,2)-B(:,2)).^2+(A(:,3)-B(:,3)).^2).*w)/denominator);
reflected = sqrt(sum(((A(:,2)-flipud(B(:,2))).^2+ ...
    (A(:,3)-flipud(B(:,3))).^2).*w)/denominator);
end

function record = empty_record()
record = struct('event_index',NaN,'period_B',NaN,'period_C',NaN, ...
    'period_mean',NaN,'period_difference',NaN,'period_cv_B',NaN, ...
    'period_cv_C',NaN,'crossing_count_B',NaN,'crossing_count_C',NaN, ...
    'first_crossing_time_B',NaN,'first_crossing_time_C',NaN, ...
    'phase_time_difference',NaN,'minimum_relative_transversality_B',NaN, ...
    'minimum_relative_transversality_C',NaN, ...
    'maximum_period_relative_deviation_B',NaN, ...
    'maximum_period_relative_deviation_C',NaN,'unique_branch_B',false, ...
    'unique_branch_C',false,'phase_pair_distance',NaN, ...
    'phase_condition_residual',NaN,'reference_distance',NaN, ...
    'flow_residual',NaN,'reflected_flow_residual',NaN, ...
    'previous_template_distance',NaN,'previous_template_direct_distance',NaN, ...
    'previous_template_reflected_distance',NaN, ...
    'previous_template_uses_reflection',false,'period_change',NaN, ...
    'audit_pass',false);
end

function summary = make_summary(records,best_event,source)
last = records(end);
later = records(max(2,numel(records)-2):end);
summary = struct('status','fixed_phase_audit_complete', ...
    'event_count',numel(records),'total_shadow_time',source.summary.total_shadow_time, ...
    'all_events_pass_crossing_audit',all([records.audit_pass]), ...
    'final_period',last.period_mean,'last_period_change',last.period_change, ...
    'last_template_distance',last.previous_template_distance, ...
    'maximum_last_three_template_distance', ...
        max([later.previous_template_distance]), ...
    'best_seed_event',best_event, ...
    'best_seed_period',records(best_event).period_mean, ...
    'best_seed_flow_residual',records(best_event).flow_residual, ...
    'best_seed_phase_residual',records(best_event).phase_condition_residual, ...
    'interpretation',['fixed-section shooting seed only; no Newton convergence ', ...
        'or Floquet multiplier has yet been established']);
end

function make_figure(outdir,records)
events = [records.event_index];
fig = figure('Color','w','Position',[100 100 1000 820]);
layout = tiledlayout(fig,4,1);
ax = nexttile(layout); plot(ax,events,[records.period_mean],'o-', ...
    'LineWidth',1.4,'MarkerFaceColor',[0.2 0.5 0.8]);
ylabel(ax,'period'); title(ax,'Fixed-hyperplane edge-template audit'); style_axes(ax);
ax = nexttile(layout); semilogy(ax,events(2:end), ...
    [records(2:end).previous_template_distance],'o-','LineWidth',1.4, ...
    'MarkerFaceColor',[0.85 0.45 0.15]); ylabel(ax,'template distance'); style_axes(ax);
ax = nexttile(layout); semilogy(ax,events,[records.flow_residual],'o-', ...
    'LineWidth',1.4,'MarkerFaceColor',[0.25 0.65 0.35]);
ylabel(ax,'flow residual'); style_axes(ax);
ax = nexttile(layout); plot(ax,events, ...
    min([[records.minimum_relative_transversality_B]; ...
         [records.minimum_relative_transversality_C]],[],1),'o-', ...
    'LineWidth',1.4,'MarkerFaceColor',[0.55 0.35 0.75]);
xlabel(ax,'post-rebracketing cycle'); ylabel(ax,'min. relative transversality');
style_axes(ax);
exportgraphics(fig,fullfile(outdir,'fixed_phase_template_convergence.png'), ...
    'Resolution',300); close(fig);
end

function style_axes(ax)
set(ax,'Color','w','XColor','k','YColor','k','GridColor',[0.75 0.75 0.75]);
grid(ax,'on'); box(ax,'on'); ax.Title.Color='k';
ax.XLabel.Color='k'; ax.YLabel.Color='k';
end
