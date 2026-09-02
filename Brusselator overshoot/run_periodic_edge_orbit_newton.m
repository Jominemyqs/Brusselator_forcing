function results = run_periodic_edge_orbit_newton(cfg)
%RUN_PERIODIC_EDGE_ORBIT_NEWTON Converge the B--C edge candidate by shooting.
%   Solves Phi_T(X)-X=0 together with a fixed time-phase condition using a
%   matrix-free finite-difference Newton--GMRES method. The nonlinear and
%   Krylov norms use trapezoidal spatial weights in the original u/v units.

if nargin < 1 || isempty(cfg)
    cfg = brusselator_periodic_edge_newton_config();
end
outdir = fullfile(cfg.output.root,cfg.experiment_name);
checkpoint_file = fullfile(outdir,'newton_checkpoint.mat');
final_file = fullfile(outdir,'periodic_edge_orbit_newton.mat');
if isfile(final_file)
    error('run_periodic_edge_orbit_newton:OutputExists', ...
        'Refusing to overwrite completed result: %s',final_file);
end
if ~isfile(cfg.newton.seed_file)
    error('run_periodic_edge_orbit_newton:MissingSeed', ...
        'Missing fixed-phase shooting seed: %s',cfg.newton.seed_file);
end
loaded = load(cfg.newton.seed_file,'seed');
seed = loaded.seed;
reference = seed.reference_state;
F_reference = seed.reference_tangent;
x = reference(:,1);
N = numel(x);
w = spatial_weights(x);
sqrt_weights = sqrt([w;w]);
reference_vector = state_to_vector(reference);
reference_tangent_vector = fields_to_vector(F_reference);
reference_state_norm = norm(sqrt_weights.*reference_vector);
reference_tangent_norm = norm(sqrt_weights.*reference_tangent_vector);
period_scale = seed.period;
flow_evaluation_count = 0;

if isfile(checkpoint_file)
    saved = load(checkpoint_file,'z','history','next_iteration', ...
        'flow_evaluation_count','status');
    z = saved.z; history = saved.history;
    next_iteration = saved.next_iteration;
    flow_evaluation_count = saved.flow_evaluation_count;
    status = saved.status;
    fprintf('Resuming Newton--Krylov at iteration %d.\n',next_iteration);
else
    if isfolder(outdir) && folder_has_contents(outdir)
        error('run_periodic_edge_orbit_newton:IncompleteOutput', ...
            'Existing output lacks a checkpoint: %s',outdir);
    end
    if ~isfolder(outdir)
        mkdir(outdir);
    end
    z = encode(seed.state,seed.period);
    iteration_template = empty_iteration_record();
    history = iteration_template([]);
    next_iteration = 1;
    status = 'newton_in_progress';
    brusselator_write_metadata(outdir,cfg, ...
        brusselator_run_metadata(cfg,mfilename,NaN));
    save_checkpoint();
end

timer = tic;
while next_iteration <= cfg.newton.maximum_iterations
    iteration = next_iteration;
    [R,metrics] = evaluate_residual(z);
    fprintf(['Newton %d: ||G||=%.6e, flow=%.6e, phase=%.3e, ', ...
        'T=%.10f.\n'],iteration,norm(R),metrics.flow_residual, ...
        metrics.phase_residual,metrics.period);
    if converged(metrics,R)
        status = 'newton_converged';
        break;
    end

    base_z = z;
    base_R = R;
    forcing_tolerance = min(cfg.newton.gmres_maximum_tolerance, ...
        max(cfg.newton.gmres_minimum_tolerance,0.5*sqrt(norm(R))));
    Afun = @(direction) jacobian_vector_product( ...
        direction,base_z,base_R);
    [step,gmres_flag,gmres_relative_residual,gmres_iterations] = gmres( ...
        Afun,-R,cfg.newton.gmres_restart,forcing_tolerance, ...
        cfg.newton.gmres_maximum_outer_iterations);
    if any(~isfinite(step))
        status = 'nonfinite_newton_step';
        save_checkpoint();
        break;
    end

    alpha = 1;
    accepted = false;
    trial_R = [];
    trial_metrics = struct();
    while alpha >= cfg.newton.line_search_minimum_alpha
        trial_z = z + alpha*step;
        if trial_z(end) <= 0
            alpha = alpha/2;
            continue;
        end
        [trial_R,trial_metrics] = evaluate_residual(trial_z);
        if norm(trial_R) <= (1-cfg.newton.line_search_sufficient_decrease*alpha)*norm(R)
            accepted = true;
            break;
        end
        alpha = alpha/2;
    end

    record = empty_iteration_record();
    record.iteration = iteration;
    record.period_before = metrics.period;
    record.residual_before = norm(R);
    record.flow_residual_before = metrics.flow_residual;
    record.phase_residual_before = metrics.phase_residual;
    record.gmres_tolerance = forcing_tolerance;
    record.gmres_flag = gmres_flag;
    record.gmres_relative_residual = gmres_relative_residual;
    record.gmres_iterations = total_gmres_iterations( ...
        gmres_iterations,cfg.newton.gmres_restart);
    record.step_norm = norm(step);
    record.line_search_alpha = alpha;
    record.accepted = accepted;
    record.flow_evaluations = flow_evaluation_count;
    if accepted
        record.period_after = trial_metrics.period;
        record.residual_after = norm(trial_R);
        record.flow_residual_after = trial_metrics.flow_residual;
        record.phase_residual_after = trial_metrics.phase_residual;
        z = trial_z;
        next_iteration = iteration+1;
    else
        record.period_after = metrics.period;
        record.residual_after = norm(R);
        record.flow_residual_after = metrics.flow_residual;
        record.phase_residual_after = metrics.phase_residual;
        status = 'line_search_failed';
    end
    history(end+1,1) = record; %#ok<AGROW>
    write_progress();
    save_checkpoint();
    if ~accepted
        break;
    end
end

[final_R,final_metrics,final_state,final_flow_state] = evaluate_residual(z);
if converged(final_metrics,final_R)
    status = 'newton_converged';
elseif strcmp(status,'newton_in_progress') && ...
        next_iteration > cfg.newton.maximum_iterations
    status = 'maximum_newton_iterations_reached';
end

period = final_metrics.period;
orbit_cfg = cfg;
orbit_cfg.solver.output_dt = cfg.newton.orbit_output_dt;
par = brusselator_make_parameters(orbit_cfg,@(t)cfg.edge_tracking.frozen_b);
[~,orbit_times,orbit_V,orbit_U] = solve_brusselator_1d_forced( ...
    final_state,par,period,0);
verification = independently_verify(final_state,period,cfg);
summary = struct('status',status,'newton_iterations',numel(history), ...
    'flow_evaluations',flow_evaluation_count,'period',period, ...
    'normalized_residual',norm(final_R), ...
    'flow_residual',final_metrics.flow_residual, ...
    'phase_residual',final_metrics.phase_residual, ...
    'state_correction_from_seed',state_distance(final_state,seed.state,cfg.grid.Lx), ...
    'period_correction_from_seed',period-seed.period, ...
    'ode45_verification_flow_residual',verification.flow_residual, ...
    'claim_scope',sprintf(['a converged result verifies a periodic orbit of ', ...
        'the N=%d semidiscrete frozen system; it does not yet give Floquet ', ...
        'stability or continuum convergence'],N));
write_progress();
writetable(struct2table(summary,'AsArray',true), ...
    fullfile(outdir,'newton_summary.csv'));
results = struct('configuration',cfg,'seed',seed,'status',status, ...
    'state',final_state,'period',period,'flow_state',final_flow_state, ...
    'orbit_times',orbit_times,'orbit_U',orbit_U,'orbit_V',orbit_V, ...
    'history',history,'final_residual',final_R, ...
    'final_metrics',final_metrics,'verification',verification,'summary',summary);
save(final_file,'results','-v7.3');
brusselator_write_metadata(outdir,cfg, ...
    brusselator_run_metadata(cfg,mfilename,toc(timer)));
disp(struct2table(history,'AsArray',true));
disp(struct2table(summary,'AsArray',true));
fprintf('Periodic edge-orbit Newton result saved in: %s\n',outdir);

    function z_value = encode(state,period_value)
        y = state_to_vector(state);
        z_value = [sqrt_weights.*y;period_value/period_scale];
    end

    function [state,period_value] = decode(z_value)
        y = z_value(1:2*N)./sqrt_weights;
        state = vector_to_state(y,x);
        period_value = period_scale*z_value(end);
    end

    function [residual,metrics,state,flow_state] = evaluate_residual(z_value)
        [state,period_value] = decode(z_value);
        flow_state = frozen_flow(state,period_value,cfg,cfg.solver);
        flow_evaluation_count = flow_evaluation_count+1;
        y = state_to_vector(state);
        y_flow = state_to_vector(flow_state);
        state_residual = sqrt_weights.*(y_flow-y)/reference_state_norm;
        phase = sum(w.*((state(:,2)-reference(:,2)).*F_reference(:,1)+ ...
            (state(:,3)-reference(:,3)).*F_reference(:,2))) / ...
            (reference_state_norm*reference_tangent_norm);
        residual = [state_residual;phase];
        metrics = struct('flow_residual',norm(state_residual), ...
            'phase_residual',phase,'period',period_value);
    end

    function product = jacobian_vector_product(direction,base_z,base_R)
        direction_norm = norm(direction);
        if direction_norm == 0
            product = zeros(size(direction));
            return;
        end
        h = cfg.newton.finite_difference_relative_step * ...
            (1+norm(base_z))/direction_norm;
        product = (evaluate_residual(base_z+h*direction)-base_R)/h;
    end

    function yes = converged(metrics,residual)
        yes = norm(residual) <= cfg.newton.residual_tolerance && ...
            abs(metrics.phase_residual) <= cfg.newton.phase_tolerance;
    end

    function save_checkpoint()
        save(checkpoint_file,'cfg','z','history','next_iteration', ...
            'flow_evaluation_count','status','-v7.3');
    end

    function write_progress()
        if ~isempty(history)
            writetable(struct2table(history,'AsArray',true), ...
                fullfile(outdir,'newton_progress.csv'));
        end
    end
end

function state = frozen_flow(state,period,cfg,solver)
flow_cfg = cfg;
flow_cfg.solver = solver;
flow_cfg.solver.output_dt = period;
par = brusselator_make_parameters(flow_cfg,@(t)cfg.edge_tracking.frozen_b);
state = solve_brusselator_1d_forced(state,par,period,0);
end

function verification = independently_verify(state,period,cfg)
final_state = frozen_flow(state,period,cfg,cfg.newton.verification_solver);
verification = struct('solver',cfg.newton.verification_solver, ...
    'flow_residual',state_distance(state,final_state,cfg.grid.Lx), ...
    'final_state',final_state, ...
    'interpretation','independent tight ode45 integration of one period');
end

function y = state_to_vector(state)
y = [state(:,2);state(:,3)];
end

function y = fields_to_vector(fields)
y = [fields(:,1);fields(:,2)];
end

function state = vector_to_state(y,x)
N = numel(x);
state = [x,y(1:N),y(N+1:end)];
end

function weights = spatial_weights(x)
dx = x(2)-x(1);
weights = dx*ones(numel(x),1);
weights([1,end]) = 0.5*dx;
end

function distance = state_distance(A,B,Lx)
N = size(A,1); dx = Lx/(N-1); w = dx*ones(N,1); w([1,end]) = 0.5*dx;
numerator = sum(((A(:,2)-B(:,2)).^2+(A(:,3)-B(:,3)).^2).*w);
energy_A = sum((A(:,2).^2+A(:,3).^2).*w);
energy_B = sum((B(:,2).^2+B(:,3).^2).*w);
distance = sqrt(numerator/max(0.5*(energy_A+energy_B),eps));
end

function total = total_gmres_iterations(iterations,restart)
if isempty(iterations)
    total = 0;
elseif isscalar(iterations)
    total = iterations;
else
    total = max(0,iterations(1)-1)*restart+iterations(2);
end
end

function record = empty_iteration_record()
record = struct('iteration',NaN,'period_before',NaN,'residual_before',NaN, ...
    'flow_residual_before',NaN,'phase_residual_before',NaN, ...
    'gmres_tolerance',NaN,'gmres_flag',NaN,'gmres_relative_residual',NaN, ...
    'gmres_iterations',NaN,'step_norm',NaN,'line_search_alpha',NaN, ...
    'accepted',false,'period_after',NaN,'residual_after',NaN, ...
    'flow_residual_after',NaN,'phase_residual_after',NaN, ...
    'flow_evaluations',NaN);
end

function yes = folder_has_contents(folder)
entries = dir(folder);
names = {entries.name};
yes = any(~ismember(names,{'.','..'}));
end
