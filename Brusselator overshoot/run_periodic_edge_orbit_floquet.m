function results = run_periodic_edge_orbit_floquet(cfg)
%RUN_PERIODIC_EDGE_ORBIT_FLOQUET Compute leading Floquet multipliers.
%   Applies the monodromy operator matrix-free by integrating the exact
%   variational equation along the Newton-converged periodic orbit. Arnoldi
%   vectors use the physical trapezoidal spatial weighting.

if nargin < 1 || isempty(cfg)
    cfg = brusselator_periodic_edge_floquet_config();
end
outdir = fullfile(cfg.output.root,cfg.experiment_name);
final_file = fullfile(outdir,'periodic_edge_orbit_floquet.mat');
if isfile(final_file)
    error('run_periodic_edge_orbit_floquet:OutputExists', ...
        'Refusing to overwrite completed result: %s',final_file);
end
if ~isfile(cfg.floquet.orbit_file)
    error('run_periodic_edge_orbit_floquet:MissingOrbit', ...
        'Missing Newton-converged orbit: %s',cfg.floquet.orbit_file);
end
loaded = load(cfg.floquet.orbit_file,'results');
orbit_result = loaded.results;
if ~strcmp(orbit_result.status,'newton_converged')
    error('run_periodic_edge_orbit_floquet:UnconvergedOrbit', ...
        'Floquet analysis requires a Newton-converged orbit.');
end
if isfolder(outdir)
    error('run_periodic_edge_orbit_floquet:OutputDirectoryExists', ...
        'Refusing to reuse output directory: %s',outdir);
end
mkdir(outdir);
timer = tic;
state = orbit_result.state;
period = orbit_result.period;
x = state(:,1); N = numel(x);
w = spatial_weights(x);
sqrt_weights = sqrt([w;w]);
y0 = state_to_vector(state);
[D2,d1,d2,a,b] = frozen_operators(state,cfg);
solver_options = odeset('RelTol',cfg.floquet.solver.RelTol, ...
    'AbsTol',cfg.floquet.solver.AbsTol,'MaxStep',cfg.floquet.solver.MaxStep, ...
    'Jacobian',@base_jacobian);
base_solution = ode15s(@base_rhs,[0,period],y0,solver_options);
base_closure_error = physical_relative_vector_distance( ...
    deval(base_solution,period),y0,w);
monodromy_actions = 0;

F0 = base_rhs(0,y0);
phase_direction = sqrt_weights.*F0;
phase_direction = phase_direction/norm(phase_direction);
phase_image = monodromy_action(phase_direction);
phase_rayleigh_multiplier = phase_direction'*phase_image;
phase_direction_residual = norm(phase_image- ...
    phase_rayleigh_multiplier*phase_direction);

options = struct('tol',cfg.floquet.arnoldi_tolerance, ...
    'maxit',cfg.floquet.arnoldi_maximum_iterations, ...
    'p',cfg.floquet.arnoldi_subspace_dimension,'issym',false, ...
    'isreal',true,'v0',phase_direction,'disp',0);
[vectors,D,arnoldi_flag] = eigs(@monodromy_action,2*N, ...
    cfg.floquet.number_multipliers,'largestabs',options);
multipliers = diag(D);
[~,order] = sort(abs(multipliers),'descend');
multipliers = multipliers(order);
vectors = vectors(:,order);
eigen_residuals = zeros(numel(multipliers),1);
phase_alignments = zeros(numel(multipliers),1);
for k = 1:numel(multipliers)
    image = monodromy_action(vectors(:,k));
    eigen_residuals(k) = norm(image-multipliers(k)*vectors(:,k)) / ...
        max(norm(vectors(:,k)),eps);
    phase_alignments(k) = abs(vectors(:,k)'*phase_direction) / ...
        max(norm(vectors(:,k)),eps);
end
[~,phase_index] = min(abs(multipliers-cfg.floquet.phase_multiplier_target));
nontrivial = true(size(multipliers)); nontrivial(phase_index) = false;
unstable_indices = find(nontrivial & ...
    abs(multipliers)>cfg.floquet.unstable_modulus_threshold);

leading_index = NaN;
if ~isempty(unstable_indices)
    [~,j] = max(abs(multipliers(unstable_indices)));
    leading_index = unstable_indices(j);
end
finite_difference_check = struct('performed',false);
if isfinite(leading_index) && ...
        abs(imag(multipliers(leading_index))) < 1e-8
    eigenvector = real(vectors(:,leading_index));
    eigenvector = eigenvector/norm(eigenvector);
    finite_difference_check = check_tangent_action(eigenvector, ...
        monodromy_action(eigenvector));
end

records = multiplier_records(multipliers,eigen_residuals,phase_alignments, ...
    phase_index,period,cfg);
summary = struct('status','floquet_computation_complete', ...
    'period',period,'base_closure_error',base_closure_error, ...
    'arnoldi_flag',arnoldi_flag,'monodromy_actions',monodromy_actions, ...
    'phase_multiplier_real',real(multipliers(phase_index)), ...
    'phase_multiplier_imaginary',imag(multipliers(phase_index)), ...
    'phase_multiplier_residual',eigen_residuals(phase_index), ...
    'phase_direction_rayleigh_multiplier',phase_rayleigh_multiplier, ...
    'phase_direction_residual',phase_direction_residual, ...
    'nontrivial_unstable_count',numel(unstable_indices), ...
    'leading_unstable_multiplier_real',conditional_component( ...
        multipliers,leading_index,'real'), ...
    'leading_unstable_multiplier_imaginary',conditional_component( ...
        multipliers,leading_index,'imag'), ...
    'predicted_multiplier',cfg.floquet.predicted_transverse_multiplier, ...
    'leading_prediction_difference',conditional_difference( ...
        multipliers,leading_index,cfg.floquet.predicted_transverse_multiplier), ...
    'interpretation',sprintf(['leading Floquet multipliers of the N=%d ', ...
        'semidiscrete periodic orbit; one trivial unit multiplier is ', ...
        'required by time phase'],N));
writetable(struct2table(records,'AsArray',true), ...
    fullfile(outdir,'floquet_multipliers.csv'));
writetable(struct2table(summary,'AsArray',true), ...
    fullfile(outdir,'floquet_summary.csv'));
make_figure(outdir,multipliers,phase_index,unstable_indices);
results = struct('configuration',cfg,'orbit_file',cfg.floquet.orbit_file, ...
    'state',state,'period',period,'multipliers',multipliers, ...
    'weighted_eigenvectors',vectors,'eigen_residuals',eigen_residuals, ...
    'phase_alignments',phase_alignments,'phase_index',phase_index, ...
    'unstable_indices',unstable_indices,'leading_unstable_index',leading_index, ...
    'phase_direction',phase_direction, ...
    'finite_difference_check',finite_difference_check,'records',records, ...
    'summary',summary);
save(final_file,'results','-v7.3');
brusselator_write_metadata(outdir,cfg, ...
    brusselator_run_metadata(cfg,mfilename,toc(timer)));
disp(struct2table(records,'AsArray',true));
disp(struct2table(summary,'AsArray',true));
fprintf('Floquet result saved in: %s\n',outdir);

    function f = base_rhs(~,y)
        u = y(1:N); v = y(N+1:end);
        f = [d1*D2*u+a-(b+1)*u+u.^2.*v; ...
             d2*D2*v+b*u-u.^2.*v];
    end

    function J = base_jacobian(~,y)
        u = y(1:N); v = y(N+1:end);
        J11 = d1*D2+spdiags(-(b+1)+2*u.*v,0,N,N);
        J12 = spdiags(u.^2,0,N,N);
        J21 = spdiags(b-2*u.*v,0,N,N);
        J22 = d2*D2-spdiags(u.^2,0,N,N);
        J = [J11,J12;J21,J22];
    end

    function q_final = monodromy_action(q_initial)
        if ~isreal(q_initial)
            q_final = monodromy_real_action(real(q_initial)) + ...
                1i*monodromy_real_action(imag(q_initial));
        else
            q_final = monodromy_real_action(q_initial);
        end
    end

    function q_final = monodromy_real_action(q_initial)
        perturbation = q_initial./sqrt_weights;
        tangent_options = odeset('RelTol',cfg.floquet.solver.RelTol, ...
            'AbsTol',cfg.floquet.solver.AbsTol, ...
            'MaxStep',cfg.floquet.solver.MaxStep, ...
            'Jacobian',@tangent_jacobian);
        [~,solution] = ode15s(@tangent_rhs,[0,period],perturbation, ...
            tangent_options);
        q_final = sqrt_weights.*solution(end,:)';
        monodromy_actions = monodromy_actions+1;
    end

    function f = tangent_rhs(t,perturbation)
        f = tangent_matrix(t)*perturbation;
    end

    function J = tangent_jacobian(t,~)
        J = tangent_matrix(t);
    end

    function J = tangent_matrix(t)
        J = base_jacobian(t,deval(base_solution,t));
    end

    function check = check_tangent_action(q_direction,q_tangent)
        amplitude = cfg.floquet.finite_difference_check_amplitude;
        state_norm = norm(sqrt_weights.*y0);
        delta_q = amplitude*state_norm*q_direction;
        perturbed_y = y0+delta_q./sqrt_weights;
        perturbed_state = vector_to_state(perturbed_y,x);
        flow_cfg = cfg;
        flow_cfg.solver = struct('name','ode15s','output_dt',period, ...
            'RelTol',cfg.floquet.solver.RelTol, ...
            'AbsTol',cfg.floquet.solver.AbsTol, ...
            'MaxStep',cfg.floquet.solver.MaxStep);
        par = brusselator_make_parameters(flow_cfg,@(t)b);
        perturbed_final = solve_brusselator_1d_forced( ...
            perturbed_state,par,period,0);
        finite_difference_image = sqrt_weights.*( ...
            state_to_vector(perturbed_final)-deval(base_solution,period)) / ...
            (amplitude*state_norm);
        check = struct('performed',true,'relative_amplitude',amplitude, ...
            'relative_error',norm(finite_difference_image-q_tangent) / ...
                max(norm(q_tangent),eps), ...
            'finite_difference_image_norm',norm(finite_difference_image), ...
            'tangent_image_norm',norm(q_tangent));
    end
end

function [D2,d1,d2,a,b] = frozen_operators(state,cfg)
N = size(state,1); dx = state(2,1)-state(1,1); e = ones(N,1);
D2 = sparse(1:N-1,[2:N-1 N],ones(N-1,1),N,N) ...
    - sparse(1:N,1:N,e,N,N);
D2 = D2+D2'; D2(1,2)=2; D2(N,N-1)=2; D2=D2/dx^2;
d2 = cfg.model.d2; d1 = cfg.model.sigma*d2; a = cfg.model.a;
b = cfg.edge_tracking.frozen_b;
end

function records = multiplier_records(multipliers,residuals,alignments, ...
        phase_index,period,cfg)
template = struct('rank',NaN,'real_part',NaN,'imaginary_part',NaN, ...
    'modulus',NaN,'floquet_exponent_real',NaN,'eigen_residual',NaN, ...
    'phase_alignment',NaN,'is_phase_multiplier',false, ...
    'is_nontrivially_unstable',false);
records = repmat(template,numel(multipliers),1);
for k = 1:numel(multipliers)
    records(k).rank = k;
    records(k).real_part = real(multipliers(k));
    records(k).imaginary_part = imag(multipliers(k));
    records(k).modulus = abs(multipliers(k));
    records(k).floquet_exponent_real = log(abs(multipliers(k))) / period;
    records(k).eigen_residual = residuals(k);
    records(k).phase_alignment = alignments(k);
    records(k).is_phase_multiplier = k==phase_index;
    records(k).is_nontrivially_unstable = k~=phase_index && ...
        abs(multipliers(k))>cfg.floquet.unstable_modulus_threshold;
end
end

function value = conditional_component(values,index,component)
if ~isfinite(index)
    value = NaN;
elseif strcmp(component,'real')
    value = real(values(index));
else
    value = imag(values(index));
end
end

function value = conditional_difference(values,index,target)
if ~isfinite(index)
    value = NaN;
else
    value = abs(values(index)-target);
end
end

function y = state_to_vector(state)
y = [state(:,2);state(:,3)];
end

function state = vector_to_state(y,x)
N = numel(x); state = [x,y(1:N),y(N+1:end)];
end

function weights = spatial_weights(x)
dx = x(2)-x(1); weights = dx*ones(numel(x),1); weights([1,end])=0.5*dx;
end

function distance = physical_relative_vector_distance(A,B,w)
N = numel(w);
numerator = sum(((A(1:N)-B(1:N)).^2+(A(N+1:end)-B(N+1:end)).^2).*w);
energy_A = sum((A(1:N).^2+A(N+1:end).^2).*w);
energy_B = sum((B(1:N).^2+B(N+1:end).^2).*w);
distance = sqrt(numerator/max(0.5*(energy_A+energy_B),eps));
end

function make_figure(outdir,multipliers,phase_index,unstable_indices)
fig = figure('Color','w','Position',[100 100 750 700]); ax = axes(fig); hold(ax,'on');
theta = linspace(0,2*pi,500); plot(ax,cos(theta),sin(theta),'k--','LineWidth',1);
scatter(ax,real(multipliers),imag(multipliers),55,[0.2 0.5 0.8],'filled');
scatter(ax,real(multipliers(phase_index)),imag(multipliers(phase_index)), ...
    90,[0.25 0.65 0.35],'filled','DisplayName','phase');
if ~isempty(unstable_indices)
    scatter(ax,real(multipliers(unstable_indices)),imag(multipliers(unstable_indices)), ...
        95,[0.85 0.35 0.15],'filled','DisplayName','unstable');
end
xline(ax,0,'Color',[0.7 0.7 0.7]); yline(ax,0,'Color',[0.7 0.7 0.7]);
axis(ax,'equal'); grid(ax,'on'); box(ax,'on');
xlabel(ax,'Re(\mu)'); ylabel(ax,'Im(\mu)'); title(ax,'Leading Floquet multipliers');
exportgraphics(fig,fullfile(outdir,'floquet_multipliers.png'),'Resolution',300);
close(fig);
end
