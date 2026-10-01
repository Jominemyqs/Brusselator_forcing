function [ic, S, V, U] = solve_brusselator_1d_forced(ic, par, T, flag)
%SOLVE_BRUSSELATOR_1D_FORCED Integrate the 1D forced Brusselator.
%   This backward-compatible method-of-lines solver uses a second-order
%   finite-difference Laplacian with homogeneous Neumann conditions.  The
%   legacy defaults are retained unless optional fields are supplied:
%
%     par.sigma, par.Bfun                         required
%     par.a, par.d2                               optional model values
%     par.solver.output_dt                        output spacing (default 0.1)
%     par.solver.RelTol, par.solver.AbsTol,
%     par.solver.MaxStep                          optional ODE tolerances
%
%   OUTPUTS
%     ic  final [x,u,v] state
%     S   requested output times
%     V   v-component at every requested output time
%     U   u-component at every requested output time (optional fourth output)
%
%   The optional solver fields allow convergence tests without silently
%   changing the historical calculations, which use ode45 defaults.

%% setup
x  = ic(:,1);
u  = [ic(:,2); ic(:,3)];
N  = length(x);
dx = x(2)-x(1);

if N < 3
    error('solve_brusselator_1d_forced:GridTooSmall', ...
        'At least three grid points are required.');
end

%% parameters
if ~isfield(par, 'sigma') || ~isfield(par, 'Bfun')
    error('solve_brusselator_1d_forced:MissingParameter', ...
        'par.sigma and par.Bfun are required.');
end

% De Wit parameterization. Optional fields make each run self-describing.
a = get_optional_field(par, 'a', 2.5);
d2 = get_optional_field(par, 'd2', 9.73);
sigma = par.sigma;
d1 = sigma * d2;

Bfun  = par.Bfun;

% % Tzou
% a   = 1.4;
% d1  = 0.2666;
% d2  = 1;

%% solver configuration
solver = get_optional_field(par, 'solver', struct());
solver_name = lower(get_optional_field(solver, 'name', 'ode45'));
output_dt = get_optional_field(solver, 'output_dt', 0.1);
if ~(isscalar(output_dt) && isfinite(output_dt) && output_dt > 0)
    error('solve_brusselator_1d_forced:InvalidOutputDt', ...
        'par.solver.output_dt must be a positive finite scalar.');
end

tspan = 0:output_dt:T;
if tspan(end) < T
    tspan = [tspan, T]; %#ok<AGROW>
end

ode_options = odeset();
ode_options = set_ode_option(ode_options, solver, 'RelTol');
ode_options = set_ode_option(ode_options, solver, 'AbsTol');
ode_options = set_ode_option(ode_options, solver, 'MaxStep');

switch solver_name
    case 'ode45'
        ode_solver = @ode45;
    case 'ode15s'
        ode_solver = @ode15s;
    otherwise
        error('solve_brusselator_1d_forced:UnsupportedSolver', ...
            'Supported solver names are ode45 and ode15s.');
end

%% compute finite difference approximation D2 of second derivative
e  = ones(N,1);
D2 = sparse(1:N-1,[2:N-1 N],ones(N-1,1),N,N) - sparse(1:N,[1:N],e,N,N);
D2 = D2 + D2';
D2(1,2)   = 2;
D2(N,N-1) = 2;   % Neumann boundary conditions
D2 = D2 / dx^2;

% Optional exact sparse Jacobian; the historical default remains unchanged.
if strcmp(solver_name,'ode15s') && get_optional_field(solver,'analytic_jacobian',false)
    ode_options = odeset(ode_options,'Jacobian',@jacobian);
end

%% define right-hand side of partial differential equation
function f = rhs(t,u)
    u1 = u(1:N);
    u2 = u(N+1:2*N);

    b = Bfun(t);

    f = [d1*D2*u1 + a - (b+1)*u1 + u1.^2.*u2; ...
         d2*D2*u2 + b*u1 - u1.^2.*u2];
end

function J = jacobian(t,y)
    u1 = y(1:N); u2 = y(N+1:2*N); b = Bfun(t);
    cross = 2*u1.*u2; square = u1.^2;
    J = [d1*D2 + spdiags(-(b+1)+cross,0,N,N), spdiags(square,0,N,N); ...
         spdiags(b-cross,0,N,N), d2*D2-spdiags(square,0,N,N)];
end

%% compute solution
[S, UV] = ode_solver(@rhs, tspan, u, ode_options);

V  = UV(:, N+1:end);
if nargout >= 4
    U = UV(:, 1:N);
end
ic = [x, UV(end,1:N)', UV(end,N+1:end)'];

%% space-time plot
if flag == 1
    figure(1); clf;
    image(V, 'CDataMapping', 'scaled');
    ax = gca;
    ax.YDir = 'normal';
%   colormap(ax, brewermap([],'RdYlBu'));
    colormap(ax, 'sky');
    colorbar(ax);
    drawnow;
end

end

function value = get_optional_field(s, field_name, default_value)
if isstruct(s) && isfield(s, field_name) && ~isempty(s.(field_name))
    value = s.(field_name);
else
    value = default_value;
end
end

function options = set_ode_option(options, solver, field_name)
if isfield(solver, field_name) && ~isempty(solver.(field_name))
    options = odeset(options, field_name, solver.(field_name));
end
end
