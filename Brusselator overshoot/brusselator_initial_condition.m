function [x, ic0] = brusselator_initial_condition(cfg)
%BRUSSELATOR_INITIAL_CONDITION Build the seeded perturbation about (a,b0/a).
%   The caller's random-number state is restored before return so that a
%   seed is explicit in the experiment configuration rather than implicit
%   in execution order.

brusselator_validate_config(cfg);

previous_rng = rng;
restore_rng = onCleanup(@() rng(previous_rng)); %#ok<NASGU>
rng(cfg.initial.seed, 'twister');

x = linspace(0, cfg.grid.Lx, cfg.grid.N)';
u_ss = cfg.model.a;
v_ss = cfg.forcing.b0 / cfg.model.a;
amplitude = cfg.initial.perturbation_amplitude;

u0 = u_ss + amplitude * randn(cfg.grid.N, 1);
v0 = v_ss + amplitude * randn(cfg.grid.N, 1);
ic0 = [x, u0, v0];
end
