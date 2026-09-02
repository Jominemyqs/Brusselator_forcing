function result = brusselator_run_pair(cfg)
%BRUSSELATOR_RUN_PAIR Run one matched control/forced experiment.
%   Both simulations receive exactly the same seeded initial condition.

brusselator_validate_config(cfg);
[x, ic0] = brusselator_initial_condition(cfg);

control_par = brusselator_make_parameters(cfg, @(t) cfg.forcing.b0); %#ok<NASGU>
forced_par = brusselator_make_parameters(cfg, @(t) overshoot_B(t, ...
    cfg.forcing.b0, cfg.forcing.bmax, cfg.forcing.Tup, ...
    cfg.forcing.Thold, cfg.forcing.Tdown));

run_timer = tic;
[ic_control, S_control, V_control] = solve_brusselator_1d_forced( ...
    ic0, control_par, cfg.time.Tfinal, 0);
[ic_forced, S_forced, V_forced] = solve_brusselator_1d_forced( ...
    ic0, forced_par, cfg.time.Tfinal, 0);
runtime_seconds = toc(run_timer);

u_control = ic_control(:,2);
v_control = ic_control(:,3);
u_forced = ic_forced(:,2);
v_forced = ic_forced(:,3);

metrics.rel_diff_u = relative_l2(u_forced, u_control);
metrics.rel_diff_v = relative_l2(v_forced, v_control);
metrics.max_diff_u = max(abs(u_forced - u_control));
metrics.max_diff_v = max(abs(v_forced - v_control));

if metrics.rel_diff_v < cfg.classification.recovery_threshold
    class_label = 'Recovery';
elseif metrics.rel_diff_v > cfg.classification.transition_threshold
    class_label = 'LastingTransition';
else
    class_label = 'Intermediate';
end

result = struct();
result.schema_version = cfg.schema_version;
result.config = cfg;
result.x = x;
result.ic0 = ic0;
result.control = struct('S', S_control, 'V', V_control, 'ic_final', ic_control);
result.forced = struct('S', S_forced, 'V', V_forced, 'ic_final', ic_forced);
result.metrics = metrics;
result.class_label = class_label;
result.runtime_seconds = runtime_seconds;
end

function value = relative_l2(candidate, reference)
value = norm(candidate - reference, 2) / max(norm(reference, 2), 1e-12);
end
