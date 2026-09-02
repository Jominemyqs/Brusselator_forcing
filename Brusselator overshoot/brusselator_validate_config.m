function brusselator_validate_config(cfg)
%BRUSSELATOR_VALIDATE_CONFIG Check constraints needed by new experiments.

required_top = {'model', 'grid', 'forcing', 'time', 'solver', 'initial', ...
    'classification', 'output'};
for k = 1:numel(required_top)
    if ~isfield(cfg, required_top{k})
        error('brusselator_validate_config:MissingField', ...
            'Missing cfg.%s.', required_top{k});
    end
end

if cfg.grid.N < 3 || cfg.grid.N ~= floor(cfg.grid.N)
    error('brusselator_validate_config:InvalidGrid', ...
        'cfg.grid.N must be an integer of at least three.');
end
if ~(isscalar(cfg.grid.Lx) && isfinite(cfg.grid.Lx) && cfg.grid.Lx > 0)
    error('brusselator_validate_config:InvalidDomain', ...
        'cfg.grid.Lx must be positive.');
end
if ~(cfg.model.d2 > 0 && cfg.model.sigma > 0)
    error('brusselator_validate_config:InvalidDiffusion', ...
        'd2 and sigma must be positive.');
end

timing = [cfg.forcing.Tup, cfg.forcing.Thold, cfg.forcing.Tdown];
if any(~isfinite(timing)) || any(timing < 0) || cfg.forcing.Tup == 0 || cfg.forcing.Tdown == 0
    error('brusselator_validate_config:InvalidForcing', ...
        'Ramp times must be positive and the hold time must be nonnegative.');
end
if ~(isscalar(cfg.time.post_forcing) && isfinite(cfg.time.post_forcing) && ...
        cfg.time.post_forcing >= 0)
    error('brusselator_validate_config:InvalidPostForcingTime', ...
        'cfg.time.post_forcing must be nonnegative.');
end
if cfg.time.Tfinal < cfg.time.forcing_end
    error('brusselator_validate_config:ForcingNotComplete', ...
        'Tfinal must reach the end of the forcing schedule.');
end
if ~ismember(lower(cfg.solver.name), {'ode45', 'ode15s'})
    error('brusselator_validate_config:UnsupportedSolver', ...
        'cfg.solver.name must be ode45 or ode15s.');
end
if cfg.classification.recovery_threshold >= cfg.classification.transition_threshold
    error('brusselator_validate_config:InvalidClassification', ...
        'Recovery threshold must be less than transition threshold.');
end
end
