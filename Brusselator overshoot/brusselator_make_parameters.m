function par = brusselator_make_parameters(cfg, Bfun)
%BRUSSELATOR_MAKE_PARAMETERS Translate a run config into solver inputs.

if ~isa(Bfun, 'function_handle')
    error('brusselator_make_parameters:InvalidForcing', ...
        'Bfun must be a function handle.');
end

par = struct();
par.sigma = cfg.model.sigma;
par.a = cfg.model.a;
par.d2 = cfg.model.d2;
par.Bfun = Bfun;
par.solver = cfg.solver;
end
