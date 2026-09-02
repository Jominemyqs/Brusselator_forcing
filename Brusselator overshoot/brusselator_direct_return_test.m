function analysis = brusselator_direct_return_test(section_state, cfg, ...
        frozen_b, period_guess, threshold)
%BRUSSELATOR_DIRECT_RETURN_TEST Verify a period by a full-state restart.
%   Restarts a phase-fixed full state under the frozen system, samples the
%   evolution finely, and locates the first close identity return. Reflection
%   returns are recorded separately. The resulting period is based on direct
%   full-state recurrence, not solely on section-crossing times.

section_state = section_state(:)';
N = numel(section_state) / 2;
if N ~= floor(N) || N ~= cfg.grid.N
    error('brusselator_direct_return_test:DimensionMismatch', ...
        'The section state must contain concatenated u and v fields on cfg.grid.N.');
end
validateattributes(period_guess, {'numeric'}, ...
    {'scalar','real','finite','positive'});
validateattributes(threshold, {'numeric'}, ...
    {'scalar','real','finite','positive'});

test_cfg = cfg;
test_cfg.solver.output_dt = min(0.01, cfg.solver.output_dt);
if isempty(test_cfg.solver.MaxStep)
    test_cfg.solver.MaxStep = 0.02;
else
    test_cfg.solver.MaxStep = min(0.02, test_cfg.solver.MaxStep);
end
duration = 2.25 * period_guess;
x = linspace(0, cfg.grid.Lx, N)';
initial_state = [x, section_state(1:N)', section_state(N+1:end)'];
par = brusselator_make_parameters(test_cfg, @(t) frozen_b);
[~, S, V, U] = solve_brusselator_1d_forced(initial_state, par, duration, 0);

u0 = section_state(1:N);
v0 = section_state(N+1:end);
normalizer = max(sqrt(sum(u0.^2) + sum(v0.^2)), eps);
direct = sqrt(sum((U - u0).^2 + (V - v0).^2, 2)) / normalizer;
reflected = sqrt(sum((fliplr(U) - u0).^2 + ...
    (fliplr(V) - v0).^2, 2)) / normalizer;
eligible = S >= 0.35 * period_guess;
[direct_period, direct_distance, direct_close] = ...
    first_close_local_minimum(S, direct, eligible, threshold);
[reflection_period, reflection_distance, reflection_close] = ...
    first_close_local_minimum(S, reflected, eligible, threshold);

analysis = struct('period_guess', period_guess, 'duration', duration, ...
    'output_dt', test_cfg.solver.output_dt, 'MaxStep', test_cfg.solver.MaxStep, ...
    'direct_period', direct_period, 'direct_distance', direct_distance, ...
    'direct_close_return_found', direct_close, ...
    'reflection_period', reflection_period, ...
    'reflection_distance', reflection_distance, ...
    'reflection_close_return_found', reflection_close, ...
    'S', S, 'direct_distance_history', direct, ...
    'reflection_distance_history', reflected, ...
    'metric', ['full-state relative L2 return from a phase-fixed restart; ', ...
        'identity and exact reflection reported separately']);
end

function [period, distance, close_found] = ...
        first_close_local_minimum(S, values, eligible, threshold)
local = false(size(values));
local(2:end-1) = values(2:end-1) <= values(1:end-2) & ...
    values(2:end-1) <= values(3:end);
candidates = find(local & eligible);
close_index = candidates(find(values(candidates) < threshold, 1));
if isempty(close_index)
    eligible_indices = find(eligible);
    [distance, local_index] = min(values(eligible_indices));
    index = eligible_indices(local_index);
    close_found = false;
else
    index = close_index;
    distance = values(index);
    close_found = true;
end
period = S(index);
end
