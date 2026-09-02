function analysis = brusselator_direct_period_recurrence(S, U, V, ...
        analysis_start, period_guess)
%BRUSSELATOR_DIRECT_PERIOD_RECURRENCE Estimate period from full trajectories.
%   Compares X(t+tau) directly with X(t) over a late-time window and searches
%   over tau near a supplied Poincare-period estimate. Identity and exact
%   spatial-reflection returns are optimized separately. This complements,
%   rather than replaces, section-crossing diagnostics.

S = S(:);
if size(U,1) ~= numel(S) || size(V,1) ~= numel(S) || ...
        size(U,2) ~= size(V,2)
    error('brusselator_direct_period_recurrence:DimensionMismatch', ...
        'S, U, and V dimensions are inconsistent.');
end
validateattributes(period_guess, {'numeric'}, ...
    {'scalar','real','finite','positive'});

% The section can return after a reflection-related half-cycle. Searching
% beyond twice the supplied estimate lets the identity return determine the
% direct full-state period instead of silently equating section and orbit
% periods.
lag_bounds = [0.35, 2.25] * period_guess;
window_start = max(analysis_start, S(end) - max(10 * period_guess, 40));
window_end = S(end) - lag_bounds(2);
if window_end <= window_start
    error('brusselator_direct_period_recurrence:InsufficientWindow', ...
        'The late trajectory is too short for direct recurrence analysis.');
end
sample_dt = max(median(diff(S)), 0.1);
base_times = (window_start:sample_dt:window_end)';
if numel(base_times) < 10
    base_times = linspace(window_start, window_end, 10)';
end
base_U = interp1(S, U, base_times, 'pchip');
base_V = interp1(S, V, base_times, 'pchip');

scan_lags = linspace(lag_bounds(1), lag_bounds(2), 381)';
direct_scan = zeros(size(scan_lags));
reflection_scan = zeros(size(scan_lags));
for k = 1:numel(scan_lags)
    [direct_scan(k), reflection_scan(k)] = lag_distances( ...
        scan_lags(k), S, U, V, base_times, base_U, base_V);
end
[direct_period, direct_distance] = refine_minimum( ...
    scan_lags, direct_scan, false, S, U, V, base_times, base_U, base_V);
[reflection_period, reflection_distance] = refine_minimum( ...
    scan_lags, reflection_scan, true, S, U, V, base_times, base_U, base_V);

if reflection_distance < direct_distance
    best_period = reflection_period;
    best_distance = reflection_distance;
    best_uses_reflection = true;
else
    best_period = direct_period;
    best_distance = direct_distance;
    best_uses_reflection = false;
end
analysis = struct( ...
    'analysis_start', analysis_start, ...
    'comparison_window_start', window_start, ...
    'comparison_window_end', window_end, ...
    'period_guess', period_guess, ...
    'lag_bounds', lag_bounds, ...
    'direct_period', direct_period, ...
    'direct_distance', direct_distance, ...
    'reflection_period', reflection_period, ...
    'reflection_distance', reflection_distance, ...
    'best_period', best_period, ...
    'best_distance', best_distance, ...
    'best_uses_reflection', best_uses_reflection, ...
    'scan_lags', scan_lags, ...
    'direct_scan', direct_scan, ...
    'reflection_scan', reflection_scan, ...
    'metric', ['symmetric relative RMS over full (u,v) trajectory window; ', ...
        'identity and x-to-L-minus-x reflection tested separately']);
end

function [period, distance] = refine_minimum(lags, values, reflect, ...
        S, U, V, base_times, base_U, base_V)
[~, index] = min(values);
left = lags(max(1, index - 1));
right = lags(min(numel(lags), index + 1));
if right <= left
    period = lags(index);
    distance = values(index);
    return;
end
objective = @(lag) selected_distance(lag, reflect, ...
    S, U, V, base_times, base_U, base_V);
[period, distance] = fminbnd(objective, left, right, ...
    optimset('TolX', 1e-7, 'Display', 'off'));
end

function value = selected_distance(lag, reflect, ...
        S, U, V, base_times, base_U, base_V)
[direct, reflected] = lag_distances(lag, S, U, V, ...
    base_times, base_U, base_V);
if reflect
    value = reflected;
else
    value = direct;
end
end

function [direct, reflected] = lag_distances(lag, S, U, V, ...
        base_times, base_U, base_V)
shifted_U = interp1(S, U, base_times + lag, 'pchip');
shifted_V = interp1(S, V, base_times + lag, 'pchip');
normalizer = sqrt(0.5 * mean(sum(base_U.^2 + base_V.^2, 2) + ...
    sum(shifted_U.^2 + shifted_V.^2, 2)));
direct = sqrt(mean(sum((shifted_U - base_U).^2 + ...
    (shifted_V - base_V).^2, 2))) / max(normalizer, eps);
reflected = sqrt(mean(sum((fliplr(shifted_U) - base_U).^2 + ...
    (fliplr(shifted_V) - base_V).^2, 2))) / max(normalizer, eps);
end
