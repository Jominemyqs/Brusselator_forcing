function analysis = brusselator_periodic_recurrence(S, U, V, analysis_start, threshold)
%BRUSSELATOR_PERIODIC_RECURRENCE Test late full-state Poincare recurrence.
%   ANALYSIS = BRUSSELATOR_PERIODIC_RECURRENCE(S,U,V,T0,TOL) constructs an
%   upward Poincare section from the spatial mean of v over S >= T0. It
%   compares full (u,v) section states at several return lags, both directly
%   and after the exact Neumann-domain reflection x -> L-x.
%
%   This function diagnoses finite-time recurrence. It does not establish
%   orbital stability or prove that the recurrent set is an attractor.

validateattributes(S, {'numeric'}, {'vector', 'real', 'finite'});
validateattributes(U, {'numeric'}, {'2d', 'real', 'finite'});
validateattributes(V, {'numeric'}, {'2d', 'real', 'finite'});
validateattributes(analysis_start, {'numeric'}, {'scalar', 'real', 'finite'});
validateattributes(threshold, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});

S = S(:);
if size(U,1) ~= numel(S) || size(V,1) ~= numel(S) || ...
        size(U,2) ~= size(V,2)
    error('brusselator_periodic_recurrence:DimensionMismatch', ...
        'S, U, and V dimensions are inconsistent.');
end

late = S >= analysis_start;
if nnz(late) < 4
    error('brusselator_periodic_recurrence:InsufficientLateSamples', ...
        'At least four samples are required at or after analysis_start.');
end

S_late = S(late);
U_late = U(late,:);
V_late = V(late,:);
observable = mean(V_late, 2);
observable = observable - mean(observable);
cross_index = find(observable(1:end-1) <= 0 & observable(2:end) > 0);

times = zeros(numel(cross_index), 1);
states = zeros(numel(cross_index), 2 * size(U_late,2));
for k = 1:numel(cross_index)
    j = cross_index(k);
    denominator = observable(j+1) - observable(j);
    if abs(denominator) <= eps
        alpha = 0;
    else
        alpha = -observable(j) / denominator;
    end
    times(k) = S_late(j) + alpha * (S_late(j+1) - S_late(j));
    u = (1 - alpha) * U_late(j,:) + alpha * U_late(j+1,:);
    v = (1 - alpha) * V_late(j,:) + alpha * V_late(j+1,:);
    states(k,:) = [u, v];
end

section = struct('times', times, 'states', states, 'periods', [], ...
    'period', NaN, 'period_cv', NaN, 'crossing_count', numel(times));
record_template = struct('return_crossings', NaN, 'return_time', NaN, ...
    'direct_median_rel_diff_full_state', NaN, ...
    'reflection_median_rel_diff_full_state', NaN, ...
    'best_median_rel_diff_full_state', NaN, 'best_uses_reflection', false, ...
    'direct_max_rel_diff_full_state', NaN, ...
    'reflection_max_rel_diff_full_state', NaN);
recurrence = struct('records', record_template([]), ...
    'best', record_template, 'first_close', record_template);

if numel(times) < 4
    analysis = struct('section', section, 'recurrence', recurrence, ...
        'status', 'insufficient_section_crossings');
    return;
end

periods = diff(times);
section.periods = periods;
section.period = median(periods);
section.period_cv = std(periods) / max(mean(periods), eps);

max_lag = min(30, size(states,1) - 2);
records = repmat(record_template, max_lag, 1);
for lag = 1:max_lag
    first = states(1:end-lag,:);
    second = states(1+lag:end,:);
    direct = rowwise_relative_distance(second, first);
    reflected = rowwise_relative_distance(reflect_full_states(second), first);
    record = record_template;
    record.return_crossings = lag;
    record.return_time = median(times(1+lag:end) - times(1:end-lag));
    record.direct_median_rel_diff_full_state = median(direct);
    record.reflection_median_rel_diff_full_state = median(reflected);
    if record.reflection_median_rel_diff_full_state < ...
            record.direct_median_rel_diff_full_state
        record.best_median_rel_diff_full_state = ...
            record.reflection_median_rel_diff_full_state;
        record.best_uses_reflection = true;
    else
        record.best_median_rel_diff_full_state = ...
            record.direct_median_rel_diff_full_state;
    end
    record.direct_max_rel_diff_full_state = max(direct);
    record.reflection_max_rel_diff_full_state = max(reflected);
    records(lag) = record;
end

[~, best_index] = min([records.best_median_rel_diff_full_state]);
best = records(best_index);
close_index = find([records.best_median_rel_diff_full_state] < threshold, 1);
if isempty(close_index)
    first_close = record_template;
else
    first_close = records(close_index);
end
recurrence = struct('records', records, 'best', best, ...
    'first_close', first_close);
if isfinite(first_close.return_crossings) && section.period_cv < 0.05
    status = 'close_recurrence_with_regular_section_timing';
elseif section.period_cv < 0.02
    status = 'regular_section_timing_without_close_recurrence';
else
    status = 'no_close_regular_recurrence';
end
analysis = struct('section', section, 'recurrence', recurrence, 'status', status);
end

function distances = rowwise_relative_distance(A, reference)
numerator = sqrt(sum((A - reference).^2, 2));
denominator = sqrt(sum(reference.^2, 2));
distances = numerator ./ max(denominator, eps);
end

function reflected = reflect_full_states(states)
N = size(states,2) / 2;
if N ~= floor(N)
    error('brusselator_periodic_recurrence:InvalidStateDimension', ...
        'Full states must concatenate equally sized u and v fields.');
end
reflected = [fliplr(states(:,1:N)), fliplr(states(:,N+1:end))];
end
