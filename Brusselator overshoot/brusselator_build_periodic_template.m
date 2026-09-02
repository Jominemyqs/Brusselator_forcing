function template = brusselator_build_periodic_template(S, U, V, period, phase_samples)
%BRUSSELATOR_BUILD_PERIODIC_TEMPLATE Sample the latest complete orbit cycle.
%   The returned rows use a normalized temporal phase in [0,1). The final
%   endpoint is excluded so circular phase shifts do not duplicate phase zero.

if nargin < 5 || isempty(phase_samples)
    phase_samples = 360;
end
S = S(:);
validateattributes(period, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(phase_samples, {'numeric'}, ...
    {'scalar','integer','>=',16});
if size(U,1) ~= numel(S) || size(V,1) ~= numel(S) || ...
        size(U,2) ~= size(V,2)
    error('brusselator_build_periodic_template:DimensionMismatch', ...
        'S, U, and V dimensions are inconsistent.');
end
start_time = S(end) - period;
if start_time < S(1)
    error('brusselator_build_periodic_template:InsufficientWindow', ...
        'The trajectory does not contain one complete requested period.');
end
phase_fraction = (0:phase_samples-1)' / phase_samples;
query_times = start_time + phase_fraction * period;
template = struct('phase_fraction', phase_fraction, ...
    'phase_times', phase_fraction * period, ...
    'period', period, 'start_time', start_time, ...
    'U', interp1(S, U, query_times, 'linear'), ...
    'V', interp1(S, V, query_times, 'linear'));
end
