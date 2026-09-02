function [record, history] = brusselator_classify_three_way_trajectory(S, U, V, refs, options)
%BRUSSELATOR_CLASSIFY_THREE_WAY_TRAJECTORY Classify a full-state trajectory.
%   The late window is labeled A, B, C, ambiguous, or unresolved using the
%   same phase/reflection-aware distances and finite thresholds as the saved
%   three-way remap. Unresolved is an explicit scientific outcome, not an
%   error or forced nearest-class assignment.

arguments
    S (:,1) double
    U (:,:) double
    V (:,:) double
    refs (1,1) struct
    options.late_window (1,1) double {mustBePositive} = 20
    options.sample_interval (1,1) double {mustBePositive} = 1
    options.median_threshold (1,1) double {mustBePositive} = 1e-2
    options.maximum_threshold (1,1) double {mustBePositive} = 2e-2
end

if size(U,1) ~= numel(S) || size(V,1) ~= numel(S) || ...
        size(U,2) ~= size(V,2)
    error('brusselator_classify_three_way_trajectory:DimensionMismatch', ...
        'S, U, and V dimensions are inconsistent.');
end
if S(end) - S(1) < options.late_window
    error('brusselator_classify_three_way_trajectory:ShortTrajectory', ...
        'The trajectory is shorter than the requested late window.');
end

dt = median(diff(S));
late_start = S(end) - options.late_window;
eligible = find(S >= late_start);
stride = max(1, round(options.sample_interval / dt));
indices = unique([eligible(1:stride:end); eligible(end)]);
times = S(indices);
distances = zeros(numel(indices), 3);
for k = 1:numel(indices)
    distances(k,:) = brusselator_three_way_distances( ...
        U(indices(k),:), V(indices(k),:), refs);
end
medians = median(distances, 1);
maxima = max(distances, [], 1);
valid = medians < options.median_threshold & maxima < options.maximum_threshold;
outcome = outcome_label(valid);

record = struct( ...
    'trajectory_end_time', S(end), ...
    'late_window', options.late_window, ...
    'sample_interval', options.sample_interval, ...
    'late_median_distance_A', medians(1), ...
    'late_maximum_distance_A', maxima(1), ...
    'late_median_distance_B', medians(2), ...
    'late_maximum_distance_B', maxima(2), ...
    'late_median_distance_C', medians(3), ...
    'late_maximum_distance_C', maxima(3), ...
    'final_distance_A', distances(end,1), ...
    'final_distance_B', distances(end,2), ...
    'final_distance_C', distances(end,3), ...
    'outcome', outcome);
history = struct('sample_times', times, 'distance_A', distances(:,1), ...
    'distance_B', distances(:,2), 'distance_C', distances(:,3));
end

function label = outcome_label(valid)
names = {'A_stationary_neighborhood', 'B_periodic_reflection_neighborhood', ...
    'C_periodic_direct_neighborhood'};
if sum(valid) == 1
    label = names{find(valid, 1)};
elseif sum(valid) > 1
    label = 'ambiguous_multiple_endpoint_neighborhoods';
else
    label = 'unresolved_within_window';
end
end
