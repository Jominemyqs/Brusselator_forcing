function [refs, cfg] = brusselator_three_way_references()
%BRUSSELATOR_THREE_WAY_REFERENCES Load phase-aware A, B, and C references.
%   A is the stationary frozen baseline. B and C are full-state periodic
%   templates. Reflection is the exact spatial symmetry used by downstream
%   distance calculations; B and C templates also quotient temporal phase.

slice_file = fullfile('experiment_outputs', 'phase_aware_basin_slice_v1', ...
    'raw_phase_aware_basin_slice.mat');
candidate_file = fullfile('experiment_outputs', ...
    'edge_candidate_frozen_confirmation_v1', ...
    'candidate_confirmation_trajectory.mat');
if ~isfile(slice_file) || ~isfile(candidate_file)
    error('brusselator_three_way_references:MissingInput', ...
        'The saved phase-aware slice and candidate-C confirmation are required.');
end

loaded = load(slice_file, 'results');
slice = loaded.results;
loaded = load(candidate_file, 'results');
candidate = loaded.results;
cfg = slice.configuration;

refs = struct('A', slice.endpoint_A, 'B', slice.template);
S = candidate.trajectory.S;
U = candidate.trajectory.U;
V = candidate.trajectory.V;
period = candidate.section.period;
crossing_times = candidate.section.times(:);
valid = find(crossing_times >= 0.55 * S(end) & crossing_times + period <= S(end));
if isempty(valid)
    valid = find(crossing_times + period <= S(end));
end
if isempty(valid)
    error('brusselator_three_way_references:TemplateWindow', ...
        'No complete late candidate-C reference period is available.');
end

start_index = valid(1);
start_time = crossing_times(start_index);
phase_times = (0:cfg.basin_slice.template_phase_step:period)';
if phase_times(end) < period
    phase_times(end+1,1) = period;
end
state = candidate.section.states(start_index,:);
N = cfg.grid.N;
refs.C = struct('phase_times', phase_times, ...
    'U', interp1(S, U, start_time + phase_times, 'linear'), ...
    'V', interp1(S, V, start_time + phase_times, 'linear'));
refs.C.U(1,:) = state(1:N);
refs.C.V(1,:) = state(N+1:end);
refs.info = struct('grid_N', N, 'domain_Lx', cfg.grid.Lx, ...
    'A_type', 'stationary', 'B_period', refs.B.phase_times(end), ...
    'C_period', period, 'C_section_time', start_time, ...
    'slice_file', slice_file, 'candidate_file', candidate_file);
end
