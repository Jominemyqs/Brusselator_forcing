function [diagnostic, phase_fixed] = brusselator_phase_section_diagnostics( ...
        template, S, U, V, Lx, phase_cfg)
%BRUSSELATOR_PHASE_SECTION_DIAGNOSTICS Audit a periodic-orbit phase section.
%   Tests transversality, crossing-count and crossing-phase stability under
%   temporal and spatial resampling, reflection equivalence of multiple
%   branches, and reproducibility across late source-trajectory cycles.

required = {'phase_fraction','period','U','V'};
for k = 1:numel(required)
    if ~isfield(template,required{k})
        error('brusselator_phase_section_diagnostics:InvalidTemplate', ...
            'Template is missing %s.', required{k});
    end
end
orbit_mean = mean(mean(template.V,2));
base = cyclic_crossings(template.phase_fraction, template.U, template.V, ...
    template.period, orbit_mean);
if base.count == 0
    error('brusselator_phase_section_diagnostics:NoCrossing', ...
        'The proposed orbit-centered mean-v section has no upward crossing.');
end
mean_v_signal = mean(template.V,2);
relative_transversality = abs(base.derivative)*template.period / ...
    max(max(mean_v_signal)-min(mean_v_signal),eps);
transversality_pass = all(relative_transversality >= ...
    phase_cfg.minimum_relative_transversality);

strides = [1,2,4];
temporal_counts = zeros(size(strides));
temporal_phase_errors = zeros(size(strides));
for k = 1:numel(strides)
    indices = 1:strides(k):size(template.U,1);
    reduced = cyclic_crossings(template.phase_fraction(indices), ...
        template.U(indices,:), template.V(indices,:), template.period, orbit_mean);
    temporal_counts(k) = reduced.count;
    temporal_phase_errors(k) = phase_set_distance( ...
        base.phase_fraction, reduced.phase_fraction);
end
crossing_count_stable = all(temporal_counts == base.count);
temporal_resampling_pass = crossing_count_stable && ...
    all(temporal_phase_errors <= phase_cfg.crossing_time_tolerance);

source_x = linspace(0,1,size(template.U,2));
coarse_x = linspace(0,1,200);
coarse_U = interp1(source_x,template.U',coarse_x,'pchip')';
coarse_V = interp1(source_x,template.V',coarse_x,'pchip')';
coarse = cyclic_crossings(template.phase_fraction, coarse_U, coarse_V, ...
    template.period, mean(mean(coarse_V,2)));
spatial_count_stable = coarse.count == base.count;
spatial_phase_error = phase_set_distance(base.phase_fraction, ...
    coarse.phase_fraction);
spatial_resampling_pass = spatial_count_stable && ...
    spatial_phase_error <= phase_cfg.crossing_time_tolerance;

[maximum_branch_distance, branches_reflection_equivalent] = ...
    compare_branches(base.U,base.V,Lx,phase_cfg.branch_equivalence_distance);
if branches_reflection_equivalent
    branch_policy = 'reflection_equivalent_branches_quotiented';
else
    branch_policy = 'retain_all_non_equivalent_branches';
end

source = source_crossings(S,U,V,orbit_mean, ...
    phase_cfg.phase_fixed_cycles_to_retain*template.period);
if source.count == 0
    error('brusselator_phase_section_diagnostics:NoSourceCrossing', ...
        'No late source-trajectory crossing is available.');
end
source_branch = zeros(source.count,1);
source_branch_distance = zeros(source.count,1);
for k = 1:source.count
    distances = zeros(base.count,1);
    for j = 1:base.count
        item = brusselator_physical_orbit_distance(source.U(k,:),source.V(k,:), ...
            base.U(j,:),base.V(j,:),Lx);
        distances(j) = item.distance;
    end
    [source_branch_distance(k),source_branch(k)] = min(distances);
end
branch_reproducibility_pass = max(source_branch_distance) < ...
    phase_cfg.branch_equivalence_distance;
expected_source_crossings = phase_cfg.phase_fixed_cycles_to_retain*base.count;
source_crossing_count_pass = abs(source.count-expected_source_crossings) <= 1;

orientation = brusselator_reflection_canonicalize(source.U,source.V,Lx);
phase_fixed = struct('times',source.times,'raw_U',source.U,'raw_V',source.V, ...
    'canonical_U',orientation.canonical_U, ...
    'canonical_V',orientation.canonical_V, ...
    'uses_reflection',orientation.uses_reflection, ...
    'branch_id',source_branch,'distance_to_template_branch',source_branch_distance);
diagnostic = struct( ...
    'period',template.period,'orbit_mean_v',orbit_mean, ...
    'number_crossings_per_period',base.count, ...
    'crossing_phase_fractions',base.phase_fraction, ...
    'crossing_derivatives',base.derivative, ...
    'relative_transversality',relative_transversality, ...
    'minimum_relative_transversality',min(relative_transversality), ...
    'transversality_pass',transversality_pass, ...
    'temporal_resampling_strides',strides, ...
    'temporal_resampling_counts',temporal_counts, ...
    'temporal_phase_errors',temporal_phase_errors, ...
    'crossing_count_stable',crossing_count_stable, ...
    'temporal_resampling_pass',temporal_resampling_pass, ...
    'spatial_resampling_count',coarse.count, ...
    'spatial_phase_error',spatial_phase_error, ...
    'spatial_resampling_pass',spatial_resampling_pass, ...
    'maximum_branch_distance_modulo_reflection',maximum_branch_distance, ...
    'branches_reflection_equivalent',branches_reflection_equivalent, ...
    'branch_policy',branch_policy, ...
    'number_late_source_crossings',source.count, ...
    'expected_late_source_crossings',expected_source_crossings, ...
    'source_crossing_count_pass',source_crossing_count_pass, ...
    'maximum_source_branch_distance',max(source_branch_distance), ...
    'branch_reproducibility_pass',branch_reproducibility_pass, ...
    'phase_section_pass',transversality_pass && temporal_resampling_pass && ...
        spatial_resampling_pass && source_crossing_count_pass && ...
        branch_reproducibility_pass);
end

function crossings = cyclic_crossings(phase,U,V,period,orbit_mean)
phase = phase(:);
q = mean(V,2)-orbit_mean;
n = numel(phase);
indices = zeros(0,1); alphas = zeros(0,1); derivatives = zeros(0,1);
phases = zeros(0,1); states_U = zeros(0,size(U,2)); states_V = states_U;
for j = 1:n
    next = mod(j,n)+1;
    if q(j) <= 0 && q(next) > 0
        if next == 1
            dt_fraction = 1-phase(j)+phase(next);
        else
            dt_fraction = phase(next)-phase(j);
        end
        alpha = -q(j)/(q(next)-q(j));
        crossing_phase = mod(phase(j)+alpha*dt_fraction,1);
        indices(end+1,1) = j; %#ok<AGROW>
        alphas(end+1,1) = alpha; %#ok<AGROW>
        derivatives(end+1,1) = (q(next)-q(j))/(dt_fraction*period); %#ok<AGROW>
        phases(end+1,1) = crossing_phase; %#ok<AGROW>
        states_U(end+1,:) = (1-alpha)*U(j,:)+alpha*U(next,:); %#ok<AGROW>
        states_V(end+1,:) = (1-alpha)*V(j,:)+alpha*V(next,:); %#ok<AGROW>
    end
end
[phases,order] = sort(phases);
crossings = struct('count',numel(phases),'phase_fraction',phases, ...
    'indices',indices(order),'alpha',alphas(order), ...
    'derivative',derivatives(order),'U',states_U(order,:), ...
    'V',states_V(order,:));
end

function source = source_crossings(S,U,V,orbit_mean,window_duration)
S = S(:);
late = S >= S(end)-window_duration;
S = S(late); U = U(late,:); V = V(late,:);
q = mean(V,2)-orbit_mean;
indices = find(q(1:end-1) <= 0 & q(2:end) > 0);
times = zeros(numel(indices),1);
states_U = zeros(numel(indices),size(U,2));
states_V = states_U;
for k = 1:numel(indices)
    j = indices(k);
    alpha = -q(j)/(q(j+1)-q(j));
    times(k) = S(j)+alpha*(S(j+1)-S(j));
    states_U(k,:) = (1-alpha)*U(j,:)+alpha*U(j+1,:);
    states_V(k,:) = (1-alpha)*V(j,:)+alpha*V(j+1,:);
end
source = struct('count',numel(times),'times',times,'U',states_U,'V',states_V);
end

function distance = phase_set_distance(reference,candidate)
reference = sort(reference(:)); candidate = sort(candidate(:));
if numel(reference) ~= numel(candidate) || isempty(reference)
    distance = Inf;
    return;
end
n = numel(reference);
distance = Inf;
for shift = 0:n-1
    shifted = circshift(candidate,shift);
    delta = abs(reference-shifted);
    delta = min(delta,1-delta);
    distance = min(distance,max(delta));
end
end

function [maximum, equivalent] = compare_branches(U,V,Lx,threshold)
if size(U,1) <= 1
    maximum = 0;
    equivalent = true;
    return;
end
maximum = 0;
for i = 1:size(U,1)
    for j = i+1:size(U,1)
        item = brusselator_physical_orbit_distance(U(i,:),V(i,:),U(j,:),V(j,:),Lx);
        maximum = max(maximum,item.distance);
    end
end
equivalent = maximum < threshold;
end
