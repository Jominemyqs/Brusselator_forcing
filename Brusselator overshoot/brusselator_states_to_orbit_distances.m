function [distances, phase_fractions, reflections] = ...
        brusselator_states_to_orbit_distances(U, V, template)
%BRUSSELATOR_STATES_TO_ORBIT_DISTANCES Distance states to a periodic orbit.
%   Minimizes relative full-state distance over sampled temporal phase and
%   exact spatial reflection. Candidate states and the template must share a
%   spatial grid.

if size(U,1) ~= size(V,1) || size(U,2) ~= size(V,2) || ...
        size(U,2) ~= size(template.U,2)
    error('brusselator_states_to_orbit_distances:DimensionMismatch', ...
        'Candidate states and the template must share a spatial grid.');
end
template_U_reflected = fliplr(template.U);
template_V_reflected = fliplr(template.V);
template_norm = max(sqrt(sum(template.U.^2 + template.V.^2, 2)), eps);
reflected_norm = max(sqrt(sum(template_U_reflected.^2 + ...
    template_V_reflected.^2, 2)), eps);
n = size(U,1);
distances = zeros(n,1);
phase_fractions = zeros(n,1);
reflections = false(n,1);
for k = 1:n
    direct = sqrt(sum((template.U - U(k,:)).^2 + ...
        (template.V - V(k,:)).^2, 2)) ./ template_norm;
    reflected = sqrt(sum((template_U_reflected - U(k,:)).^2 + ...
        (template_V_reflected - V(k,:)).^2, 2)) ./ reflected_norm;
    [direct_minimum, direct_index] = min(direct);
    [reflection_minimum, reflection_index] = min(reflected);
    if reflection_minimum < direct_minimum
        distances(k) = reflection_minimum;
        phase_fractions(k) = template.phase_fraction(reflection_index);
        reflections(k) = true;
    else
        distances(k) = direct_minimum;
        phase_fractions(k) = template.phase_fraction(direct_index);
    end
end
end
