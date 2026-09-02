function result = brusselator_physical_orbit_distance(u, v, ref_U, ref_V, Lx)
%BRUSSELATOR_PHYSICAL_ORBIT_DISTANCE Physical full-state orbit distance.
%   Uses trapezoidal spatial quadrature in the original u/v units and
%   minimizes a relative full-state distance over saved temporal phases and
%   exact reflection. No POD component scaling is used.

u = u(:)';
v = v(:)';
if numel(u) ~= numel(v) || size(ref_U,2) ~= numel(u) || ...
        ~isequal(size(ref_U), size(ref_V))
    error('brusselator_physical_orbit_distance:DimensionMismatch', ...
        'Candidate and reference fields must use the same spatial grid.');
end
N = numel(u);
dx = Lx/(N-1);
weights = dx * ones(1,N);
weights([1,end]) = 0.5*dx;
direct_num = sum(((ref_U-u).^2 + (ref_V-v).^2).*weights, 2);
reflected_U = fliplr(ref_U);
reflected_V = fliplr(ref_V);
reflected_num = sum(((reflected_U-u).^2 + ...
    (reflected_V-v).^2).*weights, 2);
denominator = sum((ref_U.^2 + ref_V.^2).*weights, 2);
direct = sqrt(direct_num ./ max(denominator,eps));
reflected_denominator = sum((reflected_U.^2 + reflected_V.^2).*weights, 2);
reflected = sqrt(reflected_num ./ max(reflected_denominator,eps));
[direct_minimum, direct_index] = min(direct);
[reflection_minimum, reflection_index] = min(reflected);
if reflection_minimum < direct_minimum
    result = struct('distance', reflection_minimum, ...
        'uses_reflection', true, 'reference_index', reflection_index, ...
        'direct_distance', direct_minimum, ...
        'reflection_distance', reflection_minimum, ...
        'norm', 'quadrature_weighted_relative_full_state_L2_original_units');
else
    result = struct('distance', direct_minimum, ...
        'uses_reflection', false, 'reference_index', direct_index, ...
        'direct_distance', direct_minimum, ...
        'reflection_distance', reflection_minimum, ...
        'norm', 'quadrature_weighted_relative_full_state_L2_original_units');
end
end
