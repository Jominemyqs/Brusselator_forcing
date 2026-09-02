function distances = brusselator_three_way_distances(u, v, refs)
%BRUSSELATOR_THREE_WAY_DISTANCES Distances to A, B, and C modulo symmetry.
%   A is minimized over identity/reflection. B and C are minimized over all
%   saved temporal phases and identity/reflection. The full (u,v) state and
%   the reference state provide the relative-L2 normalization.

u = u(:)';
v = v(:)';
if numel(u) ~= numel(v) || numel(u) ~= numel(refs.A.U)
    error('brusselator_three_way_distances:DimensionMismatch', ...
        'Candidate and reference fields must use the same grid.');
end
distances = [distance_stationary(u, v, refs.A), ...
    distance_orbit(u, v, refs.B), distance_orbit(u, v, refs.C)];
end

function distance = distance_stationary(u, v, reference)
direct = sqrt(sum((u - reference.U).^2) + sum((v - reference.V).^2));
reflected = sqrt(sum((u - fliplr(reference.U)).^2) + ...
    sum((v - fliplr(reference.V)).^2));
normalizer = max(sqrt(sum(reference.U.^2) + sum(reference.V.^2)), eps);
distance = min(direct, reflected) / normalizer;
end

function distance = distance_orbit(u, v, template)
direct = sqrt(sum((template.U - u).^2 + (template.V - v).^2, 2)) ./ ...
    max(sqrt(sum(template.U.^2 + template.V.^2, 2)), eps);
reflected = sqrt(sum((fliplr(template.U) - u).^2 + ...
    (fliplr(template.V) - v).^2, 2)) ./ ...
    max(sqrt(sum(fliplr(template.U).^2 + fliplr(template.V).^2, 2)), eps);
distance = min([direct; reflected]);
end
