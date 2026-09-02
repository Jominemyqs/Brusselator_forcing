function result = brusselator_reflection_canonicalize(U, V, Lx)
%BRUSSELATOR_REFLECTION_CANONICALIZE Choose an exact-reflection orientation.
%   Uses the first sufficiently nonzero coefficient in the ordered sequence
%   (u,k=1), (v,k=1), (u,k=3), (v,k=3), (u,k=5), (v,k=5). All raw,
%   reflected, and selected states are returned so quotient decisions remain
%   auditable.

if isvector(U), U = U(:)'; end
if isvector(V), V = V(:)'; end
if ~isequal(size(U), size(V))
    error('brusselator_reflection_canonicalize:DimensionMismatch', ...
        'U and V must have identical dimensions.');
end
N = size(U,2);
x = linspace(0, Lx, N);
modes = [1, 3, 5];
features = zeros(size(U,1), 2*numel(modes));
for k = 1:numel(modes)
    basis = cos(modes(k) * pi * x / Lx);
    features(:,2*k-1) = trapz(x, U .* basis, 2);
    features(:,2*k) = trapz(x, V .* basis, 2);
end
energy = sqrt(trapz(x, U.^2 + V.^2, 2));
tolerance = 1e-10 * max(energy * sqrt(Lx), 1);
selected_feature = zeros(size(U,1),1);
selected_index = zeros(size(U,1),1);
uses_reflection = false(size(U,1),1);
for row = 1:size(U,1)
    index = find(abs(features(row,:)) > tolerance(row), 1);
    if isempty(index)
        selected_feature(row) = 0;
        selected_index(row) = 0;
    else
        selected_feature(row) = features(row,index);
        selected_index(row) = index;
        uses_reflection(row) = selected_feature(row) < 0;
    end
end
canonical_U = U;
canonical_V = V;
canonical_U(uses_reflection,:) = fliplr(U(uses_reflection,:));
canonical_V(uses_reflection,:) = fliplr(V(uses_reflection,:));
result = struct('raw_U', U, 'raw_V', V, ...
    'reflected_U', fliplr(U), 'reflected_V', fliplr(V), ...
    'canonical_U', canonical_U, 'canonical_V', canonical_V, ...
    'uses_reflection', uses_reflection, 'features', features, ...
    'selected_feature_index', selected_index, ...
    'selected_feature_value', selected_feature, ...
    'feature_order', {{'u_cos1','v_cos1','u_cos3','v_cos3','u_cos5','v_cos5'}}, ...
    'tolerance', tolerance);
end
