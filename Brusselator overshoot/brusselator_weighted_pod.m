function model = brusselator_weighted_pod(U, V, sample_weights, Lx, ...
        separate_component_scaling, maximum_components)
%BRUSSELATOR_WEIGHTED_POD Fit an explicitly weighted full-state POD basis.
%   Samples are centered with their normalized weights. Each centered row is
%   multiplied by sqrt(sample_weight) before SVD. Spatial trapezoidal weights
%   are always applied; optional separate u/v RMS scales affect POD geometry
%   only and are saved in the returned model.

if ~isequal(size(U),size(V)) || size(U,1) ~= numel(sample_weights)
    error('brusselator_weighted_pod:DimensionMismatch', ...
        'U, V, and sample weights have inconsistent sizes.');
end
weights = sample_weights(:);
if any(~isfinite(weights)) || any(weights <= 0)
    error('brusselator_weighted_pod:InvalidWeights', ...
        'Every POD sample weight must be positive and finite.');
end
weights = weights/sum(weights);
N = size(U,2);
dx = Lx/(N-1);
spatial_weights = dx*ones(1,N);
spatial_weights([1,end]) = 0.5*dx;
sqrt_spatial_weights = sqrt(spatial_weights);
mean_U = weights'*U;
mean_V = weights'*V;
centered_U = U-mean_U;
centered_V = V-mean_V;
if separate_component_scaling
    scale_U = sqrt(sum(weights .* sum(centered_U.^2.*spatial_weights,2)));
    scale_V = sqrt(sum(weights .* sum(centered_V.^2.*spatial_weights,2)));
else
    scale_U = 1;
    scale_V = 1;
end
transformed = [(centered_U.*sqrt_spatial_weights)/max(scale_U,eps), ...
    (centered_V.*sqrt_spatial_weights)/max(scale_V,eps)];
weighted = transformed.*sqrt(weights);
[~, singular_matrix, modes] = svd(weighted,'econ');
singular_values = diag(singular_matrix);
component_count = min([maximum_components,size(modes,2),numel(singular_values)]);
modes = modes(:,1:component_count);
singular_values = singular_values(1:component_count);
all_energy = diag(singular_matrix).^2;
explained = 100*singular_values.^2/max(sum(all_energy),eps);
scores = transformed*modes;
model = struct('mean_U',mean_U,'mean_V',mean_V, ...
    'sample_weights',weights,'spatial_weights',spatial_weights, ...
    'sqrt_spatial_weights',sqrt_spatial_weights, ...
    'separate_component_scaling',separate_component_scaling, ...
    'scale_U',scale_U,'scale_V',scale_V,'modes',modes, ...
    'singular_values',singular_values,'explained_percent',explained, ...
    'training_scores',scores,'maximum_components',component_count, ...
    'fit_rule', ['weighted mean; centered rows multiplied by sqrt(sample weight) ', ...
        'before SVD']);
end
