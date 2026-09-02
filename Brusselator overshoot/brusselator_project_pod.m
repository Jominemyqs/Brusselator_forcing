function [scores, errors] = brusselator_project_pod(U, V, model, component_counts)
%BRUSSELATOR_PROJECT_POD Project full states and report scaled-norm errors.

if ~isequal(size(U),size(V)) || size(U,2) ~= numel(model.mean_U)
    error('brusselator_project_pod:DimensionMismatch', ...
        'States do not match the POD model spatial grid.');
end
centered_U = U-model.mean_U;
centered_V = V-model.mean_V;
transformed = [(centered_U.*model.sqrt_spatial_weights)/max(model.scale_U,eps), ...
    (centered_V.*model.sqrt_spatial_weights)/max(model.scale_V,eps)];
scores = transformed*model.modes;
if nargin < 4 || isempty(component_counts)
    component_counts = model.maximum_components;
end
component_counts = component_counts(component_counts <= model.maximum_components);
errors = zeros(size(U,1),numel(component_counts));
denominator = sqrt(sum(transformed.^2,2));
for k = 1:numel(component_counts)
    count = component_counts(k);
    reconstruction = scores(:,1:count)*model.modes(:,1:count)';
    errors(:,k) = sqrt(sum((transformed-reconstruction).^2,2)) ./ ...
        max(denominator,eps);
end
end
