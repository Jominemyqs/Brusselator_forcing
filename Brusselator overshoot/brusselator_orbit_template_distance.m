function result = brusselator_orbit_template_distance(A, B, varargin)
%BRUSSELATOR_ORBIT_TEMPLATE_DISTANCE Compare full periodic orbit templates.
%   Minimizes a symmetric relative RMS full-state distance over sampled
%   temporal phase and exact spatial reflection. Different spatial grids are
%   interpolated in normalized x/L to a common comparison grid.

parser = inputParser;
addParameter(parser, 'common_spatial_points', 801, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 3);
parse(parser, varargin{:});
common_N = parser.Results.common_spatial_points;
if size(A.U,1) ~= size(B.U,1)
    error('brusselator_orbit_template_distance:PhaseSampleMismatch', ...
        'Templates must use the same number of normalized phase samples.');
end
validate_template(A);
validate_template(B);

coordinate = linspace(0, 1, common_N);
A_U = interp1(linspace(0,1,size(A.U,2)), A.U', coordinate, 'pchip')';
A_V = interp1(linspace(0,1,size(A.V,2)), A.V', coordinate, 'pchip')';
B_U = interp1(linspace(0,1,size(B.U,2)), B.U', coordinate, 'pchip')';
B_V = interp1(linspace(0,1,size(B.V,2)), B.V', coordinate, 'pchip')';
X = [A_U, A_V];
B_direct = [B_U, B_V];
B_reflected = [fliplr(B_U), fliplr(B_V)];
normalizer = sqrt(0.5 * (mean(sum(X.^2,2)) + ...
    mean(sum(B_direct.^2,2))));
n = size(X,1);
direct = zeros(n,1);
reflected = zeros(n,1);
for shift = 0:n-1
    direct(shift+1) = sqrt(mean(sum((X - ...
        circshift(B_direct, shift, 1)).^2, 2))) / max(normalizer, eps);
    reflected(shift+1) = sqrt(mean(sum((X - ...
        circshift(B_reflected, shift, 1)).^2, 2))) / max(normalizer, eps);
end
[direct_minimum, direct_index] = min(direct);
[reflection_minimum, reflection_index] = min(reflected);
if reflection_minimum < direct_minimum
    distance = reflection_minimum;
    uses_reflection = true;
    phase_fraction = (reflection_index - 1) / n;
else
    distance = direct_minimum;
    uses_reflection = false;
    phase_fraction = (direct_index - 1) / n;
end
result = struct('distance', distance, ...
    'uses_reflection', uses_reflection, ...
    'phase_fraction', phase_fraction, ...
    'direct_distance', direct_minimum, ...
    'direct_phase_fraction', (direct_index - 1) / n, ...
    'reflection_distance', reflection_minimum, ...
    'reflection_phase_fraction', (reflection_index - 1) / n, ...
    'common_spatial_points', common_N, ...
    'phase_samples', n);
end

function validate_template(template)
if ~isstruct(template) || ~isfield(template, 'U') || ~isfield(template, 'V') || ...
        size(template.U,1) ~= size(template.V,1) || ...
        size(template.U,2) ~= size(template.V,2)
    error('brusselator_orbit_template_distance:InvalidTemplate', ...
        'Each template must contain consistently sized U and V arrays.');
end
end
