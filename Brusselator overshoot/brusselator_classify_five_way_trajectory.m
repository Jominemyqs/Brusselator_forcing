function [record,history] = brusselator_classify_five_way_trajectory( ...
        S,U,V,library,varargin)
%BRUSSELATOR_CLASSIFY_FIVE_WAY_TRAJECTORY Classify against A/B/C/R5/R6.
%   Every distance is a quadrature-weighted relative full-state L2 distance
%   in the original u/v variables. Periodic references are minimized over
%   saved temporal phase and exact reflection. The function retains
%   unresolved and ambiguous outcomes rather than choosing the nearest class.

parser=inputParser;
addParameter(parser,'late_window',20,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(parser,'sample_interval',1,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(parser,'median_threshold',1e-2,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(parser,'maximum_threshold',2e-2,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(parser,'Lx',40,@(x)isnumeric(x)&&isscalar(x)&&x>0);
parse(parser,varargin{:});
options=parser.Results;

if size(U,1)~=numel(S) || ~isequal(size(U),size(V))
    error('brusselator_classify_five_way_trajectory:DimensionMismatch', ...
        'S, U, and V dimensions are inconsistent.');
end
if S(end)-S(1)<options.late_window
    error('brusselator_classify_five_way_trajectory:ShortTrajectory', ...
        'The trajectory is shorter than the requested late window.');
end
required={'label','type','U','V'};
for j=1:numel(library)
    if ~all(isfield(library(j),required)) || ...
            size(library(j).U,2)~=size(U,2) || ...
            ~isequal(size(library(j).U),size(library(j).V))
        error('brusselator_classify_five_way_trajectory:InvalidLibrary', ...
            'Reference %d is incomplete or uses another spatial grid.',j);
    end
end

dt=median(diff(S));
stride=max(1,round(options.sample_interval/dt));
indices=unique([1:stride:numel(S),numel(S)]);
times=S(indices);
distances=zeros(numel(indices),numel(library));
reflections=false(size(distances));
phase_indices=zeros(size(distances));
for k=1:numel(indices)
    for j=1:numel(library)
        item=brusselator_physical_orbit_distance( ...
            U(indices(k),:),V(indices(k),:),library(j).U,library(j).V,options.Lx);
        distances(k,j)=item.distance;
        reflections(k,j)=item.uses_reflection;
        phase_indices(k,j)=item.reference_index;
    end
end
late=times>=S(end)-options.late_window;
medians=median(distances(late,:),1);
maxima=max(distances(late,:),[],1);
valid=medians<options.median_threshold & maxima<options.maximum_threshold;
if sum(valid)==1
    outcome=library(find(valid,1)).label;
elseif sum(valid)>1
    outcome='ambiguous_multiple_outcomes';
else
    outcome='unresolved';
end
record=struct('outcome',outcome,'trajectory_end_time',S(end), ...
    'late_window',options.late_window,'sample_interval',options.sample_interval, ...
    'labels',{{library.label}},'late_median_distances',medians, ...
    'late_maximum_distances',maxima,'final_distances',distances(end,:), ...
    'valid_outcomes',valid,'median_threshold',options.median_threshold, ...
    'maximum_threshold',options.maximum_threshold, ...
    'norm','quadrature_weighted_relative_full_state_L2_original_units');
history=struct('sample_indices',indices,'sample_times',times, ...
    'labels',{{library.label}},'distances',distances, ...
    'uses_reflection',reflections,'reference_phase_index',phase_indices);
end
