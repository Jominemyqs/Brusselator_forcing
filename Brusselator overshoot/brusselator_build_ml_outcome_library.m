function [library,provenance] = brusselator_build_ml_outcome_library(cfg)
%BRUSSELATOR_BUILD_ML_OUTCOME_LIBRARY Exact B/C plus validated other outcomes.
%   B and C come from independently Newton-converged N=400 periodic orbits.
%   A, R5, and R6 come from the validated Stage-2 library. No edge-orbit or
%   held-out challenge data are loaded.

required={'stage2_file','exact_B_file','exact_C_file'};
for k=1:numel(required)
    file=cfg.ml_labeling.(required{k});
    if ~isfile(file)
        error('brusselator_build_ml_outcome_library:MissingSource', ...
            'Missing %s: %s',required{k},file);
    end
end
stage2_loaded=load(cfg.ml_labeling.stage2_file,'results');
stage2=stage2_loaded.results;
B_loaded=load(cfg.ml_labeling.exact_B_file,'results'); B=B_loaded.results;
C_loaded=load(cfg.ml_labeling.exact_C_file,'results'); C=C_loaded.results;
validate_exact_orbit(B,'B',cfg); validate_exact_orbit(C,'C',cfg);

template=struct('label','','type','','period',NaN,'U',[],'V',[]);
library=repmat(template,5,1);
provenance_template=struct('label','','source_kind','','source_file','', ...
    'period',NaN,'phase_sample_count',NaN,'newton_residual',NaN, ...
    'independent_verification_residual',NaN);
provenance=repmat(provenance_template,5,1);

library(1)=stage2_reference(stage2,'A');
provenance(1)=stage2_provenance(library(1),cfg.ml_labeling.stage2_file);
[library(2),provenance(2)]=exact_reference(B,'B', ...
    cfg.ml_labeling.exact_B_file,cfg.ml_labeling.exact_periodic_phase_samples);
[library(3),provenance(3)]=exact_reference(C,'C', ...
    cfg.ml_labeling.exact_C_file,cfg.ml_labeling.exact_periodic_phase_samples);
library(4)=stage2_reference(stage2,'R5');
provenance(4)=stage2_provenance(library(4),cfg.ml_labeling.stage2_file);
library(5)=stage2_reference(stage2,'R6');
provenance(5)=stage2_provenance(library(5),cfg.ml_labeling.stage2_file);

if ~isequal({library.label},cfg.ml_labeling.raw_library_order)
    error('brusselator_build_ml_outcome_library:LibraryOrder', ...
        'The five-way library order differs from the frozen labeling contract.');
end
if any(arrayfun(@(x)size(x.U,2)~=cfg.grid.N || ...
        ~isequal(size(x.U),size(x.V)),library))
    error('brusselator_build_ml_outcome_library:GridMismatch', ...
        'At least one reference does not use the N=%d grid.',cfg.grid.N);
end
end

function validate_exact_orbit(orbit,label,cfg)
if ~strcmp(orbit.status,'newton_converged') || size(orbit.state,1)~=cfg.grid.N || ...
        abs(orbit.state(end,1)-cfg.grid.Lx)>1e-12 || ...
        orbit.summary.normalized_residual>1e-8 || ...
        orbit.verification.flow_residual>1e-8
    error('brusselator_build_ml_outcome_library:ExactOrbitGate', ...
        'Exact %s orbit failed the saved convergence or grid gate.',label);
end
end

function reference = stage2_reference(stage2,label)
index=find(strcmp({stage2.outcomes.label},label),1);
if isempty(index)
    error('brusselator_build_ml_outcome_library:MissingStage2Outcome', ...
        'Stage 2 has no validated %s outcome.',label);
end
source=stage2.outcomes(index);
if strcmp(source.type,'stationary')
    U=source.orbit_U; V=source.orbit_V;
else
    U=source.dense_template.U; V=source.dense_template.V;
end
reference=struct('label',label,'type',source.type,'period',source.period, ...
    'U',U,'V',V);
end

function record = stage2_provenance(reference,file)
record=struct('label',reference.label,'source_kind','validated_stage2_reference', ...
    'source_file',file,'period',reference.period, ...
    'phase_sample_count',size(reference.U,1),'newton_residual',NaN, ...
    'independent_verification_residual',NaN);
end

function [reference,record] = exact_reference(source,label,file,count)
phase=(0:count-1)'/count;
theta=source.orbit_times(:)/source.period;
[theta,indices]=unique(theta,'stable');
U=interp1(theta,source.orbit_U(indices,:),phase,'pchip');
V=interp1(theta,source.orbit_V(indices,:),phase,'pchip');
reference=struct('label',label,'type','periodic','period',source.period, ...
    'U',U,'V',V);
record=struct('label',label,'source_kind','exact_newton_periodic_orbit', ...
    'source_file',file,'period',source.period,'phase_sample_count',count, ...
    'newton_residual',source.summary.normalized_residual, ...
    'independent_verification_residual',source.verification.flow_residual);
end
