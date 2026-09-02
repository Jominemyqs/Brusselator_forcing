function repair_periodic_edge_resolution_provenance()
%REPAIR_PERIODIC_EDGE_RESOLUTION_PROVENANCE Correct inherited N=400 labels.
%   The continuation calculations used the correct state sizes and spatial
%   grids. This idempotent repair changes only descriptive configuration and
%   summary strings in the newly generated N=200 and N=600 case files.

root = fullfile('experiment_outputs', ...
    'periodic_edge_orbit_resolution_continuation_v1');
resolutions = [200,600];
for N = resolutions
    newton_dir = fullfile(root,sprintf('N%04d_newton',N));
    newton_file = fullfile(newton_dir,'periodic_edge_orbit_newton.mat');
    loaded = load(newton_file,'results'); results = loaded.results;
    results.configuration.newton.interpretation = ...
        'periodic orbit of the configured semidiscrete frozen system';
    results.summary.claim_scope = sprintf(['a converged result verifies a ', ...
        'periodic orbit of the N=%d semidiscrete frozen system; it does not ', ...
        'yet give Floquet stability or continuum convergence'],N);
    save(newton_file,'results','-v7.3');
    repair_metadata(newton_dir,results.configuration);

    floquet_dir = fullfile(root,sprintf('N%04d_floquet',N));
    floquet_file = fullfile(floquet_dir,'periodic_edge_orbit_floquet.mat');
    loaded = load(floquet_file,'results'); results = loaded.results;
    results.configuration.floquet.interpretation = ...
        'leading Floquet multipliers of the configured semidiscrete orbit';
    results.summary.interpretation = sprintf(['leading Floquet multipliers ', ...
        'of the N=%d semidiscrete periodic orbit; one trivial unit ', ...
        'multiplier is required by time phase'],N);
    save(floquet_file,'results','-v7.3');
    repair_metadata(floquet_dir,results.configuration);
end
fprintf('Corrected descriptive resolution provenance for N=200 and N=600.\n');
end

function repair_metadata(folder,cfg)
file = fullfile(folder,'configuration_and_metadata.mat');
loaded = load(file,'metadata'); metadata = loaded.metadata;
brusselator_write_metadata(folder,cfg,metadata);
end
