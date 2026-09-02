function brusselator_write_metadata(outdir, cfg, metadata)
%BRUSSELATOR_WRITE_METADATA Save machine-readable run configuration/provenance.

if ~isfolder(outdir)
    mkdir(outdir);
end

save(fullfile(outdir, 'configuration_and_metadata.mat'), 'cfg', 'metadata');

fid = fopen(fullfile(outdir, 'configuration_and_metadata.json'), 'w');
if fid < 0
    error('brusselator_write_metadata:FileOpenFailed', ...
        'Could not create metadata JSON in %s.', outdir);
end
cleanup_file = onCleanup(@() fclose(fid)); %#ok<NASGU>
fwrite(fid, jsonencode(struct('configuration', cfg, 'metadata', metadata)), 'char');
end
