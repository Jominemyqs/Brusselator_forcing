function metadata = brusselator_run_metadata(cfg, script_name, runtime_seconds)
%BRUSSELATOR_RUN_METADATA Collect provenance for a configured experiment.

if nargin < 3
    runtime_seconds = NaN;
end

metadata = struct();
metadata.schema_version = cfg.schema_version;
metadata.script = script_name;
metadata.created_utc = char(datetime('now', 'TimeZone', 'UTC', ...
    'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z'''));
metadata.matlab_version = version;
metadata.matlab_release = version('-release');
metadata.platform = computer;
metadata.working_directory = pwd;
metadata.runtime_seconds = runtime_seconds;

[status, commit] = system('git rev-parse HEAD');
metadata.git_available = (status == 0);
if metadata.git_available
    metadata.git_commit = strtrim(commit);
else
    metadata.git_commit = '';
end
end
