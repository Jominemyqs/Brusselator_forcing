function results = run_periodic_edge_orbit_unstable_branches()
%RUN_PERIODIC_EDGE_ORBIT_UNSTABLE_BRANCHES Follow both Floquet branches.
%   Perturbs the Newton-converged orbit along the real unstable Floquet
%   eigenvector at two amplitudes and classifies every trajectory against
%   A/B/C/R5/R6 with the established physical classifier.

cfg = brusselator_periodic_edge_branch_config();
outdir = fullfile(cfg.output.root,cfg.experiment_name);
checkpoint = fullfile(outdir,'branch_checkpoint.mat');
final_file = fullfile(outdir,'periodic_edge_orbit_unstable_branches.mat');
if isfile(final_file)
    error('run_periodic_edge_orbit_unstable_branches:OutputExists', ...
        'Refusing to overwrite completed result: %s',final_file);
end
if ~isfile(cfg.branches.floquet_file) || ...
        ~isfile(cfg.dynamic_edge.parent_library_file)
    error('run_periodic_edge_orbit_unstable_branches:MissingSource', ...
        'The Floquet result or five-way orbit library is missing.');
end
loaded = load(cfg.branches.floquet_file,'results'); floquet = loaded.results;
loaded = load(cfg.dynamic_edge.parent_library_file,'results');
library = loaded.results.library;
if numel(floquet.unstable_indices) ~= 1
    error('run_periodic_edge_orbit_unstable_branches:UnstableDimension', ...
        'Branch test requires exactly one computed nontrivial unstable multiplier.');
end
index = floquet.unstable_indices(1);
if abs(imag(floquet.multipliers(index))) > 1e-8
    error('run_periodic_edge_orbit_unstable_branches:ComplexUnstableMode', ...
        'The leading unstable multiplier is not real.');
end
state = floquet.state; x = state(:,1); N = numel(x);
w = spatial_weights(x); sqrt_weights = sqrt([w;w]);
y = [state(:,2);state(:,3)]; state_norm = norm(sqrt_weights.*y);
direction_q = real(floquet.weighted_eigenvectors(:,index));
direction_q = direction_q/norm(direction_q);
direction_y = direction_q./sqrt_weights;

case_specs = build_case_specs(cfg);
record_template = empty_record();
if isfile(checkpoint)
    saved = load(checkpoint,'records','raw_files','next_case');
    records = saved.records; raw_files = saved.raw_files; next_case = saved.next_case;
    fprintf('Resuming unstable-branch test at case %d.\n',next_case);
else
    if isfolder(outdir)
        error('run_periodic_edge_orbit_unstable_branches:IncompleteOutput', ...
            'Existing output lacks a checkpoint: %s',outdir);
    end
    mkdir(outdir); mkdir(fullfile(outdir,'raw_cases'));
    records = record_template([]); raw_files = {}; next_case = 1;
    brusselator_write_metadata(outdir,cfg, ...
        brusselator_run_metadata(cfg,mfilename,NaN));
    save(checkpoint,'cfg','records','raw_files','next_case','-v7.3');
end
timer = tic;

for k = next_case:numel(case_specs)
    spec = case_specs(k);
    delta_y = spec.sign*spec.relative_amplitude*state_norm*direction_y;
    initial_state = state;
    initial_state(:,2) = initial_state(:,2)+delta_y(1:N);
    initial_state(:,3) = initial_state(:,3)+delta_y(N+1:end);
    raw_file = fullfile(outdir,'raw_cases',[spec.case_id,'.mat']);
    partial_file = fullfile(outdir,'raw_cases',[spec.case_id,'_partial.mat']);
    if isfile(raw_file)
        loaded = load(raw_file,'edge_case'); edge_case = loaded.edge_case;
    else
        edge_case = brusselator_evolve_and_classify_state( ...
            spec.case_id,initial_state,library,cfg,partial_file);
        save(raw_file,'edge_case','-v7.3');
    end
    record = record_template;
    record.case_id = spec.case_id;
    record.relative_amplitude = spec.relative_amplitude;
    record.sign = spec.sign;
    record.initial_physical_relative_distance = ...
        state_distance(initial_state,state,cfg.grid.Lx);
    record.outcome = edge_case.classification.outcome;
    record.tested_duration = edge_case.duration;
    record.decision_time = edge_case.decision_time;
    record.runtime_seconds = edge_case.runtime_seconds;
    record.late_median_distance_A = edge_case.classification.late_median_distances(1);
    record.late_median_distance_B = edge_case.classification.late_median_distances(2);
    record.late_median_distance_C = edge_case.classification.late_median_distances(3);
    record.late_median_distance_R5 = edge_case.classification.late_median_distances(4);
    record.late_median_distance_R6 = edge_case.classification.late_median_distances(5);
    record.raw_file = raw_file;
    records(end+1,1) = record; %#ok<AGROW>
    raw_files{end+1,1} = raw_file; %#ok<AGROW>
    next_case = k+1;
    save(checkpoint,'cfg','records','raw_files','next_case','-v7.3');
    writetable(struct2table(records,'AsArray',true), ...
        fullfile(outdir,'branch_progress.csv'));
end

summary = summarize(records,floquet,cfg);
writetable(struct2table(records,'AsArray',true), ...
    fullfile(outdir,'unstable_branch_outcomes.csv'));
writetable(struct2table(summary,'AsArray',true), ...
    fullfile(outdir,'unstable_branch_summary.csv'));
make_figure(outdir,records,library);
results = struct('configuration',cfg,'floquet_file',cfg.branches.floquet_file, ...
    'multiplier',floquet.multipliers(index), ...
    'weighted_unstable_eigenvector',direction_q, ...
    'physical_unstable_eigenvector',direction_y, ...
    'records',records,'raw_files',{raw_files},'summary',summary);
save(final_file,'results','-v7.3');
brusselator_write_metadata(outdir,cfg, ...
    brusselator_run_metadata(cfg,mfilename,toc(timer)));
disp(struct2table(records,'AsArray',true));
disp(struct2table(summary,'AsArray',true));
fprintf('Unstable Floquet branches saved in: %s\n',outdir);
end

function specs = build_case_specs(cfg)
template = struct('case_id','','relative_amplitude',NaN,'sign',NaN);
specs = repmat(template,numel(cfg.branches.relative_amplitudes)*2,1); k = 0;
for amplitude = cfg.branches.relative_amplitudes
    for sign_value = cfg.branches.signs
        k = k+1;
        specs(k) = struct('case_id',sprintf('epsilon_%s_sign_%s', ...
            strrep(sprintf('%.0e',amplitude),'-','m'),sign_token(sign_value)), ...
            'relative_amplitude',amplitude,'sign',sign_value);
    end
end
end

function token = sign_token(value)
if value < 0, token = 'minus'; else, token = 'plus'; end
end

function record = empty_record()
record = struct('case_id','','relative_amplitude',NaN,'sign',NaN, ...
    'initial_physical_relative_distance',NaN,'outcome','', ...
    'tested_duration',NaN,'decision_time',NaN,'runtime_seconds',NaN, ...
    'late_median_distance_A',NaN,'late_median_distance_B',NaN, ...
    'late_median_distance_C',NaN,'late_median_distance_R5',NaN, ...
    'late_median_distance_R6',NaN,'raw_file','');
end

function summary = summarize(records,floquet,cfg)
amplitudes = unique([records.relative_amplitude]);
opposite = true;
mapping_consistent = true;
reference_minus = ''; reference_plus = '';
for k = 1:numel(amplitudes)
    subset = records([records.relative_amplitude]==amplitudes(k));
    minus = subset([subset.sign]<0).outcome;
    plus = subset([subset.sign]>0).outcome;
    opposite = opposite && isequal(sort({minus,plus}),sort(cfg.branches.required_outcomes));
    if k == 1
        reference_minus = minus; reference_plus = plus;
    else
        mapping_consistent = mapping_consistent && ...
            strcmp(minus,reference_minus) && strcmp(plus,reference_plus);
    end
end
if opposite && mapping_consistent
    status = 'opposite_BC_unstable_branches_supported';
else
    status = 'unstable_branch_outcomes_not_cleanly_BC';
end
summary = struct('status',status,'multiplier_real', ...
    real(floquet.multipliers(floquet.unstable_indices(1))), ...
    'amplitude_count',numel(amplitudes), ...
    'opposite_BC_at_every_amplitude',opposite, ...
    'sign_mapping_consistent',mapping_consistent, ...
    'minus_outcome',reference_minus,'plus_outcome',reference_plus, ...
    'maximum_decision_time',max([records.decision_time]), ...
    'interpretation',cfg.branches.interpretation);
end

function distance = state_distance(A,B,Lx)
N = size(A,1); dx = Lx/(N-1); w = dx*ones(N,1); w([1,end])=0.5*dx;
numerator = sum(((A(:,2)-B(:,2)).^2+(A(:,3)-B(:,3)).^2).*w);
energy_A = sum((A(:,2).^2+A(:,3).^2).*w);
energy_B = sum((B(:,2).^2+B(:,3).^2).*w);
distance = sqrt(numerator/max(0.5*(energy_A+energy_B),eps));
end

function weights = spatial_weights(x)
dx = x(2)-x(1); weights = dx*ones(numel(x),1); weights([1,end])=0.5*dx;
end

function make_figure(outdir,records,library)
fig = figure('Color','w','Position',[100 100 1050 760]);
layout = tiledlayout(fig,numel(records),1);
colors = lines(numel(library));
for j = 1:numel(records)
    loaded = load(records(j).raw_file,'edge_case'); data = loaded.edge_case;
    ax = nexttile(layout); hold(ax,'on');
    for k = 1:numel(library)
        semilogy(ax,data.history.sample_times,data.history.distances(:,k), ...
            'Color',colors(k,:),'LineWidth',1,'DisplayName',library(k).label);
    end
    title(ax,sprintf('%s: %s',records(j).case_id,records(j).outcome), ...
        'Interpreter','none'); ylabel(ax,'orbit distance'); grid(ax,'on'); box(ax,'on');
    if j==1, legend(ax,'Location','eastoutside'); end
end
xlabel(nexttile(layout,numel(records)),'frozen time');
exportgraphics(fig,fullfile(outdir,'unstable_branch_distances.png'),'Resolution',300);
close(fig);
end
