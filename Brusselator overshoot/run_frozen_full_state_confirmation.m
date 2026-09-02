function results = run_frozen_full_state_confirmation()
%RUN_FROZEN_FULL_STATE_CONFIRMATION Resolve the late frozen b=10 dynamics.
%   Restarts the saved control and bmax=11.14 final states from the long
%   trajectory study, then evolves both under the autonomous baseline system
%   with fine temporal output. The run stores both u and v, so the comparison
%   is no longer limited to the v-field or to unit-spaced output times.
%
%   This is a finite-time confirmation of a candidate recurrent outcome. It
%   does not by itself prove the existence or stability of an attractor.

study_id = 'frozen_full_state_confirmation_v1';
parent_dir = fullfile('experiment_outputs', 'extended_postforcing_reference_v2');
parent_file = fullfile(parent_dir, 'raw_extended_trajectories.mat');
outdir = fullfile('experiment_outputs', study_id);
if ~isfile(parent_file)
    error('run_frozen_full_state_confirmation:MissingParentData', ...
        'Expected long-trajectory data at %s.', parent_file);
end
if isfolder(outdir)
    error('run_frozen_full_state_confirmation:OutputExists', ...
        'Refusing to overwrite existing experiment directory: %s', outdir);
end

parent = load(parent_file, 'raw');
raw = parent.raw;
case_index = find(cellfun(@(s) abs(s.bmax - 11.14) < 1e-12, ...
    raw.forced_trajectories), 1);
if isempty(case_index)
    error('run_frozen_full_state_confirmation:MissingDepartureCase', ...
        'The parent output has no bmax = 11.14 trajectory.');
end

confirmation_duration = 200;
solver_settings = struct('name', 'ode15s', 'output_dt', 0.05, ...
    'RelTol', 1e-8, 'AbsTol', 1e-10, 'MaxStep', 0.05);
cfg = brusselator_default_config(struct( ...
    'experiment_name', study_id, ...
    'model', raw.config.model, ...
    'grid', raw.config.grid, ...
    'forcing', raw.config.forcing, ...
    'initial', raw.config.initial, ...
    'solver', solver_settings));
cfg.confirmation = struct( ...
    'parent_raw_file', parent_file, ...
    'parent_final_time', raw.config.time.Tfinal, ...
    'frozen_b', raw.config.forcing.b0, ...
    'source_bmax', raw.forced_trajectories{case_index}.bmax, ...
    'duration', confirmation_duration, ...
    'full_state_metric', 'relative_l2_over_concatenated_u_v');

mkdir(outdir);
brusselator_write_metadata(outdir, cfg, ...
    brusselator_run_metadata(cfg, mfilename, NaN));

frozen_par = brusselator_make_parameters(cfg, @(t) cfg.forcing.b0); %#ok<NASGU>
fprintf('Restarting control at b = %.2f for %.0f time units ...\n', ...
    cfg.confirmation.frozen_b, confirmation_duration);
control_timer = tic;
[control_final, S, V_control, U_control] = solve_brusselator_1d_forced( ...
    raw.control.ic_final, frozen_par, confirmation_duration, 0);
control_runtime = toc(control_timer);

fprintf('Restarting bmax = %.2f departure state at b = %.2f for %.0f time units ...\n', ...
    cfg.confirmation.source_bmax, cfg.confirmation.frozen_b, confirmation_duration);
forced_timer = tic;
[forced_final, S_forced, V_forced, U_forced] = solve_brusselator_1d_forced( ...
    raw.forced_trajectories{case_index}.ic_final, frozen_par, ...
    confirmation_duration, 0);
forced_runtime = toc(forced_timer);
if ~isequal(S, S_forced)
    error('run_frozen_full_state_confirmation:TimeMismatch', ...
        'Control and forced restart output times must match.');
end

alignment = full_state_alignment(S, U_forced, V_forced, U_control, V_control, 5);
control_variation = full_state_temporal_variation(U_control, V_control);
forced_variation = full_state_temporal_variation(U_forced, V_forced);
control_section = full_state_poincare(S, U_control, V_control);
forced_section = full_state_poincare(S, U_forced, V_forced);
profile = compare_mean_v_profiles(V_forced, V_control);

record = struct();
record.parent_final_time = cfg.confirmation.parent_final_time;
record.confirmation_duration = confirmation_duration;
record.output_dt = cfg.solver.output_dt;
record.MaxStep = cfg.solver.MaxStep;
record.source_bmax = cfg.confirmation.source_bmax;
record.same_time_rel_diff_full_state = alignment.same_time_distance;
record.phase_reflection_rel_diff_full_state = alignment.best_distance;
record.best_lag = alignment.best_lag;
record.best_reflection = alignment.best_reflection;
record.control_relative_temporal_variation = control_variation;
record.forced_relative_temporal_variation = forced_variation;
record.control_section_period = control_section.period;
record.forced_section_period = forced_section.period;
record.forced_section_period_cv = forced_section.period_cv;
record.forced_section_recurrence_rel_diff_full_state = forced_section.recurrence;
record.forced_section_direct_recurrence_rel_diff_full_state = ...
    forced_section.direct_recurrence;
record.forced_section_reflection_recurrence_rel_diff_full_state = ...
    forced_section.reflection_recurrence;
record.forced_section_recurrence_uses_reflection = ...
    forced_section.recurrence_uses_reflection;
record.mean_profile_rel_diff_v = profile.distance;
record.profile_uses_reflection = profile.uses_reflection;
record.interpretation = classify_confirmation(record, cfg.classification);

summary_table = struct2table(record, 'AsArray', true);
writetable(summary_table, fullfile(outdir, 'frozen_full_state_summary.csv'));
write_section_table(outdir, control_section, forced_section);

make_confirmation_figures(outdir, S, U_control, V_control, U_forced, V_forced, ...
    alignment, control_section, forced_section, cfg);

results = struct();
results.configuration = cfg;
results.summary = record;
results.alignment = alignment;
results.control_section = control_section;
results.forced_section = forced_section;
results.profile = profile;
results.trajectories = struct('S', S, 'U_control', U_control, ...
    'V_control', V_control, 'U_forced', U_forced, 'V_forced', V_forced, ...
    'control_final', control_final, 'forced_final', forced_final);
save(fullfile(outdir, 'frozen_full_state_trajectories.mat'), 'results', '-v7.3');

metadata = brusselator_run_metadata(cfg, mfilename, control_runtime + forced_runtime);
brusselator_write_metadata(outdir, cfg, metadata);
disp(summary_table);
fprintf('Frozen full-state confirmation complete. Results saved in: %s\n', outdir);
end

function alignment = full_state_alignment(S, Uf, Vf, Uc, Vc, max_lag_time)
dt = median(diff(S));
max_lag_steps = round(max_lag_time / dt);
coarse_stride = 4;
lags = (-max_lag_steps:max_lag_steps)';
no_reflection = zeros(size(lags));
reflection = zeros(size(lags));
for j = 1:numel(lags)
    no_reflection(j) = full_state_shifted_distance(Uf, Vf, Uc, Vc, lags(j), ...
        false, coarse_stride);
    reflection(j) = full_state_shifted_distance(Uf, Vf, Uc, Vc, lags(j), ...
        true, coarse_stride);
end
[minimum_identity, identity_idx] = min(no_reflection);
[minimum_reflection, reflection_idx] = min(reflection);
if minimum_reflection < minimum_identity
    reflection_flag = true;
    coarse_best = lags(reflection_idx);
else
    reflection_flag = false;
    coarse_best = lags(identity_idx);
end

fine_lags = max(-max_lag_steps, coarse_best - coarse_stride): ...
    min(max_lag_steps, coarse_best + coarse_stride);
fine_distances = zeros(size(fine_lags));
for j = 1:numel(fine_lags)
    fine_distances(j) = full_state_shifted_distance(Uf, Vf, Uc, Vc, ...
        fine_lags(j), reflection_flag, 1);
end
[best_distance, best_idx] = min(fine_distances);

alignment = struct();
alignment.lags = lags * dt;
alignment.no_reflection = no_reflection;
alignment.reflection = reflection;
alignment.best_distance = best_distance;
alignment.best_lag = fine_lags(best_idx) * dt;
alignment.best_reflection = reflection_flag;
alignment.same_time_distance = full_state_shifted_distance(Uf, Vf, Uc, Vc, 0, false, 1);
end

function distance = full_state_shifted_distance(Uf, Vf, Uc, Vc, lag, reflect_control, stride)
count = size(Uf, 1);
if lag >= 0
    idx_f = 1:(count - lag);
    idx_c = (1 + lag):count;
else
    idx_f = (1 - lag):count;
    idx_c = 1:(count + lag);
end
idx_f = idx_f(1:stride:end);
idx_c = idx_c(1:stride:end);
Uref = Uc(idx_c,:);
Vref = Vc(idx_c,:);
if reflect_control
    Uref = fliplr(Uref);
    Vref = fliplr(Vref);
end
numerator = sum((Uf(idx_f,:) - Uref).^2, 'all') + ...
    sum((Vf(idx_f,:) - Vref).^2, 'all');
denominator = sum(Uref.^2, 'all') + sum(Vref.^2, 'all');
distance = sqrt(numerator / max(denominator, eps));
end

function variation = full_state_temporal_variation(U, V)
Uc = U - mean(U, 1);
Vc = V - mean(V, 1);
numerator = mean(sum(Uc.^2, 2) + sum(Vc.^2, 2));
denominator = mean(sum(U.^2, 2) + sum(V.^2, 2));
variation = sqrt(numerator / max(denominator, eps));
end

function section = full_state_poincare(S, U, V)
q = mean(V, 2);
q = q - mean(q);
cross_idx = find(q(1:end-1) <= 0 & q(2:end) > 0);
times = zeros(numel(cross_idx), 1);
states = zeros(numel(cross_idx), 2 * size(U,2));
for j = 1:numel(cross_idx)
    i = cross_idx(j);
    alpha = -q(i) / (q(i+1) - q(i));
    times(j) = S(i) + alpha * (S(i+1) - S(i));
    u = (1 - alpha) * U(i,:) + alpha * U(i+1,:);
    v = (1 - alpha) * V(i,:) + alpha * V(i+1,:);
    states(j,:) = [u, v];
end
if numel(times) >= 2
    periods = diff(times);
    period = median(periods);
    period_cv = std(periods) / max(mean(periods), eps);
    previous = states(1:end-1,:);
    next = states(2:end,:);
    norms = sqrt(sum(previous.^2, 2));
    direct_by_crossing = sqrt(sum((next - previous).^2, 2)) ./ max(norms, eps);
    N = size(U,2);
    reflected_next = [fliplr(next(:,1:N)), fliplr(next(:,N+1:end))];
    reflection_by_crossing = sqrt(sum((reflected_next - previous).^2, 2)) ./ ...
        max(norms, eps);
    direct_recurrence = median(direct_by_crossing);
    reflection_recurrence = median(reflection_by_crossing);
    if reflection_recurrence < direct_recurrence
        recurrence = reflection_recurrence;
        recurrence_uses_reflection = true;
        recurrence_by_crossing = reflection_by_crossing;
    else
        recurrence = direct_recurrence;
        recurrence_uses_reflection = false;
        recurrence_by_crossing = direct_by_crossing;
    end
else
    period = NaN;
    period_cv = NaN;
    recurrence = NaN;
    direct_recurrence = NaN;
    reflection_recurrence = NaN;
    recurrence_uses_reflection = false;
    recurrence_by_crossing = zeros(0,1);
end
section = struct('q', q, 'times', times, 'states', states, ...
    'period', period, 'period_cv', period_cv, 'recurrence', recurrence, ...
    'direct_recurrence', direct_recurrence, ...
    'reflection_recurrence', reflection_recurrence, ...
    'recurrence_uses_reflection', recurrence_uses_reflection, ...
    'recurrence_by_crossing', recurrence_by_crossing);
end

function profile = compare_mean_v_profiles(Vf, Vc)
mean_forced = mean(Vf, 1);
mean_control = mean(Vc, 1);
direct = norm(mean_forced - mean_control) / max(norm(mean_control), eps);
reflected = norm(mean_forced - fliplr(mean_control)) / max(norm(mean_control), eps);
if reflected < direct
    profile = struct('distance', reflected, 'uses_reflection', true, ...
        'mean_forced', mean_forced, 'mean_control', mean_control);
else
    profile = struct('distance', direct, 'uses_reflection', false, ...
        'mean_forced', mean_forced, 'mean_control', mean_control);
end
end

function label = classify_confirmation(record, classification)
if record.phase_reflection_rel_diff_full_state > classification.transition_threshold && ...
        record.forced_relative_temporal_variation > 1e-3 && ...
        record.forced_section_recurrence_rel_diff_full_state < classification.recovery_threshold && ...
        record.forced_section_period_cv < 0.05
    label = 'persistent_distinct_recurrent_full_state_candidate';
else
    label = 'unresolved_finite_time_full_state_dynamics';
end
end

function write_section_table(outdir, control, forced)
records = struct('trajectory', {}, 'crossing_time', {});
for j = 1:numel(control.times)
    records(end+1) = struct('trajectory', 'control', ... %#ok<AGROW>
        'crossing_time', control.times(j));
end
for j = 1:numel(forced.times)
    records(end+1) = struct('trajectory', 'forced', ... %#ok<AGROW>
        'crossing_time', forced.times(j));
end
if ~isempty(records)
    writetable(struct2table(records, 'AsArray', true), ...
        fullfile(outdir, 'full_state_poincare_crossings.csv'));
end
end

function make_confirmation_figures(outdir, S, Uc, Vc, Uf, Vf, alignment, ...
        control_section, forced_section, cfg)
fig = figure('Color', 'w', 'Position', [100 100 1000 650]);
subplot(2,2,1);
plot(alignment.lags, alignment.no_reflection, 'LineWidth', 1.2, ...
    'DisplayName', 'identity'); hold on;
plot(alignment.lags, alignment.reflection, 'LineWidth', 1.2, ...
    'DisplayName', 'reflection');
xline(alignment.best_lag, '--', 'best lag', 'HandleVisibility', 'off');
xlabel('Control time lag'); ylabel('Relative full-state distance');
title('Phase/reflection-aware comparison'); legend('Location', 'best'); grid on;

subplot(2,2,2);
distance = sqrt(sum((Uf-Uc).^2 + (Vf-Vc).^2, 2)) ./ ...
    max(sqrt(sum(Uc.^2 + Vc.^2, 2)), eps);
plot(S, distance, 'LineWidth', 1.2); grid on;
xlabel('Frozen restart time'); ylabel('Relative full-state distance');
title('Same-time separation under frozen b=10');

subplot(2,2,3);
plot(S, control_section.q, 'LineWidth', 1.0, 'DisplayName', 'control'); hold on;
plot(S, forced_section.q, 'LineWidth', 1.0, 'DisplayName', 'forced');
xlabel('Frozen restart time'); ylabel('Centered spatial-mean v');
title('Physical Poincare observable'); legend('Location', 'best'); grid on;

subplot(2,2,4);
if numel(forced_section.times) >= 2
    semilogy(2:numel(forced_section.times), forced_section.recurrence_by_crossing, ...
        '-o', 'LineWidth', 1.0);
    xlabel('Poincare crossing number'); ylabel('Relative full-state return distance');
    if forced_section.recurrence_uses_reflection
        title('Successive section-state recurrence (with reflection)');
    else
        title('Successive section-state recurrence');
    end
    grid on;
else
    axis off; text(0.5, 0.5, 'No forced section crossings', ...
        'HorizontalAlignment', 'center');
end
sgtitle(sprintf('Frozen b=%.2f full-state confirmation from t=%.0f', ...
    cfg.confirmation.frozen_b, cfg.confirmation.parent_final_time));
exportgraphics(fig, fullfile(outdir, 'full_state_phase_and_recurrence.png'), ...
    'Resolution', 300);
close(fig);

fig = figure('Color', 'w', 'Position', [100 100 900 500]);
imagesc(linspace(0, cfg.grid.Lx, size(Vf,2)), S, Vf);
axis xy; colorbar; xlabel('x'); ylabel('Frozen restart time');
title(sprintf('Fine-output v(x,t) after b_{max}=%.2f forcing', ...
    cfg.confirmation.source_bmax));
exportgraphics(fig, fullfile(outdir, 'forced_fine_spacetime_v.png'), 'Resolution', 300);
close(fig);
end
