function plot_extended_postforcing_results(outdir)
%PLOT_EXTENDED_POSTFORCING_RESULTS Recreate trajectory figures from raw data.
%   This intentionally reads only saved output from
%   RUN_EXTENDED_POSTFORCING_TRAJECTORIES, so figures can be regenerated
%   without rerunning the costly simulations.

raw_file = fullfile(outdir, 'raw_extended_trajectories.mat');
if ~isfile(raw_file)
    error('plot_extended_postforcing_results:MissingRawData', ...
        'Expected saved trajectory data at %s.', raw_file);
end

loaded = load(raw_file, 'raw');
raw = loaded.raw;
cfg = raw.config;
S_control = raw.control.S;
x = linspace(0, cfg.grid.Lx, cfg.grid.N);
num_cases = numel(raw.forced_trajectories);

fig = figure('Color', 'w', 'Position', [100 100 850 500]);
hold on;
handles = gobjects(num_cases, 1);
for k = 1:num_cases
    bmax = raw.forced_trajectories{k}.bmax;
    handles(k) = semilogy(S_control, max(raw.distance_to_control{k}, 1e-16), ...
        'LineWidth', 1.5, 'DisplayName', sprintf('b_{max}=%.2f', bmax));
end
xline(cfg.time.forcing_end, '--', 'forcing ends', 'HandleVisibility', 'off');
yline(cfg.classification.recovery_threshold, '--', 'recovery threshold', ...
    'HandleVisibility', 'off');
ylabel('Relative distance to matched control in v');
xlabel('Time');
title('Extended post-forcing separation from the control trajectory');
legend(handles, 'Location', 'best');
grid on; box on;
exportgraphics(fig, fullfile(outdir, 'distance_to_control_over_time.png'), ...
    'Resolution', 300);
close(fig);

late_start = cfg.time.Tfinal - 1000;
idx = S_control >= late_start;
for k = 1:num_cases
    bmax = raw.forced_trajectories{k}.bmax;
    fig = figure('Color', 'w', 'Position', [100 100 900 500]);
    imagesc(x, S_control(idx), raw.forced_trajectories{k}.V(idx,:));
    axis xy; colorbar;
    xlabel('x'); ylabel('t');
    title(sprintf('Late-time v(x,t), b_{max}=%.2f', bmax));
    exportgraphics(fig, fullfile(outdir, ...
        sprintf('late_time_spacetime_v_bmax_%0.2f.png', bmax)), ...
        'Resolution', 300);
    close(fig);
end
end
