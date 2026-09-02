function output_file = plot_ml_BC_boundary_initial_labels(result_file)
%PLOT_ML_BC_BOUNDARY_INITIAL_LABELS Rebuild the labeling figure from data.

if nargin<1
    cfg=brusselator_ml_boundary_labeling_config();
    result_file=fullfile(cfg.output.root,cfg.experiment_name, ...
        'ml_BC_boundary_initial_labels.mat');
end
if ~isfile(result_file)
    error('plot_ml_BC_boundary_initial_labels:MissingResult', ...
        'Missing saved labeling result: %s',result_file);
end
loaded=load(result_file,'results'); results=loaded.results;
records=results.records; cfg=results.configuration;
outdir=fileparts(result_file);

fig=figure('Visible','off','Color','w','Position',[100,100,1250,420]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
labels=cfg.ml_labeling.reported_labels;
colors=[0 0.45 0.74; 0.85 0.33 0.10; 0.49 0.18 0.56; 0.35 0.35 0.35];
ax=nexttile; hold(ax,'on');
for k=1:numel(labels)
    selected=strcmp({records.final_label},labels{k});
    scatter(ax,[records(selected).alpha_1],[records(selected).alpha_2],70, ...
        colors(k,:),'filled','DisplayName',labels{k});
end
xlabel(ax,'\alpha_1'); ylabel(ax,'\alpha_2'); title(ax,'Adaptive PDE labels');
grid(ax,'on'); axis(ax,'equal'); legend(ax,'Location','best');
ax=nexttile;
scatter(ax,[records.alpha_1],[records.alpha_2],75,[records.tested_duration],'filled');
xlabel(ax,'\alpha_1'); ylabel(ax,'\alpha_2'); title(ax,'Tested frozen time');
grid(ax,'on'); axis(ax,'equal'); colorbar(ax);
ax=nexttile;
scatter(ax,[records.late_median_B],[records.late_median_C],70, ...
    [records.tested_duration],'filled');
xline(ax,cfg.edge_tracking.median_distance_threshold,':k', ...
    'HandleVisibility','off');
yline(ax,cfg.edge_tracking.median_distance_threshold,':k', ...
    'HandleVisibility','off');
set(ax,'XScale','log','YScale','log'); grid(ax,'on'); colorbar(ax);
xlabel(ax,'late median distance to B'); ylabel(ax,'late median distance to C');
title(ax,'Physical classifier distances');
sgtitle(fig,'Frozen b=10 initial-seed labeling; U denotes unresolved');
output_file=fullfile(outdir,'initial_seed_labeling.png');
exportgraphics(fig,output_file,'Resolution',300); close(fig);
end
