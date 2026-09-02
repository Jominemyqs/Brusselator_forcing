function plot_ml_BC_edge_recovery()
%PLOT_ML_BC_EDGE_RECOVERY Rebuild the recovery diagnostic from saved data.

outdir=fullfile('experiment_outputs','ml_BC_edge_recovery_v1');
loaded=load(fullfile(outdir,'ml_BC_edge_recovery.mat'),'results');
models=loaded.results.models;cfg=loaded.results.configuration;
fig=figure('Visible','off','Color','w','Position',[100 100 1000 560]);
ax=axes(fig);hold(ax,'on');
for k=1:numel(models)
    if isempty(models(k).distance_times),continue;end
    semilogy(ax,models(k).distance_times,max(models(k).edge_orbit_distances,eps), ...
        'LineWidth',1.5,'DisplayName',models(k).model);
end
yline(ax,cfg.ml_edge.edge_neighborhood_threshold,'--k', ...
    'edge neighborhood','HandleVisibility','off','LabelHorizontalAlignment','right');
xlabel(ax,'frozen time');ylabel(ax,'distance to exact E_{BC} orbit');
title_handle=title(ax,'Dynamical usefulness of label-blind ML brackets');
set(ax,'Color','w','XColor','k','YColor','k','GridColor',[0.75 0.75 0.75]);
title_handle.Color='k';grid(ax,'on');box(ax,'on');
legend_handle=legend(ax,'Location','northeast');
set(legend_handle,'Color','w','TextColor','k','EdgeColor',[0.2 0.2 0.2]);
text(ax,0.02,0.93,['active: no verified B/C pair among ', ...
    'four frozen proposals'],'Units','normalized','Color',[0.55 0.1 0.1], ...
    'FontWeight','bold');
exportgraphics(fig,fullfile(outdir,'ML_bracket_edge_recovery.png'),'Resolution',300);
close(fig);
end
