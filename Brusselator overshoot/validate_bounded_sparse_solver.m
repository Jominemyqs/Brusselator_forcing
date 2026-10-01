function validate_bounded_sparse_solver()
root=fullfile('experiment_outputs','bounded_BC_transfer_20261001_v1');fd=fullfile(root,'flow_cosine_b10');
G=load(fullfile(fd,'geometry.mat'));cfg=G.cfg;cfg.solver.analytic_jacobian=true;
ids=[6,8];records=struct([]);
for id=ids
    q=load(fullfile(fd,'labels',sprintf('state_%03d.mat',id)));old=q.case_data;
    par=brusselator_make_parameters(cfg,@(t)10);timer=tic;
    [st,S,V,U]=solve_brusselator_1d_forced(old.initial_state,par,old.duration,0);elapsed=toc(timer);
    [c,~]=brusselator_classify_five_way_trajectory(S,U,V,G.library);
    err=norm(st(:,2:3)-old.final_state(:,2:3),'fro')/norm(old.final_state(:,2:3),'fro');
    r=struct('state_index',id,'duration',old.duration,'dense_label',old.classification.outcome, ...
        'sparse_label',c.outcome,'relative_final_state_difference',err,'dense_runtime_seconds',old.runtime_seconds, ...
        'sparse_integration_seconds',elapsed,'passed',strcmp(c.outcome,old.classification.outcome)&&err<1e-3);
    if isempty(records),records=r;else,records(end+1)=r;end;disp(r); %#ok<AGROW>
end
f=fopen(fullfile(root,'sparse_solver_validation.json'),'w');fprintf(f,'%s\n',jsonencode(records,PrettyPrint=true));fclose(f);
assert(all([records.passed]),'Sparse solver validation failed.');
end
