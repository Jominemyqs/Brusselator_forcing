function run_bounded_BC_transfer(action,family,part,total)
% Independent bounded transfer audit; never writes manuscript inputs.
if nargin<3,part=1;end
if nargin<4,total=1;end
root=fullfile('experiment_outputs','bounded_BC_transfer_20261001_v1');
protocol=jsondecode(fileread(fullfile(root,'protocol.json')));
fd=fullfile(root,family); if ~isfolder(fd),mkdir(fd);end
cfg=brusselator_ml_boundary_labeling_config();
b=10;if strcmp(family,'parameter_b10p02'),b=10.02;end
cfg.edge_tracking.frozen_b=b;cfg.model.b0=b;
if strcmp(action,'prepare')
    if isfile(fullfile(fd,'geometry.mat')) || isfile(fullfile(fd,'gate_failed.json'))
        fprintf('Preparation already recorded for %s.\n',family);return;
    end
    try
        prepare(fd,family,cfg,protocol);
    catch e
        writejson(fullfile(fd,'gate_failed.json'),struct('family',family, ...
            'status','preparation_gate_failed','identifier',e.identifier, ...
            'message',e.message,'report',getReport(e,'extended')));
        fprintf('%s preparation gate failed: %s\n',family,e.message);
    end
    return;
end
G=load(fullfile(fd,'geometry.mat'));cfg=G.cfg;library=G.library;
cfg.solver.analytic_jacobian=true; % Execution amendment, validated against dense runs.
if strcmp(action,'initial') || strcmp(action,'endpoints')
    selection=readtable(fullfile(fd,[action,'_selection.csv']));
    for k=part:total:height(selection)
        id=selection.state_index(k);state=state_at(G,id);
        label_cached(fd,sprintf('state_%03d',id),state,library,cfg);
    end
elseif strcmp(action,'audit')
    rep=part;P=readtable(fullfile(fd,sprintf('proposals_%d.csv',rep)),'TextType','string');
    D=readtable(fullfile(fd,'designs.csv'));ids=D.state_index(D.replicate==rep);
    labels=cell(size(ids));times=zeros(size(ids));
    for k=1:numel(ids)
        d=label_cached(fd,sprintf('state_%03d',ids(k)),state_at(G,ids(k)),library,cfg);
        labels{k}=d.classification.outcome;times(k)=d.runtime_seconds;
    end
    methods={'pair_aware','pointwise_uncertainty','space_filling','direct_refinement'};
    records=struct([]);
    for mi=1:4
        name=methods{mi};out=fullfile(fd,sprintf('audit_%d_%s',rep,name));
        final=fullfile(out,'summary.json');
        if isfile(final),continue;end
        if ~isfolder(out),mkdir(out);end
        rec=struct('family',family,'replicate',rep,'method',name,'status','no_short_BC_pair', ...
            'initial_labels',numel(ids),'initial_runtime_seconds',sum(times), ...
            'acquisition_labels',0,'acquisition_runtime_seconds',0,'pair_yield',0, ...
            'labels_to_first_short_pair',NaN,'short_pair_found',false, ...
            'initial_normalized_separation',NaN,'initial_physical_separation',NaN, ...
            'guard_labels',0,'guard_runtime_seconds',0,'guard_complete',false, ...
            'minimum_edge_distance',NaN,'shadow_duration',0,'period',NaN, ...
            'period_relative_error',NaN,'edge_recovered',false,'diagnostic_runtime_seconds',0);
        left=[];right=[];
        if mi<4
            rows=P(P.method==string(name),:);rows=sortrows(rows,'proposal_rank');
            for j=1:height(rows)
                a=rows.first_state_index(j);c=rows.second_state_index(j);
                da=label_cached(fd,sprintf('state_%03d',a),state_at(G,a),library,cfg);
                dc=label_cached(fd,sprintf('state_%03d',c),state_at(G,c),library,cfg);
                rec.acquisition_labels=rec.acquisition_labels+2;
                rec.acquisition_runtime_seconds=rec.acquisition_runtime_seconds+da.runtime_seconds+dc.runtime_seconds;
                labs={da.classification.outcome,dc.classification.outcome};
                if all(ismember({'B','C'},labs))
                    rec.pair_yield=rec.pair_yield+1;
                    if isempty(left)
                        if strcmp(labs{1},'B'),bi=a;ci=c;else,bi=c;ci=a;end
                        left=state_at(G,bi);right=state_at(G,ci);
                        rec.labels_to_first_short_pair=2*j;
                        rec.initial_normalized_separation=norm(G.Z(bi,:)-G.Z(ci,:));
                    end
                end
            end
        else
            Bids=ids(strcmp(labels,'B'));Cids=ids(strcmp(labels,'C'));
            best=inf;bi=NaN;ci=NaN;
            for a=Bids(:)'
                for c=Cids(:)'
                    d=norm(G.Z(a,:)-G.Z(c,:));
                    if d<best,best=d;bi=a;ci=c;end
                end
            end
            if isfinite(best)
                left=state_at(G,bi);right=state_at(G,ci);zb=G.Z(bi,:);zc=G.Z(ci,:);
                while norm(zb-zc)>protocol.short_pair_max && rec.acquisition_labels<protocol.endpoint_budget
                    step=rec.acquisition_labels+1;mid=interpolate(left,right);zm=(zb+zc)/2;
                    data=label_cached(out,sprintf('direct_%02d',step),mid,library,cfg);
                    rec.acquisition_labels=step;rec.acquisition_runtime_seconds=rec.acquisition_runtime_seconds+data.runtime_seconds;
                    if strcmp(data.classification.outcome,'B'),left=mid;zb=zm;
                    elseif strcmp(data.classification.outcome,'C'),right=mid;zc=zm;
                    else,left=[];right=[];rec.status=['direct_obstructed_',data.classification.outcome];break;end
                end
                if ~isempty(left) && norm(zb-zc)<=protocol.short_pair_max
                    rec.pair_yield=1;rec.labels_to_first_short_pair=rec.acquisition_labels;
                    rec.initial_normalized_separation=norm(zb-zc);
                else,left=[];right=[];end
            end
        end
        if ~isempty(left)
            rec.short_pair_found=true;rec.status='short_pair_found';
            rec.initial_physical_separation=physical_separation(left,right);
            for step=1:protocol.guard_steps
                mid=interpolate(left,right);data=label_cached(out,sprintf('guard_%02d',step),mid,library,cfg);
                rec.guard_labels=step;rec.guard_runtime_seconds=rec.guard_runtime_seconds+data.runtime_seconds;
                if strcmp(data.classification.outcome,'B'),left=mid;
                elseif strcmp(data.classification.outcome,'C'),right=mid;
                else,rec.status=['guard_obstructed_',data.classification.outcome];break;end
            end
            rec.guard_complete=rec.guard_labels==protocol.guard_steps && strcmp(rec.status,'short_pair_found');
            if rec.guard_complete
                candidate=interpolate(left,right);candidate_file=fullfile(out,'candidate.mat');
                if isfile(candidate_file),q=load(candidate_file);else
                    timer=tic;par=brusselator_make_parameters(cfg,@(t)b);
                    [~,S,V,U]=solve_brusselator_1d_forced(candidate,par,protocol.candidate_duration,0);
                    ts=(0:0.5:S(end))';Us=interp1(S,U,ts,'pchip');Vs=interp1(S,V,ts,'pchip');ds=zeros(size(ts));
                    for k=1:numel(ts)
                        d=brusselator_physical_orbit_distance(Us(k,:),Vs(k,:),G.edge.orbit_U,G.edge.orbit_V,cfg.grid.Lx);ds(k)=d.distance;
                    end
                    change=diff([false;ds<protocol.edge_radius;false]);starts=find(change==1);ends=find(change==-1)-1;
                    shadow=0;period=NaN;
                    if ~isempty(starts)
                        [shadow,j]=max(ts(ends)-ts(starts));
                        if shadow>=protocol.minimum_shadow
                            mask=S>=ts(starts(j)) & S<=ts(ends(j));
                            recurrence=brusselator_direct_period_recurrence(S(mask),U(mask,:),V(mask,:),ts(starts(j)),G.edge.period);
                            period=recurrence.best_period;
                        end
                    end
                    runtime=toc(timer);save(candidate_file,'candidate','S','U','V','ts','ds','shadow','period','runtime','-v7.3');q=load(candidate_file);
                end
                rec.minimum_edge_distance=min(q.ds);rec.shadow_duration=q.shadow;rec.period=q.period;
                rec.period_relative_error=abs(q.period-G.edge.period)/G.edge.period;
                rec.diagnostic_runtime_seconds=q.runtime;
                rec.edge_recovered=rec.minimum_edge_distance<protocol.edge_radius && ...
                    rec.shadow_duration>=protocol.minimum_shadow && rec.period_relative_error<protocol.period_tolerance;
                rec.status='audit_complete';
            end
        end
        rec.total_label_integrations=rec.initial_labels+rec.acquisition_labels+rec.guard_labels;
        rec.total_integration_runtime_seconds=rec.initial_runtime_seconds+rec.acquisition_runtime_seconds+rec.guard_runtime_seconds+rec.diagnostic_runtime_seconds;
        writejson(final,rec);disp(rec);
    end
else,error('Unknown action %s',action);end
end

function prepare(fd,family,cfg,p)
timer=tic;
if strcmp(family,'flow_cosine_b10')
    [library,~]=brusselator_build_ml_outcome_library(cfg);
    q=load(fullfile('experiment_outputs','periodic_edge_orbit_newton_v1','periodic_edge_orbit_newton.mat'));edge=q.results;
    q=load(fullfile('experiment_outputs','physical_ramp_BA_dynamic_edge_v1','raw_cases','event01_step08.mat'));B=q.edge_case.initial_state;
    q=load(fullfile('experiment_outputs','physical_ramp_BA_dynamic_edge_v1','raw_cases','event01_step10.mat'));C=q.edge_case.initial_state;
    par=brusselator_make_parameters(cfg,@(t)10);
    B=solve_brusselator_1d_forced(B,par,p.flow_advance,0);C=solve_brusselator_1d_forced(C,par,p.flow_advance,0);
    pilot_labels=2;pilot_runtime=0;
    for lab={'B','C'}
        if strcmp(lab{1},'B'),st=B;else,st=C;end
        d=label_cached(fd,['anchor_',lab{1}],st,library,cfg);pilot_runtime=pilot_runtime+d.runtime_seconds;
        assert(strcmp(d.classification.outcome,lab{1}),'Evolved anchor did not retain its label.');
    end
else
    names={'B','C','edge'};orbits=cell(3,1);floquets=cell(3,1);
    for k=1:3
        if k<3,file=fullfile('experiment_outputs','BC_edge_resolution_geometry_v1',[names{k},'_N0400_newton'],'periodic_edge_orbit_newton.mat');
        else,file=fullfile('experiment_outputs','periodic_edge_orbit_newton_v1','periodic_edge_orbit_newton.mat');end
        q=load(file);source=q.results;ncfg=brusselator_periodic_edge_newton_config();ncfg.edge_tracking.frozen_b=cfg.edge_tracking.frozen_b;
        ncfg.model.b0=cfg.edge_tracking.frozen_b;ncfg.output.root=fd;ncfg.experiment_name=[names{k},'_newton'];
        seed=struct('state',source.state,'period',source.period,'reference_state',source.state, ...
            'reference_tangent',brusselator_frozen_rhs(source.state,ncfg,cfg.edge_tracking.frozen_b),'event_index',NaN);
        seedfile=fullfile(fd,[names{k},'_seed.mat']);if ~isfile(seedfile),save(seedfile,'seed');end
        ncfg.newton.seed_file=seedfile;of=fullfile(fd,ncfg.experiment_name,'periodic_edge_orbit_newton.mat');
        if isfile(of),q=load(of);orbit=q.results;else,orbit=run_periodic_edge_orbit_newton(ncfg);end
        assert(strcmp(orbit.status,'newton_converged') && orbit.verification.flow_residual<1e-8,'Target-parameter orbit validation failed.');
        fcfg=brusselator_periodic_edge_floquet_config();fcfg.edge_tracking.frozen_b=cfg.edge_tracking.frozen_b;
        fcfg.model.b0=cfg.edge_tracking.frozen_b;fcfg.output.root=fd;fcfg.experiment_name=[names{k},'_floquet'];fcfg.floquet.orbit_file=of;
        ff=fullfile(fd,fcfg.experiment_name,'periodic_edge_orbit_floquet.mat');
        if isfile(ff),q=load(ff);fl=q.results;else,fl=run_periodic_edge_orbit_floquet(fcfg);end
        assert(fl.summary.arnoldi_flag==0 && max(fl.eigen_residuals)<1e-7 && abs(fl.summary.phase_multiplier_real-1)<1e-6,'Floquet accuracy gate failed.');
        expected=double(k==3);assert(fl.summary.nontrivial_unstable_count==expected,'Attractor/edge stability gate failed.');
        orbits{k}=orbit;floquets{k}=fl;
    end
    library=[reference(orbits{1},'B'),reference(orbits{2},'C')];edge=orbits{3};
    % At this changed parameter only validated B/C are named; all nonmatches stay unresolved.
    B=orbits{1}.state;C=orbits{2}.state;pilot_labels=0;pilot_runtime=0;
    for lab={'B','C'}
        if strcmp(lab{1},'B'),st=B;else,st=C;end
        d=label_cached(fd,['anchor_',lab{1}],st,library,cfg);pilot_labels=pilot_labels+1;pilot_runtime=pilot_runtime+d.runtime_seconds;
        assert(strcmp(d.classification.outcome,lab{1}),'Target-parameter anchor validation failed.');
    end
    for step=1:p.pilot_steps
        mid=interpolate(B,C);d=label_cached(fd,sprintf('pilot_%02d',step),mid,library,cfg);
        pilot_labels=pilot_labels+1;pilot_runtime=pilot_runtime+d.runtime_seconds;
        if strcmp(d.classification.outcome,'B'),B=mid;elseif strcmp(d.classification.outcome,'C'),C=mid;
        else,error('BoundedTransfer:PilotObstruction','Pilot bracket obstructed by %s at step %d.',d.classification.outcome,step);end
    end
    % Signed edge perturbations independently verify B/edge/C organization.
    fl=floquets{3};idx=fl.unstable_indices(1);x=edge.state(:,1);w=weights(x);
    vec=real(fl.weighted_eigenvectors(:,idx))./sqrt([w;w]);vec=vec/norm(sqrt([w;w]).*vec);
    normstate=norm(sqrt([w;w]).*[edge.state(:,2);edge.state(:,3)]);outcomes=cell(1,2);
    for signindex=1:2
        sg=[-1,1];st=edge.state;st(:,2:3)=st(:,2:3)+sg(signindex)*1e-4*normstate*reshape(vec,[],2);
        d=label_cached(fd,sprintf('edge_sign_%d',signindex),st,library,cfg);outcomes{signindex}=d.classification.outcome;
        pilot_labels=pilot_labels+1;pilot_runtime=pilot_runtime+d.runtime_seconds;
    end
    assert(all(ismember({'B','C'},outcomes)),'Target edge signed-branch gate failed.');
end
x=B(:,1);w=weights(x);N=numel(x);sw=sqrt([w;w]);
delta=[C(:,2)-B(:,2);C(:,3)-B(:,3)].*sw;len=norm(delta);e1=delta/len;
raw=[cos(3*pi*x/cfg.grid.Lx);-0.5*cos(2*pi*x/cfg.grid.Lx)].*sw;
e2=raw-e1*(e1'*raw);e2=e2/norm(e2);basis=[e1,e2]./sw;
center=interpolate(B,C);h=[2/3,.2]*len;
[xx,yy]=meshgrid(linspace(-1,1,17));Z=[reshape(xx',[],1),reshape(yy',[],1)];
states=zeros(N,3,289);
for k=1:289
    states(:,:,k)=center;states(:,2:3,k)=center(:,2:3)+reshape(basis*(h(:).*Z(k,:)'),N,2);
end
assert(min(states(:,2:3,:),[],'all')>1e-6,'Positivity gate failed.');
% Geometry is frozen without querying candidate labels or using edge coordinates.
save(fullfile(fd,'geometry.mat'),'states','Z','center','basis','h','cfg','library','edge','-v7.3');
T=table((1:289)',compose('state_%03d',(1:289)'),Z(:,1),Z(:,2), ...
    'VariableNames',{'state_index','sample_id','normalized_alpha_1','normalized_alpha_2'});
writetable(T,fullfile(fd,'geometry.csv'));
writejson(fullfile(fd,'preparation_summary.json'),struct('family',family,'b',cfg.edge_tracking.frozen_b, ...
    'status','passed','pilot_labels',pilot_labels,'pilot_label_runtime_seconds',pilot_runtime, ...
    'setup_wall_seconds',toc(timer),'chord_physical_length',len,'halfwidths',h, ...
    'candidate_labels_used',false,'exact_edge_coordinates_used_for_geometry',false, ...
    'exact_edge_used_only_for_parameter_validation_and_recovery',true));
end

function r=reference(o,label)
phase=(0:719)'/720;theta=o.orbit_times/o.period;
r=struct('label',label,'type','periodic','period',o.period, ...
    'U',interp1(theta,o.orbit_U,phase,'pchip'),'V',interp1(theta,o.orbit_V,phase,'pchip'));
end
function s=state_at(G,id),s=G.states(:,:,id);end
function s=interpolate(B,C),s=B;s(:,2:3)=(B(:,2:3)+C(:,2:3))/2;end
function w=weights(x),w=ones(size(x))*(x(2)-x(1));w([1,end])=w([1,end])/2;end
function d=physical_separation(B,C)
w=weights(B(:,1));d=sqrt(sum(w.*sum((B(:,2:3)-C(:,2:3)).^2,2)));
end
function data=label_cached(fd,id,state,library,cfg)
raw=fullfile(fd,'labels');if ~isfolder(raw),mkdir(raw);end
file=fullfile(raw,[id,'.mat']);
if isfile(file)
    q=load(file);data=q.case_data;assert(max(abs(data.initial_state-state),[],'all')<1e-12,'Cached state mismatch.');return;
end
case_data=brusselator_evolve_and_classify_state(id,state,library,cfg);save(file,'case_data','-v7.3');data=case_data;
writejson(fullfile(raw,[id,'.json']),struct('case_id',id,'final_label',data.classification.outcome, ...
    'duration',data.duration,'runtime_seconds',data.runtime_seconds,'decision_time',data.decision_time));
end
function writejson(file,s)
f=fopen(file,'w');assert(f>=0);fprintf(f,'%s\n',jsonencode(s,PrettyPrint=true));fclose(f);
end
