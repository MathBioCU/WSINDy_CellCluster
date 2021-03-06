%% Get single-cell models (local machine)

data_dr = '~/Desktop/data_dr/111/';

save_dr = '~/Desktop/save_dr/';

input_data='sim08-Feb-2022_000-111_1_comb4000.mat';

load([data_dr,input_data],'Xscell','Vscell','t')

singlecell_inputs;
precomp_learningenvironment;

ninds = home_cell{expr};
X = Xscell_obs{expr};
V = Vscell_obs{expr};
algout=cell(length(ninds),1);
simdat=cell(length(ninds),1);
neighbs=cell(length(ninds),1);

parfor i1=1:length(ninds)

    ind_cell = ninds(i1);

    tic;
    %%% Build linear system
    [Gs,bs,Cfs] = wsindy_sde2nd_fun_nonspatial(X,V,tobs,pt,mt,1,ind_cell,neighbs_cell{expr},tinds,fx_fcn_cell,fv_fcn_cell,hx_fcn_cell,hv_fcn_cell,dx_fcn_cell,dv_fcn_cell);    
    
    %%% Solve sparse regression problem
    [W,its,lossvals,lambda_hat] = sparsifyDynamics_seq(lambda,gamma*norm(Gs.*(M')),[Gs.*(M');G_append],[bs;b_append],M,Aineq.*(M'),bineq,Aeq.*(M'),beq,excl_inds,const_tol,max_its_qp,disp_opt,max_its_stls,alpha);
    resid = (bs-Gs*W)/norm(bs);
        
    %%% Compute learned forces
    [f_learned,h_learned,d_learned] = gen_force_fcn_xv(W,fx_fcn_cell,fv_fcn_cell,hx_fcn_cell,hv_fcn_cell,dx_fcn_cell,dv_fcn_cell);
    
    %%% Compute learned diffusion - not utilized in study
%     nu_learned = rms(Gs*W-bs)^2/2;
    
    %%% Compute neighbors of homing cell
    valid_cells = get_neighbs(X, ind_cell, home_cell, V, expr, nearestKLLneighbs, alphaKL, knnp1, tobs, numbins);
    
    %%% Validate model on neighbor cells
    [Xspred,errs] = test_neighbs(X,V,tobs,0,nufac_x,nufac_v,f_learned,h_learned,d_learned,valid_cells,opts,subdt,avg_v0,test_tinds_frac,verbose,0);
    disp(['Time to learn single-cell model: ',num2str(toc), 's; ind: ',num2str(i1)])

    %%% Save single-cell model
    algouttemp={f_learned,h_learned,d_learned,W,resid,lossvals,lambda_hat};

    algout{i1}=algouttemp;
    simdat{i1}={Xspred,errs};
    neighbs{i1}=valid_cells;

end

% poolobj = gcp('nocreate');
% delete(poolobj)

clear Xscell Vscell Vnormcol
consol_data = [save_dr,'singlecell_',input_data];
save(consol_data)

%% classify

single_cell_data = [save_dr,'singlecell_',input_data];
load(single_cell_data,'fx_fcn_cell','fv_fcn_cell','hx_fcn_cell','hv_fcn_cell','dx_fcn_cell','dv_fcn_cell','input_data','ninds','Xscell_obs','Vscell_obs','tobs','neighbs','simdat','algout');

classify_inputs;
classify_models;

% poolobj = gcp('nocreate');
% delete(poolobj)

classified_data = [save_dr,'classify_',input_data];
save(classified_data,'species_inds','species_models','valpairs_cell');

%% Generate trajectories from learned models

load([save_dr,'singlecell_',input_data])
load([save_dr,'classify_',input_data])

specs=1:sum(cellfun(@(x) ~isempty(x),species_models));
X_preds = cell(length(specs),1);

for spec = 1:length(specs)

    f_learned = species_models{specs(spec)}{2};
    h_learned = species_models{specs(spec)}{3};
    d_learned = species_models{specs(spec)}{4};

    [~,a]=sort(cellfun(@(x) x{2}(1),simdat(species_inds{spec})));
    a=a(min(2,end));
    M=dist(squeeze(X(ninds(species_inds{specs(spec)}),:,1))');
    [b,cell_inds]=sort(M(a,:));
    valid_cells = ninds(species_inds{specs(spec)}(cell_inds));

    [X_pred,~] = test_neighbs(X,V,tobs,nu_learned,nufac_x,nufac_v,f_learned,h_learned,d_learned,valid_cells,opts,subdt,avg_v0,1,2,1);
    X_preds{spec} = {X_pred,valid_cells};

end

save([save_dr,'learnedtraj_',input_data]);