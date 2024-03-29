%% choose experiment and load variables
clear all; clc

exper = '111';

data_dr = ['~/Desktop/data_dr/',exper,'/'];
save_dr = data_dr;
input_data = findfilestrloc(save_dr,'sim',1);
load([save_dr,findfile(save_dr,'classify_',[])]);
load([save_dr,input_data],'f_xv_true_cell','h_xv_true_cell','d_xv_true_cell','inds_cell_true');
ftrue_cell = {f_xv_true_cell,h_xv_true_cell,d_xv_true_cell};
load([save_dr,findfile(save_dr,'singlecell_',[])],'tobs','ninds','Xscell_obs','Vscell_obs');
tcap=1; test_tinds_frac= 0.25;
errfun = @(X,Y,V,W) norm(reshape(V(:,:,1:floor(tcap*end))-W(:,:,1:floor(tcap*end)),[],1))/...
    norm(reshape(W(:,:,1:floor(tcap*end)),[],1));
test_tinds = 1:floor(length(tobs)*test_tinds_frac);

% choose cluster and force to view

for cluster_num = 1:length(find(cellfun(@(x) ~isempty(x), species_models)))

    disp(['----cluster ', num2str(cluster_num), '----'])

    for force_ind = 1:3
    
        % set plotting arguments
        
        n = 1000;
        if force_ind<3
            r = 2;
            v =[0 1];
            [rr,th,xx,yy] = build_polar_grid(n,r,v);
        else
            r = max(reshape(vecnorm(Vscell_obs{1},2,2),[],1));
            [rr,th] = meshgrid(linspace(0,r,n),linspace(0,pi,n));
            xx=rr;
            yy=th;
        end
        x = xx(1,:);
        y = yy(:,1);
        
        % evaluate true force on grid
        
        forceDNA = cellfun(@(x) ~isempty(x), species_models{cluster_num}(2:4));
        if isequal(forceDNA,[1 1 1])
            learned_spec = 1;
        elseif isequal(forceDNA,[1 0 1])
            learned_spec = 2;
        elseif isequal(forceDNA,[0 1 1])
            learned_spec = 3;
        end
        F_true = ftrue_cell{force_ind}{learned_spec};
        if isempty(F_true)
            F_true = rr*0;
        else
            F_true = F_true(rr,th,0);
        end
        
        % evaluate learned force on grid
        
        f_learned = species_models{cluster_num}{force_ind+1};
        if ~isempty(f_learned) 
            F_learned = f_learned(rr,th);
        else
            F_learned = 0*rr;
        end
        
        force_str = {'a-r','align','drag'};
        disp(['L2 ',force_str{force_ind},' force error = ',num2str(norm(F_learned(:)-F_true(:))/norm(F_true(:)))])
    
    end

        
    CS_A = length(intersect(species_inds{cluster_num},inds_cell_true{1}))/length(intersect(inds_cell_true{1},ninds));
    CS_B = length(intersect(species_inds{cluster_num},inds_cell_true{2}))/length(intersect(inds_cell_true{2},ninds));
    CS_C = length(intersect(species_inds{cluster_num},inds_cell_true{3}))/length(intersect(inds_cell_true{3},ninds));

    disp(['CS(A) CS(B) CS(C) = ',num2str([CS_A CS_B CS_C])])

    subinds = oppinds(cell2mat(species_inds(1:cluster_num-1)),length(ninds));
    valid_cells=ninds(subinds);
    member_inds = ismember(subinds,species_inds{cluster_num});
    valpairs = valpairs_cell{cluster_num}(member_inds);
    relerr = mean(compute_errs(valid_cells(member_inds),valpairs,{Xscell_obs{1}(:,:,test_tinds)},{Vscell_obs{1}(:,:,test_tinds)},[],errfun));
    disp(['Delta V = ',num2str(relerr)])

end