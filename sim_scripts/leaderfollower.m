%%%INPUTS: 
nu_cell = {0,0}; include_spec=[1 1]; v_opt = 2;

f_xv_true_cell{1} = [];
f_xv_true_cell{2} = @(x,v,t) 15*(1+2/3*cos(2*v)).*(exp(-x/0.05) - 0.25*exp(-x/0.1));

h_xv_true_cell{1}  = [];
h_xv_true_cell{2}  = @(x,v,t) -8*exp(-8*x).*(1+cos(v));

d_xv_true_cell{1}  = [];
d_xv_true_cell{2}  = @(v,xv,t) -5*v;

opts_cell={[1 0 0],[1 0 0],[1 0 0]};

N = 500;
d = 2;
dt = 8/60/10;
subdt = 8;
M = 4000;
verbose = 1;

xmus = 0*ones(1,d)+4;
xsigs = 0.5*eye(size(xmus,2));
x0_dist = 3;

if v_opt == 1
	vmus= [-0.0205    0.0346];                                                                                         
	vsigs =   1.0e-03*[[0.1393    0.0502];[0.0502    0.1212]]; 
elseif v_opt == 2
	vmus= [0 0.01];                                                                                 
	vsigs =   [[0 0];[0 0]]; 
end

v0_dist = 0;

bdry = inf*[[-1 1];[-1 1]];
per=range(imag(bdry),2);

inc_frac = 0.2;
inc = floor(inc_frac*N);
inds_cell_true = {1:floor(inc_frac*N),floor(inc_frac*N)+1:N};

Xtemp = gen_X0(xmus,xsigs,N,x0_dist);
Vtemp = gen_X0(vmus,vsigs,N,v0_dist);
Vtemp(inds_cell_true{2},:) = 0;

RHS = @(Z,t) rhs_lean_het(f_xv_true_cell, h_xv_true_cell, d_xv_true_cell, Z, inds_cell_true, d, t, opts_cell, per);

t = 0:dt:M*dt;

rng('shuffle');
rng_seed = rng().Seed; 
rng(rng_seed);

Ztemp = reshape([Xtemp Vtemp],[],1);
Ntot = length(Ztemp);
Z = zeros(M,Ntot);
Z(1,:) = Ztemp;

m=subdt;
for mm=2:length(t)
    if verbose>0
        tic;
    end
    t_temp = t(mm-1)+dt/m:dt/m:t(mm);
    for j=1:length(t_temp)
        Ztemp =Ztemp+dt/m*RHS(Ztemp,t_temp(j));
        for k=1:length(inds_cell_true)
            nu=nu_cell{k};
            inds=inds_cell_true{k};
            Nfoc=length(inds);
            if nu>0
                Btemp = sqrt(2*nu*(dt/m))*randn(Nfoc*2,1);
                Ztemp([2*N+inds(:);3*N+inds(:)]) = Ztemp([2*N+inds(:);3*N+inds(:)]) + Btemp;
            end
        end
        xv = Ztemp(end/2+1:3*end/4);
        yv = Ztemp(3*end/4+1:end);
        if per(1) == 0
            indsl=Ztemp(1:end/4)<bdry(1,1);
            indsr=Ztemp(1:end/4)>bdry(1,2);
            x = min(max(bdry(1,1),Ztemp(1:end/4)),bdry(1,2)); 
            xv(indsl) = max(0,xv(indsl));
            xv(indsr) = min(0,xv(indsr));
        elseif per(1) == 1
            x = bdry(1,1)+mod(Ztemp(1:end/4)-bdry(1,1),bdry(1,2)-bdry(1,1)); 
        end
        if per(2) == 0
            indsl=Ztemp(end/4+1:end/2)<bdry(2,1);
            indsr=Ztemp(end/4+1:end/2)>bdry(2,2);
            y = min(max(bdry(2,1),Ztemp(end/4+1:end/2)),bdry(2,2));
            yv(indsl) = max(0,yv(indsl));
            yv(indsr) = min(0,yv(indsr));
        elseif per(2) == 1
            y = bdry(2,1)+mod(Ztemp(end/4+1:end/2)-bdry(2,1),bdry(2,2)-bdry(2,1)); 
        end            
        Ztemp=[x;y;xv;yv];
    end
    Z(mm,:) = Ztemp;
    if verbose>0
        disp([mm toc]);
    end
end

Xscell_temp = zeros(N,d,length(t)); Vscell_temp = zeros(N,d,length(t));
for m=1:length(t)
    Xscell_temp(:,:,m) = reshape(Z(m,1:end/2)',N,d);
    Vscell_temp(:,:,m) = reshape(Z(m,end/2+1:end)',N,d);
end
Xscell = {Xscell_temp}; 
Vscell = {Vscell_temp};

%% View simulation

runs=1;
toggle_plot = 10;%ceil(M/200);
io=1;iend=4000;
memlength=inf;%40*toggle_plot;
comp_inds=[1];
if toggle_plot>0
    clf;
    for nn=1:runs
        X = Xscell{nn};
        V = Vscell{nn};
        Xscell_comp=X(comp_inds,:,:);
        Xmin = min(min(X(:,1,:)));
        Xmax = max(max(X(:,1,:)));
        Ymin = min(min(X(:,2,:)));
        Ymax = max(max(X(:,2,:)));
        for m=max(1,io):toggle_plot:min(length(t),iend)
            for j=1:length(inds_cell_true)
                Xmean=mean(reshape(permute(X(inds_cell_true{j},:,m),[1 3 2]),[],2));
                plot(X(inds_cell_true{j},1,m),X(inds_cell_true{j},2,m),'o')
                hold on
            end
            % Xmean=mean(reshape(permute(X(:,:,m),[1 3 2]),[],2));
            % plot(X(:,1,m),X(:,2,m),'o')
            hold on
            if ~isempty(comp_inds)
                if m>memlength
                    plot(squeeze(Xscell_comp(:,1,m-memlength:toggle_plot:m)),squeeze(Xscell_comp(:,2,m-memlength:toggle_plot:m)),'g.','linewidth',0.5,'markersize',5)
                elseif memlength==inf
                    plot(squeeze(Xscell_comp(:,1,1:toggle_plot:m)),squeeze(Xscell_comp(:,2,1:toggle_plot:m)),'g.','linewidth',0.5,'markersize',5)
                end
            end
            hold off
            title(num2str([m mean(vecnorm(V(:,:,m),2,2))]))
            axis equal
            xlim([Xmin Xmax])
            ylim([Ymin Ymax])
            drawnow
        end
    end
end

%% Functions

function X0 = gen_X0(mus,sigs,N,x0_dist)
    [m,dim] = size(mus);
    if length(sigs(:))==1
        sigs = repmat(sigs*eye(dim),m,1);
    end
    X0 = zeros(N,dim);
    inc = floor(N/m);
    for mm=1:m-1             
        if x0_dist == 0
            X0( (mm-1)*inc+(1:inc),:) = mvnrnd(mus(mm,:),sigs(1+(mm-1)*dim:mm*dim,:),inc);
        elseif x0_dist == 1
            X0( (mm-1)*inc+(1:inc),:) = (2*rand(inc,dim)-1)*sigs(1+(mm-1)*dim:mm*dim,:)+mus(mm,:);
        elseif x0_dist == 2
            X0( (mm-1)*inc+(1:inc),:) = samplellipse(inc,sigs(1+(mm-1)*dim,1),sigs(mm*dim,2),mus(mm,1),mus(mm,2));
        elseif x0_dist == 3
            X0( (mm-1)*inc+(1:inc),:) = (2*lhsdesign(inc,dim)-1)*sigs(1+(mm-1)*dim:mm*dim,:)+mus(mm,:);
        end
    end
    inds = (m-1)*inc+1:N;
    if x0_dist == 0
        X0(inds,:) = mvnrnd(mus(m,:),sigs((m-1)*dim+1:end,:),length(inds));
    elseif x0_dist == 1
        X0(inds,:) = (2*rand(length(inds),dim)-1)*sigs(1+(m-1)*dim:end,:)+mus(m,:);
    elseif x0_dist == 2
        X0(inds,:) = samplellipse(length(inds),sigs(1+(m-1)*dim,1),sigs(end,2),mus(m,1),mus(m,2));
    elseif x0_dist == 3
        X0(inds,:) = (2*lhsdesign(inc,dim)-1)*sigs(1+(m-1)*dim:end,:)+mus(m,:);
    end
end

function X = samplellipse(N,rx,ry,ax,ay)
    init_r = rand(N,1); 
    init_theta = 2*pi*rand(N,1); 
    X = [rx*sqrt(init_r).*cos(init_theta)+ax ry*sqrt(init_r).*sin(init_theta)+ay];
end

function dZ = rhs_lean_het(f_xv_true_cell, h_xv_true_cell, d_xv_true_cell, Ztemp, inds_cell, d, t, opts_cell, per)
    Ztemp = reshape(Ztemp,[],2*d);
    agents = Ztemp(:,1:2);
    vels = Ztemp(:,3:4);
    updateV = vels*0;    
	[N,d] = size(agents);
    num_pop = length(inds_cell);
    
    for i=1:num_pop
        inds=inds_cell{i};
        f_xv_true=f_xv_true_cell{i};
        h_xv_true=h_xv_true_cell{i};
        d_xv_true=d_xv_true_cell{i};
        opts=opts_cell{i};
        if isempty(opts)
            opts = zeros(3,1);
        end
        Nfoc=length(inds);
        dX = reshape(agents(inds,:),Nfoc,1,d) - reshape(agents,1,N,d);
        if ~isempty(per)
            dX = periodize(dX,per,3);
        end
        rr = squeeze(vecnorm(dX,2,3));

        dV = repmat(reshape(vels(inds,:),Nfoc,1,d),1,N,1);    
        vv = dot(dV,-dX,3);
        vv = atan2(squeeze(dV(:,:,2).*dX(:,:,1)-dV(:,:,1).*dX(:,:,2)),vv);

        if ~isempty(f_xv_true)
            forces = f_xv_true(rr,vv,t);
            if opts(1)==1
                updateV(inds,:) = updateV(inds,:) + reshape(mean((forces./max(squeeze(vecnorm(dX,2,3)),eps)).*dX,2),Nfoc,d);
            else
                updateV(inds,:)  = updateV(inds,:)  + reshape(mean(dX.*forces,2),Nfoc,d);
            end
        end

        if ~isempty(h_xv_true)
            forces = h_xv_true(rr,vv,t);
            dV = reshape(vels(inds,:),Nfoc,1,d) - reshape(vels,1,N,d);
            if opts(2)==1
                updateV(inds,:)  = updateV(inds,:)  + reshape(mean((forces./max(squeeze(vecnorm(dV,2,3)),eps)).*dV,2),Nfoc,d);
            else
                updateV(inds,:)  = updateV(inds,:)  + reshape(mean(forces.*dV,2),Nfoc,d);
            end
        end

        if ~isempty(d_xv_true)
            rr = repmat(vecnorm(vels(inds,:),2,2),1,N);
            forces = d_xv_true(rr,vv,t);
            dV = reshape(vels(inds,:),[],1,d);
            if opts(3)==1
                updateV(inds,:)  = updateV(inds,:)  + reshape(mean(forces.*dV./max(squeeze(vecnorm(vels,2,3)),eps),2),Nfoc,d);
            else
                updateV(inds,:)  = updateV(inds,:)  + reshape(mean(forces.*dV,2),Nfoc,d);
            end
        end
    end
    dZ = reshape([vels updateV],[],1);    
end

function X = periodize(X,L,dim)
    dims=1:length(size(X));
    L = permute(L(:),dims([dim 2:dim-1 1 dim+1:end])); 
    X = X + L.*max(sign(abs(X)-L/2),0).*sign(X);
end
