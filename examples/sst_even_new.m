%%%% dmd of noaa data, equispaced

clear all

filename = 'sst_even_new.mat'

%load('~/Dropbox/data/noaa-prepped/sst_fig_data.mat','data1','indices1')
load('~/Dropbox/data/noaa-prepped/sst_fig_data_big.mat','data1','indices1')

t = indices1;

%% 

nw = 521;

data1 = data1(:,1:nw);
t = indices1(1:nw);
[mx,nx] = size(data1);

%% 

dt = t(2)-t(1);

tic
[u,sd,~] = svd(data1,'econ');
timesvd = toc

sd = diag(sd);

beta = min(nx*1.0/mx,mx*1.0/nx);
gdfac = optimal_SVHT_coef(beta,0);
            
rgd = sum(sd > gdfac*median(sd))

rs = [2, 4, 8, 16 ,32, 64, 128];

rtot = max(max(rs),rgd);

u = u(:,1:rtot);

nr = length(rs);

esexact = zeros(rtot,nr);
estls = zeros(rtot,nr);
esopt = zeros(rtot,nr);
timesexact = zeros(nr,1);
timestls = zeros(nr,1);
timesopt = zeros(nr,1);
rerrsexact = zeros(nr,1);
rerrstls = zeros(nr,1);
rerrsopt = zeros(nr,1);

%%

for i = 1:length(rs)

    r = rs(i)

    % exact dmd

    imode = 2;
    tic
    [w,e,~] = dmdexact(data1,r,imode,u(:,1:r));
    timeexact = toc
    e = log(e)/dt;

    tic
    [~,rerr] = best_reconstruction(data1,e,t);
    timer = toc
    
    esexact(1:r,i) = e;
    timesexact(i) = timeexact;
    rerrsexact(i) = rerr;

    % total least squares dmd

    imode = 2; % projected version (see dawson et. al)
    tic
    [w,e,~] = tlsdmd(data1,r,imode,u(:,1:r));
    timetls = toc
    e = log(e)/dt;

    tic
    [~,rerr] = best_reconstruction(data1,e,t);
    timer = toc
    
    estls(1:r,i) = e;
    timestls(i) = timetls;
    rerrstls(i) = rerr;


    % optimized dmd

    maxiter = 30;
    tol = 1e-6;
    opts = varpro_opts('maxiter',maxiter,'tol',tol,'eps_stall',1e-6);

    imode = 2; % projected version
    tic
    [w,e,b] = optdmd(data1,t,r,imode,opts,[],u(:,1:r));
    timeopt = toc
    
    tic
    [~,rerr] = best_reconstruction(data1,e,t);
    timer = toc
    
    esopt(1:r,i) = e;
    timesopt(i) = timeopt;
    rerrsopt(i) = rerr;
    
    if (r == 16)
        eoptsave = e;
        woptsave = w;
        boptsave = b;
    end

end
%%

save(filename,'rs','rgd','sd','mx','nx','dt', ...
    'esexact','timesexact',...
    'estls','timestls',...
    'esopt','timesopt',...
    'rerrsexact','rerrstls','rerrsopt',...
    'timesvd','woptsave','eoptsave','boptsave');
