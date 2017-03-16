%%%% dmd of noaa data, unevenly spaced

clear all

filename = 'sst_uneven_new.mat'

%load('~/Dropbox/data/noaa-prepped/sst_fig_data.mat','data2','indices1')
load('~/Dropbox/data/noaa-prepped/sst_fig_data_big.mat','data2','indices2')

t_u = indices2;

%% 

nw = 521;

data2 = data2(:,1:nw);
t_u = indices2(1:nw);
[mx,nx] = size(data2);

%% 

tic
[u,sd_u,~] = svd(data2,'econ');
timesvd_u = toc

sd_u = diag(sd_u);

beta = min(nx*1.0/mx,mx*1.0/nx);
gdfac = optimal_SVHT_coef(beta,0);
            
rgd_u = sum(sd_u > gdfac*median(sd_u))

u = u(:,1:rgd_u);

%% first rank setting

r1_u = 16;
r = r1_u;

% optimized dmd

maxiter = 30;
tol = 1e-6;
opts = varpro_opts('maxiter',maxiter,'tol',tol,'eps_stall',1e-6);

imode = 2; % projected version
tic
[woptsave_u,eoptsave_u,boptsave_u] = optdmd(data2,t_u,r,imode,opts,[],u(:,1:r));
timeoptsave_u = toc

[~,rerroptsave_u] = best_reconstruction(data2,eoptsave_u,t_u);


mx_u = mx;
nx_u = nx;

%%

save(filename,'r1_u','rgd_u','sd_u','mx_u','nx_u', ...
    'woptsave_u','eoptsave_u','boptsave_u','timeoptsave_u',...
    'rerroptsave_u','t_u','timesvd_u')
