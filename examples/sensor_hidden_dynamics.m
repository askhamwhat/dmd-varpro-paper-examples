%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Generate data for system with 
% both growing and decaying dynamics,
% with additive sensor noise
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialize

% Seed random number generator
iseed = 8675309;
rng(iseed); % consistent runs
%rng('shuffle'); % based on current time

% test settings

% range of number of snapshots
moffset = 6;
nm = 3;
maxm = 2^(moffset+nm);

nsigs = 8; % noise levels
ntrials = 100; % number of trials for each

% system settings

nx = 300;
x = linspace(0,15,nx);

t = linspace(0,2*pi,maxm);
dt = t(2)-t(1);

Xclean = zeros(nx,maxm);

% system matrix

gamma1 = 1;
gamma2 = -.2;
k1 = 1;
k2 = .4;
omega1 = 1;
omega2 = 3.7;

evals = [gamma1+1i*omega1; gamma1-1i*omega1; gamma2+1i*omega2; ...
    gamma2-1i*omega2];

mrank = 4;

% clean dynamics

[tt,xx] = meshgrid(t,x);
Xclean = sin(k1*xx-omega1*tt).*exp(gamma1*tt) + ...
    sin(k2*xx-omega2*tt).*exp(gamma2*tt);

%% storage

sigs = zeros(nsigs,1);
ms = zeros(nm,1);

% save eigenvalues as computed 
% by various methods

% 1 - exact dmd
% 2 - forward-backward dmd
% 3 - total least squares dmd
% 4 - dmdu

eigsave1 = zeros(mrank,ntrials,nsigs,nm);
eigsave2 = zeros(mrank,ntrials,nsigs,nm);
eigsave3 = zeros(mrank,ntrials,nsigs,nm);
eigsave4 = zeros(mrank,ntrials,nsigs,nm);

rerrsave1 = zeros(ntrials,nsigs,nm);
rerrsave2 = zeros(ntrials,nsigs,nm);
rerrsave3 = zeros(ntrials,nsigs,nm);
rerrsave4 = zeros(ntrials,nsigs,nm);


for i = 1:nm
    m = 2^(moffset+i);
    tdmd = t(1:m);
    ms(i) = m;
    for iii = 1:nsigs
        sigma = 2^(-iii);
        sigs(iii) = sigma;
        parfor jjj = 1:ntrials

            xrand = randn(size(Xclean(:,1:m)));
            X = Xclean(:,1:m) + sigma*xrand;

            r = mrank;
            
            [u,~,~] = svd(X,'econ');
            u = u(:,1:r);
            
            %% exact dmd 

            imode = 1;
            [w,e,atilde] = dmdexact(X,r,imode,u);

            e = log(e)/dt;
            
            [Xr,rerr] = best_reconstruction(X,e,tdmd);

            indices = match_vectors(e,evals(1:mrank));
            
            eigsave1(:,jjj,iii,i) = e(indices);
            rerrsave1(jjj,iii,i) = rerr;

            %% forward-backward dmd 

            imode = 3; % modified version ...
            [w,e,atilde] = fbdmd(X,r,imode,u);

            e = log(e)/dt;
            
            [Xr,rerr] = best_reconstruction(X,e,tdmd);

            indices = match_vectors(e,evals(1:mrank));
            eigsave2(:,jjj,iii,i) = e(indices);
            rerrsave2(jjj,iii,i) = rerr;

            %% total-least squares dmd 

            imode = 2; % projected version (see dawson et. al)
            [w,e,atilde] = tlsdmd(X,r,imode,u);
            e = log(e)/dt;
            
            [Xr,rerr] = best_reconstruction(X,e,tdmd);

            indices = match_vectors(e,evals(1:mrank));
            eigsave3(:,jjj,iii,i) = e(indices);
            rerrsave3(jjj,iii,i) = rerr;

            %% optdmd

            maxiter = 30;
            tol = sigma/10;

            opts = varpro_opts('maxiter',maxiter,'tol',tol,'eps_stall',1e-9);
            
            imode = 2; % projected version
            
            [w,e,b,atilde] = optdmd(X,tdmd,r,imode,opts,[],u);
            
            [Xr,rerr] = best_reconstruction(X,e,tdmd);

            indices = match_vectors(e,evals(1:mrank));

            eigsave4(:,jjj,iii,i) = e(indices);
            rerrsave4(jjj,iii,i) = rerr;

            sigs(iii,jjj,i) = sigma;

            %s = sprintf('iii %i jjj %i',iii,jjj);
            %disp(s)
            %fprintf('iii %i jjj %i i %i\n',iii,jjj,i)
            fprintf('%e percent complete\n',...
                100*((i-1)*nsigs*ntrials+(iii-1)*ntrials+jjj)/...
                (1.0*nsigs*ntrials*nm))

        end
    end
end

%% 

filename = 'sensor_hidden_dynamics.mat';
save(filename,'eigsave1','eigsave2','eigsave3','eigsave4', ...
    'sigs','ms','dt','Xclean','evals','tt','xx','rerrsave1',...
    'rerrsave2','rerrsave3','rerrsave4');
