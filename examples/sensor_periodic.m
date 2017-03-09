%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Generate data for peirodic system example 
% with additive sensor noise
%
% Run time:  6 mins 40 s on Intel i7 laptop
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialize

% Seed random number generator
iseed = 8675309;
rng(iseed); % consistent runs
%rng('shuffle'); % based on current time

% test settings

% range of number of snapshots
moffset = 5;
nm = 8;
maxm = 2^(moffset+nm);

nsigs = 8; % noise levels
ntrials = 100; % number of trials for each

% system settings

dt = 0.1;
mrank = 2;

Xclean = zeros(2,maxm);

% system matrix

A = [1 -2; 1 -1];
Admd = expm(A*dt);

[abasis,evals] = eig(A);

Xclean(:,1) = [1; .1];

Xcoeffs = abasis\Xclean(:,1);

% clean dynamics

for i = 2:maxm
    Xclean(:,i) = abasis*diag(exp(diag(evals*dt*(i-1))))*Xcoeffs;
end

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

asave1 = zeros(mrank,mrank,ntrials,nsigs,nm);
asave2 = zeros(mrank,mrank,ntrials,nsigs,nm);
asave3 = zeros(mrank,mrank,ntrials,nsigs,nm);
asave4 = zeros(mrank,mrank,ntrials,nsigs,nm);

for i = 1:nm
    m = 2^(moffset+i);
    ms(i) = m;
    t = (0:(m-1))*dt;
    for iii = 1:nsigs
        sigma = sqrt(10^(-(1+2*(iii-1))));
        sigs(iii) = sigma;
        for jjj = 1:ntrials

            xrand = randn(size(Xclean(:,1:m)));
            X = Xclean(:,1:m) + sigma*xrand;

            r = mrank;
            
            %% exact dmd 

            imode = 1;
            [w,e,atilde,afull] = dmdexact(X,r,imode);

            e = log(e)/dt;

            indices = match_vectors(e,evals(1:mrank));
            
            eigsave1(1:mrank,jjj,iii,i) = e(indices);
            asave1(:,:,jjj,iii,i) = afull;
            
            %% forward-backward dmd 

            imode = 2; % modified version ...
            [w,e,atilde] = fbdmd(X,r,imode);
            
            e = log(e)/dt;

            indices = match_vectors(e,evals(1:mrank));
            eigsave2(1:mrank,jjj,iii,i) = e(indices);
            asave2(:,:,jjj,iii,i) = atilde;
            
            %% total-least squares dmd 

            imode = 2; % projected version (see dawson et. al)
            [w,e,atilde,afull] = tlsdmd(X,r,imode);
            e = log(e)/dt;

            indices = match_vectors(e,evals(1:mrank));

            eigsave3(1:mrank,jjj,iii,i) = e(indices);
            asave3(:,:,jjj,iii,i) = afull;            

            %% optdmd

            maxiter = 30;
            tol = sigma/10;

            opts = varpro_opts('maxiter',maxiter,'tol',tol,'eps_stall',1e-9);
            
            imode = 2; % projected version
            [w,e,b,atilde,~,afull] = optdmd(X,t,r,imode,opts);

            indices = match_vectors(e,evals(1:mrank));

            eigsave4(1:mrank,jjj,iii,i) = e(indices);
            asave4(:,:,jjj,iii,i) = afull;

            sigs(iii,jjj,i) = sigma;

            %s = sprintf('iii %i jjj %i',iii,jjj);
            %disp(s)

        end
    end
end

filename = 'sensor_periodic.mat';
save(filename,'eigsave1','eigsave2','eigsave3','eigsave4', ...
    'asave1','asave2','asave3','asave4','sigs','ms','A','Admd','dt', ...
    'Xclean','evals');
