%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Run calculations similar to those
% in sensor_hidden_dynamics.m but
% in a way that's suitable for estimating
% timing
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

nsigs = 6; % noise levels
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
% 4 - optimized dmd

timesave1 = zeros(nsigs,nm);
timesave2 = zeros(nsigs,nm);
timesave3 = zeros(nsigs,nm);
timesave4 = zeros(nsigs,nm);

for i = 1:nm
    m = 2^(moffset+i);
    tdmd = t(1:m);
    ms(i) = m;
    for iii = 1:nsigs
        
        sigma = 2^(-iii);
        sigs(iii) = sigma;
        
        r = mrank;
        
        tic
        for j = 1:ntrials
            
            %% exact dmd 

            xrand = randn(size(Xclean(:,1:m)));
            X = Xclean(:,1:m) + sigma*xrand;
            
            imode = 1;
            [w,e,atilde] = dmdexact(X,r,imode);

        end
        temp = toc;
        
        timesave1(iii,i) = temp/ntrials;

        tic
        for j = 1:ntrials

            xrand = randn(size(Xclean(:,1:m)));
            X = Xclean(:,1:m) + sigma*xrand;

            %% forward-backward dmd 

            imode = 3; % modified version ...
            [w,e,atilde] = fbdmd(X,r,imode);

        end
        temp = toc;
        
        timesave2(iii,i) = temp/ntrials;
        
        tic
        for j = 1:ntrials

            xrand = randn(size(Xclean(:,1:m)));
            X = Xclean(:,1:m) + sigma*xrand;

            %% total-least squares dmd 

            imode = 2; % projected version (see dawson et. al)
            [w,e,atilde] = tlsdmd(X,r,imode);

        end
        temp = toc;
        
        timesave3(iii,i) = temp/ntrials;
        
        tic
        for j = 1:ntrials

            xrand = randn(size(Xclean(:,1:m)));
            X = Xclean(:,1:m) + sigma*xrand;

            %% optdmd

            maxiter = 30;
            tol = sigma/10;

            opts = varpro_opts('maxiter',maxiter,'tol',tol,'eps_stall',1e-9);
            imode = 2; % projected version
            [w,e,b,atilde] = optdmd(X,tdmd,r,imode,opts);

        end
        temp = toc;
        
        timesave4(iii,i) = temp/ntrials;

    end
end

%% 

filename = 'sensor_hidden_dynamics_times.mat';
save(filename,'timesave1','timesave2','timesave3','timesave4')
