
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate figures for 
% hidden dynamics system example with
% additive sensor noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;

prettyfigs_dvp

% load data

% eigenvalues saved as (:,ntrials,nsigs,nm)
% recovered A saved as (:,:,ntrials,nsigs,nm)
% method 1 - dmd
% method 2 - fbdmd
% method 3 - tlsdmd
% method 4 - optdmd

filename = 'sensor_hidden_dynamics.mat';
load(filename);

devals = evals; % evals 1 and 2 are dominant dynamics
                % evals 3 and 4 are "hidden"
sigs = sigs(:,1);
[~,ntrials,nsigs,nm] = size(eigsave1);

%[u,s,v] = svd(Xclean,'econ');

%% errors in eigenvalues (dominant)

erreig1dom = zeros(nm,nsigs);
erreig2dom = zeros(nm,nsigs);
erreig3dom = zeros(nm,nsigs);
erreig4dom = zeros(nm,nsigs);

idom = [1;2];

for j = 1:nm
    for i = 1:nsigs
        err1 = 0;
        err2 = 0;
        err3 = 0;
        err4 = 0;
        for k = 1:ntrials
            err1 = err1 + norm(devals(idom)-eigsave1(idom,k,i,j))/norm(devals(idom));
            err2 = err2 + norm(devals(idom)-eigsave2(idom,k,i,j))/norm(devals(idom));
            err3 = err3 + norm(devals(idom)-eigsave3(idom,k,i,j))/norm(devals(idom));
            err4 = err4 + norm(devals(idom)-eigsave4(idom,k,i,j))/norm(devals(idom));            
        end
        erreig1dom(j,i) = err1/ntrials;
        erreig2dom(j,i) = err2/ntrials;
        erreig3dom(j,i) = err3/ntrials;
        erreig4dom(j,i) = err4/ntrials;
    end
end

%%
% plot errors

figure(2)

for j = 1:4
    subplot(2,2,j)
    loglog(ms,erreig1dom(:,j),'-.kd')
    hold on
    loglog(ms,erreig2dom(:,j),'--b+')
    loglog(ms,erreig3dom(:,j),'--rs')
    loglog(ms,erreig4dom(:,j),'-.gx')
end

%% errors in eigenvalues (HIDDEN)

erreig1hid = zeros(nm,nsigs);
erreig2hid = zeros(nm,nsigs);
erreig3hid = zeros(nm,nsigs);
erreig4hid = zeros(nm,nsigs);

ihid = [3;4];

for j = 1:nm
    for i = 1:nsigs
        err1 = 0;
        err2 = 0;
        err3 = 0;
        err4 = 0;
        for k = 1:ntrials
            err1 = err1 + norm(devals(ihid)-eigsave1(ihid,k,i,j))/norm(devals(ihid));
            err2 = err2 + norm(devals(ihid)-eigsave2(ihid,k,i,j))/norm(devals(ihid));
            err3 = err3 + norm(devals(ihid)-eigsave3(ihid,k,i,j))/norm(devals(ihid));
            err4 = err4 + norm(devals(ihid)-eigsave4(ihid,k,i,j))/norm(devals(ihid));            
        end
        erreig1hid(j,i) = err1/ntrials;
        erreig2hid(j,i) = err2/ntrials;
        erreig3hid(j,i) = err3/ntrials;
        erreig4hid(j,i) = err4/ntrials;
    end
end

%%
% plot errors

figure(3)

for j = 1:4
    subplot(2,2,j)
    loglog(ms,erreig1hid(:,j),'-.kd')
    hold on
    loglog(ms,erreig2hid(:,j),'--b+')
    loglog(ms,erreig3hid(:,j),'--rs')
    loglog(ms,erreig4hid(:,j),'-.gx')
end

