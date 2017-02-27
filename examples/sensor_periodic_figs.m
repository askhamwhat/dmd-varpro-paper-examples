
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate figures for 
% periodic system example with
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
% method 4 - dmdef

filename = 'sensor_periodic.mat';
load(filename);

devals = diag(evals);
[~,ntrials,nsigs,nm] = size(eigsave1);

[u,s,v] = svd(Xclean,'econ');

%% errors in reconstructed A

erramat1 = zeros(nm,nsigs);
erramat2 = zeros(nm,nsigs);
erramat3 = zeros(nm,nsigs);
erramat4 = zeros(nm,nsigs);

for i = 1:nsigs
    for j = 1:nm
        err1 = 0;
        err2 = 0;
        err3 = 0;
        err4 = 0;
        for k = 1:ntrials
            err1 = err1 + norm(Admd-asave1(:,:,k,i,j),'fro')/norm(Admd,'fro');
            err2 = err2 + norm(Admd-asave2(:,:,k,i,j),'fro')/norm(Admd,'fro');
            err3 = err3 + norm(Admd-asave3(:,:,k,i,j),'fro')/norm(Admd,'fro');
            err4 = err4 + norm(Admd-expm(asave4(:,:,k,i,j)*dt),'fro')/norm(Admd,'fro');
        end
        erramat1(j,i) = err1/ntrials;
        erramat2(j,i) = err2/ntrials;
        erramat3(j,i) = err3/ntrials;
        erramat4(j,i) = err4/ntrials;
    end
end

% plot errors

figure(1)

for j = 1:4
    subplot(2,2,j) 
    loglog(ms,erramat1(:,j),'-.kd')
    hold on
    loglog(ms,erramat2(:,j),'--b+')
    loglog(ms,erramat3(:,j),'--rs')
    loglog(ms,erramat4(:,j),'-.gx')
end

%% errors in eigenvalues

erreig1 = zeros(nm,nsigs);
erreig2 = zeros(nm,nsigs);
erreig3 = zeros(nm,nsigs);
erreig4 = zeros(nm,nsigs);

for i = 1:nsigs
    for j = 1:nm
        err1 = 0;
        err2 = 0;
        err3 = 0;
        err4 = 0;
        for k = 1:ntrials
            err1 = err1 + norm(devals-eigsave1(:,k,i,j))/norm(devals);
            err2 = err2 + norm(devals-eigsave2(:,k,i,j))/norm(devals);
            err3 = err3 + norm(devals-eigsave3(:,k,i,j))/norm(devals);
            err4 = err4 + norm(devals-eigsave4(:,k,i,j))/norm(devals);            
        end
        erreig1(j,i) = err1/ntrials;
        erreig2(j,i) = err2/ntrials;
        erreig3(j,i) = err3/ntrials;
        erreig4(j,i) = err4/ntrials;
    end
end

% plot errors

figure(2)

for j = 1:4
    subplot(2,2,j) 
    loglog(ms,erreig1(:,j),'-.kd')
    hold on
    loglog(ms,erreig2(:,j),'--b+')
    loglog(ms,erreig3(:,j),'--rs')
    loglog(ms,erreig4(:,j),'-.gx')
end

