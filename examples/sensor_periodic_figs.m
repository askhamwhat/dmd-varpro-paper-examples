
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
% method 4 - optdmd

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

% plot confidence ellipses for highest noise, 
% fewest snapshots case for the first 
% eigenvalue

figure(3)

ieig = 2;
jsig = 2;
km = 1;

zdat1 = transpose(eigsave1(ieig,:,jsig,km));
zdat2 = transpose(eigsave2(ieig,:,jsig,km));
zdat3 = transpose(eigsave3(ieig,:,jsig,km));
zdat4 = transpose(eigsave4(ieig,:,jsig,km));

[c1,ax1,lens1] = zconf95(zdat1);
[c2,ax2,lens2] = zconf95(zdat2);
[c3,ax3,lens3] = zconf95(zdat3);
[c4,ax4,lens4] = zconf95(zdat4);

ts = linspace(0,2*pi,1000);
xys1 = c1' + ax1(:,1)*lens1(1)*cos(ts)+ax1(:,2)*lens1(2)*sin(ts);
xys2 = c2' + ax2(:,1)*lens2(1)*cos(ts)+ax2(:,2)*lens2(2)*sin(ts);
xys3 = c3' + ax3(:,1)*lens3(1)*cos(ts)+ax3(:,2)*lens3(2)*sin(ts);
xys4 = c4' + ax4(:,1)*lens4(1)*cos(ts)+ax4(:,2)*lens4(2)*sin(ts);

plot(xys1(1,:),xys1(2,:),'-.k')
hold on
%scatter(real(zdat1),imag(zdat1),'k')
%scatter(real(zdat2),imag(zdat2),'b')
%scatter(real(zdat3),imag(zdat3),'r')
%scatter(real(zdat4),imag(zdat4),'g')
plot(xys2(1,:),xys2(2,:),'-b')
plot(xys3(1,:),xys3(2,:),'--r')
plot(xys4(1,:),xys4(2,:),'-.g')
scatter(c1(1),c1(2),'kd')
scatter(c2(1),c2(2),'b+')
scatter(c3(1),c3(2),'rs')
scatter(c4(1),c4(2),'gx')
scatter(real(devals(ieig)),imag(devals(ieig)),'mo')

% reconstruction errors

figure(4)

temp = mean(rerrsave1);
rerr1 = transpose(reshape(temp,nsigs,nm));
temp = mean(rerrsave2);
rerr2 = transpose(reshape(temp,nsigs,nm));
temp = mean(rerrsave3);
rerr3 = transpose(reshape(temp,nsigs,nm));
temp = mean(rerrsave4);
rerr4 = transpose(reshape(temp,nsigs,nm));


for j = 1:4
    subplot(2,2,j) 
    loglog(ms,rerr1(:,j),'-.kd')
    hold on
    loglog(ms,rerr2(:,j),'--b+')
    loglog(ms,rerr3(:,j),'--rs')
    loglog(ms,rerr4(:,j),'-.gx')
end
