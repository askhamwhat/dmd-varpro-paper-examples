
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

figure(1)

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

figure(2)

for j = 1:4
    subplot(2,2,j)
    loglog(ms,erreig1hid(:,j),'-.kd')
    hold on
    loglog(ms,erreig2hid(:,j),'--b+')
    loglog(ms,erreig3hid(:,j),'--rs')
    loglog(ms,erreig4hid(:,j),'-.gx')
end

% plot confidence ellipses for highest noise, 
% fewest snapshots case for the first 
% eigenvalue of dominant type

figure(3)

ieig = idom(1);
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


% plot confidence ellipses for highest noise, 
% fewest snapshots case for the first 
% eigenvalue of hidden type

figure(4)

ieig = ihid(1);
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

figure(5)

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


