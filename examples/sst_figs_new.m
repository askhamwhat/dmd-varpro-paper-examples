
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SST figure generation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load precomputed data

close all
clear all

prettyfigs_dvp

filename = 'sst_even_new.mat';
load(filename)
filename = 'sst_uneven_new.mat';
load(filename)

load('sst_mask.mat')

%% singular values plot

figure(1)

clf

r999 = sum(cumsum(sd) < .999*sum(sd))+1;
r99 = sum(cumsum(sd) < .99*sum(sd))+1;
r90 = sum(cumsum(sd) < .9*sum(sd))+1;

semilogy(sd,'b-')
hold on
semilogy(rgd,sd(rgd),'ms')
semilogy(r90,sd(r90),'rx')
semilogy(r99,sd(r99),'gd')
%semilogy(r999,sd(r999),'cx')

save_sing_plot(1,'sst_sing_vals');

%% plot eigenvalues obtained for gavish donoho cut-off

figure(2)

ind1 = 7;
ind2 = 5;
ind3 = 3;

alleigs = [esexact(1:rs(ind1),ind1); 
    estls(1:rs(ind1),ind1); 
    esopt(1:rs(ind1),ind1);
    esexact(1:rs(ind2),ind2); 
    estls(1:rs(ind2),ind2); 
    esopt(1:rs(ind2),ind2);
    esexact(1:rs(ind3),ind3); 
    estls(1:rs(ind3),ind3); 
    esopt(1:rs(ind3),ind3)];

xmin = min(real(alleigs(:)));
xmax = max(real(alleigs(:)));
ymin = min(imag(alleigs(:)));
ymax = max(imag(alleigs(:)));

subplot(2,2,1)

eexact = esexact(1:rs(ind1),ind1);
etls = estls(1:rs(ind1),ind1);
eopt = esopt(1:rs(ind1),ind1);

scatter(real(eexact),imag(eexact),'rx')
hold on
scatter(real(etls),imag(etls),'bd')
scatter(real(eopt),imag(eopt),'go')
s = sprintf('$r = %i$',rs(ind1));
title(s,'interpreter','LaTeX')

xlim( [xmin xmax] )

subplot(2,2,2)

eexact = esexact(1:rs(ind2),ind2);
etls = estls(1:rs(ind2),ind2);
eopt = esopt(1:rs(ind2),ind2);

scatter(real(eexact),imag(eexact),'rx')
hold on
scatter(real(etls),imag(etls),'bd')
scatter(real(eopt),imag(eopt),'go')
s = sprintf('$r = %i$',rs(ind2));
title(s,'interpreter','LaTeX')

xlim( [xmin xmax] )

subplot(2,2,3)

eexact = esexact(1:rs(ind3),ind3);
etls = estls(1:rs(ind3),ind3);
eopt = esopt(1:rs(ind3),ind3);

scatter(real(eexact),imag(eexact),'rx')
hold on
scatter(real(etls),imag(etls),'bd')
scatter(real(eopt),imag(eopt),'go')
s = sprintf('$r = %i$',rs(ind3));
title(s,'interpreter','LaTeX')

xlim( [xmin xmax] )


% save figure

save_eig_plot_v1(2,'sst_eig_vs_rank');

%% create table --- all methods on evenly spaced

[bsort,isort] = sort(abs(boptsave),'descend');

eoptsort = eoptsave(isort);

ind = 4;

eexact = esexact(1:rs(ind),ind);
etls = estls(1:rs(ind),ind);

indexact = match_vectors(imag(eexact),imag(eoptsort));
indtls = match_vectors(imag(etls),imag(eoptsort));

% call latexTable to output a latex table

v1 = eexact(indexact);
v2 = etls(indtls);
v3 = eoptsort;

clear input

input.data = [boptsave(isort), 2*pi./imag(v3), 2*pi./imag(v2), 2*pi./imag(v1)];
input.dataFormat = {'%+6.4e'};
input.tableColLabels = {'Coefficient','optimized DMD','tlsDMD','exact DMD'};
latex_all3 = latexTable(input);

%% create table --- compare wavelengths for even and uneven

ind_u = match_vectors(eoptsave_u,eoptsort);

v1 = eoptsort;
v2 = eoptsave_u(ind_u);

projs = zeros(size(eoptsort));

for i = 1:length(eoptsort)
    w1 = woptsave(:,isort(i));
    w2 = woptsave_u(:,ind_u(i));
    projs(i) = abs( w1'*w2 )/( norm(w1)*norm(w2) );
end

clear input

input.data = [boptsave(isort), boptsave_u(ind_u), 2*pi./imag(v1), 2*pi./imag(v2), projs];
input.dataFormat = {'%+6.4e'};
input.tableColLabels = {'even --- $b$','uneven --- $b$','even --- $\lambda$','uneven --- $\lambda$','projection'};
latex_evsu = latexTable(input);

%% create table --- compare timings

clear input

input.data = [rs(:), timesexact(:), timestls(:), timesopt(:)];

input.dataFormat = {'%i',1,'%6.4e',3};
input.tableColLabels = {'rank','$t_{DMD}$','$t_{tlsDMD}$','$t_{optDMD}$'};

latex_times = latexTable(input);

%% create table --- compare reconstruction quality

fig = figure(4)
clf 
rspod = zeros(size(rs));

for i = 1:length(rs)
    rspod(i) = sqrt(sum(sd(rs(i)+1:end).^2))/sqrt(sum(sd.^2));
end

rsall = [rerrsexact(:); rerrstls(1:end-1); rerrsopt(:); rspod(:)];
ymin = min(rsall);
ymax = max(rsall);

loglog(rs,rerrsexact(:),'-.kd')
hold on
loglog(rs(1:end-1),rerrstls(1:end-1),'--rs')
loglog(rs,rerrsopt(:),'-.gx')
loglog(rs,rspod(:),'-.mo')

ylim([ymin ymax])
xlim( [ min(rs) max(rs) ] )

save_conv_plot_v3(4,'sst_rec_plot')

clear input


input.data = [rs(:), rerrsexact(:), rerrstls(:), rerrsopt(:), rspod(:)];

input.dataFormat = {'%i',1,'%6.4e',4};
input.tableColLabels = {'rank','$\rho_{DMD}$','$\rho_{tlsDMD}$', ... 
    '$\rho_{optDMD}$','$\rho_{POD}$'};

latex_times = latexTable(input);


%% figure --- compare modes for even and uneven

half_year = [1 2];
year = [3 4];
static = [13 14];

ind_u = match_vectors(imag(eoptsave_u),imag(eoptsave));

fig = figure(5)

subplot(3,4,1)
plot_with_mask(real(woptsave_u(:,ind_u(static(2)))),sst_mask)
set(gca,'YTick',[],'XTick',[])
subplot(3,4,3)
plot_with_mask(real(woptsave(:,static(2))),sst_mask)
set(gca,'YTick',[],'XTick',[])
subplot(3,4,5)
plot_with_mask(real(woptsave_u(:,ind_u(year(1)))),sst_mask)
set(gca,'YTick',[],'XTick',[])
subplot(3,4,6)
plot_with_mask(imag(woptsave_u(:,ind_u(year(1)))),sst_mask)
set(gca,'YTick',[],'XTick',[])
subplot(3,4,7)
plot_with_mask(real(woptsave(:,year(1))),sst_mask)
set(gca,'YTick',[],'XTick',[])
subplot(3,4,8)
plot_with_mask(imag(woptsave(:,year(1))),sst_mask)
set(gca,'YTick',[],'XTick',[])
subplot(3,4,9)
plot_with_mask(real(woptsave_u(:,ind_u(half_year(2)))),sst_mask)
set(gca,'YTick',[],'XTick',[])
subplot(3,4,10)
plot_with_mask(imag(woptsave_u(:,ind_u(half_year(2)))),sst_mask)
set(gca,'YTick',[],'XTick',[])
subplot(3,4,11)
plot_with_mask(real(woptsave(:,half_year(1))),sst_mask)
set(gca,'YTick',[],'XTick',[])
subplot(3,4,12)
plot_with_mask(-imag(woptsave(:,half_year(1))),sst_mask)
set(gca,'YTick',[],'XTick',[])

print(fig,'sst_modes_evsu','-dpng')
