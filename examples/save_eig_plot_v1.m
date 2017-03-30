
function save_eig_plot_v1(n,s)

fig = figure(n);

hL = legend('Exact DMD','tlsDMD','optimized DMD');

posVec = [0 0 900 600];
set(fig,'Units','points','Position',posVec,'PaperPositionMode','auto');


for i = 1:3
    subplot(2,2,i)
    set(gca,'FontSize',16)
end


newPosition = [0.65 0.15 0.2 0.2];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
set(hL,'Box','on')

print(fig,s,'-depsc')
