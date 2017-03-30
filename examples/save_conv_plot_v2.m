
function save_conv_plot_v2(n,s)

fig = figure(n);

hL = legend('Exact DMD','fbDMD','tlsDMD','optimized DMD');

posVec = [0 0 900 600];
set(fig,'Units','points','Position',posVec,'PaperPositionMode','auto');

set(gca,'FontSize',16)

newPosition = [0.65 0.65 0.2 0.2];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
set(hL,'Box','on')

set(gca,'TickLength',[.01,.025])

print(fig,s,'-depsc')
