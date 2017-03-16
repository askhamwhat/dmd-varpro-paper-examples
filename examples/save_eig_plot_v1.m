
function save_eig_plot_v1(n,s)

fig = figure(n);

hL = legend('Exact DMD','tlsDMD','optimized DMD');

newPosition = [0.65 0.15 0.2 0.2];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
set(hL,'Box','on')

set(fig,'Position',[100 100 1049 895])

print(fig,s,'-depsc')
