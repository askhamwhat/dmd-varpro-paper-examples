
function save_conv_plot_v2(n,s)

fig = figure(n);

hL = legend('Exact DMD','fbDMD','tlsDMD','optimized DMD');

newPosition = [0.65 0.65 0.2 0.2];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
set(hL,'Box','on')

set(fig,'Position',[100 100 1049 895])

print(fig,s,'-depsc')
