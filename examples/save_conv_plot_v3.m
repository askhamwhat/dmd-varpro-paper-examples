
function save_conv_plot_v3(n,s)

fig = figure(n);

hL = legend('Exact DMD','tlsDMD','optimized DMD','POD');

newPosition = [0.65 0.65 0.2 0.2];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
set(hL,'Box','on')

set(fig,'Position',[100 100 1049 895])

print(fig,s,'-depsc')
