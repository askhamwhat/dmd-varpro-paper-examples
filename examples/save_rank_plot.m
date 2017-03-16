
function save_rank_plot(n,s)

fig = figure(n);

hL = legend('Gavish-Donoho','99.9 percent','99 percent','90 percent');

newPosition = [0.75 0.25 0.1 0.1];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
set(hL,'Box','on')

set(fig,'Position',[100 100 1049 895])

print(fig,s,'-depsc')
