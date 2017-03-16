function save_sing_plot(n,s)

fig = figure(n);

hL = legend('Singular values','Gavish-Donoho','90 percent','99 percent');

newPosition = [0.65 0.65 0.2 0.2];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
set(hL,'Box','on')


set(fig,'Position',[100 100 1049 895])

print(fig,s,'-depsc')
