function save_conf_ellipses(n,s)

fig = figure(n);

hL = legend('Exact DMD','fbDMD','tlsDMD','optimized DMD',...
    'Exact DMD','fbDMD','tlsDMD','optimized DMD','Answer');

axis equal

newPosition = [0.75 0.25 0.1 0.1];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
set(hL,'Box','on')


set(fig,'Position',[100 100 1049 895])

print(fig,s,'-depsc')
