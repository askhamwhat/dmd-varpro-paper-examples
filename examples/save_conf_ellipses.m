function save_conf_ellipses(n,s)

fig = figure(n);

hL = legend('Exact DMD','fbDMD','tlsDMD','optimized DMD',...
    'Exact DMD mean','fbDMD mean','tlsDMD mean','optimized DMD mean', ...
    'Answer');

axis equal

posVec = [0 0 900 600];
set(fig,'Units','points','Position',posVec,'PaperPositionMode','auto');

newPosition = [0.75 0.25 0.1 0.1];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
set(hL,'Box','on')

set(gca,'TickLength',[.01,.025])
set(gca,'FontSize',16)

print(fig,s,'-depsc')
