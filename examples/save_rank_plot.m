
function save_rank_plot(n,s)

fig = figure(n);

hL = legend('Gavish-Donoho','99.9 percent','99 percent','90 percent');

posVec = [0 0 900 600];
set(fig,'Units','points','Position',posVec,'PaperPositionMode','auto');

for i = 1:5
    subplot(2,3,i)
    set(gca,'FontSize',16)
end

newPosition = [0.75 0.25 0.1 0.1];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
set(hL,'Box','on')


print(fig,s,'-depsc')
