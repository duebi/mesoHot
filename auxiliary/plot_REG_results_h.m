function plot_REG_results_h(depths, varnames, rho,h)

fig = figure('Units','inches');
fig.Position(3:4) = [6.5 3];
axes('Position',[.06 .05 .93 .7]);
imagesc(rho,[-.6 .6]);
hold on;
[tmpi,tmpj] = find(h);
plot(tmpj,tmpi,'ks','MarkerSize',10);
hold off;
colorbar();
set(gca,'XAxisLocation','top','XTick',1:length(varnames),'XTickLabel',varnames,'XTickLabelRotation',90);
set(gca,'YTick',1:length(depths),'YTickLabel',depths);
ylabel('depth (m)');
colormap(zerocmap(redbluecmap))

end