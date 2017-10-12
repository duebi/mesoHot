function plot_REG_phase_h(depths, xdata, ydata, rho, pval,h)

fig = figure('Units','inches');
fig.Position(3:4) = [6.5 4.5];
for j = 1:length(depths)
    
    ns = length(xdata{j});
    
    subplot(2,4,j)
    plot(xdata{j},ydata{j},'o');
    hold on;
    plot(xdata{j},xdata{j}*rho(j),'r','LineWidth',1.5);
    hold off;
    title_str = sprintf('%dm\n(p=%.1g, n=%d)',depths(j),pval(j),ns);
    if h(j)
        title(title_str,'FontWeight','bold','FontSize',12);
        set(gca,'LineWidth',2);
    else
        title(title_str,'FontWeight','normal','FontSize',12);
    end
    axis square;
    ax = gca;
    cmax = max(abs([ax.XLim ax.YLim]));
    xlim([-cmax cmax]);
    ylim([-cmax cmax]);
    if j<=4; ax.Position(2) = ax.Position(2)-.05; end
end

end