function subplot_labels(title_str, x_str, y_str)

gcf;
ax = axes('Position',[0 0 1 1]);
ax.Visible = 'off';
text(.5,1,title_str,'VerticalAlignment','top','HorizontalAlignment','center','FontSize',14,'FontWeight','bold');
text(.5,0,x_str,'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',14);
text(0,.5,y_str,'Rotation',90,'VerticalAlignment','top','HorizontalAlignment','center','FontSize',14);

end
