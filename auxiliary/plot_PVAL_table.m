function fig = plot_PVAL_table(depths, varnames, pvals)

for i = 1:length(depths)
    for j = 1:length(varnames)
           pval_cell{j,i} = sprintf('%1.1e',pvals(i,j));
    end
end
depth_str = strcat(cellstr(num2str(depths')),' m');
fig = figure('Units','inches');
fig.Position(3:4) = [18 3];
ut = uitable(fig,'Data',pval_cell','RowName',depth_str,'ColumnName',varnames,...
    'ColumnWidth',{70},'FontSize',16,'Units','normalized','Position',[0 0 1 1]);
%T = array2table(data,'VariableNames',varnames,'RowNames',depth_str);

end