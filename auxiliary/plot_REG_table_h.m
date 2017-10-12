function fig = plot_REG_table_h(depths, varnames, rho, h)

for i = 1:length(depths)
    for j = 1:length(varnames)
        if h(i,j)
            rho_cell{i,j} = sprintf('%+.2f',rho(i,j));
        else
            rho_cell{i,j} = '  NS';
        end
    end
end
depth_str = strcat(cellstr(num2str(depths')),' m');
fig = figure('Units','inches');
fig.Position(3:4) = [6.5 4.1];
ut = uitable(fig,'Data',rho_cell','RowName',varnames,'ColumnName',depth_str,...
    'ColumnWidth',{50},'FontSize',10,'Units','normalized','Position',[0 0 1 1]);
%T = array2table(data,'VariableNames',varnames,'RowNames',depth_str);

end