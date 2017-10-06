function varnames = variables_in_folder(dir_str)
% returns a cell with variable names in the directory specified by dir_str

file_str = dir(dir_str);
varnames = {};
for j = 1:length(file_str)
    tmplength = strfind(file_str(j).name,'_0_200m.txt');
    if ~isempty(tmplength)
        varnames{end+1} = file_str(j).name(1:tmplength-1);
    end
end

end