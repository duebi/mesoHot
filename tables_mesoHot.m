% script tables_mesoHot.m
%
% Calls the function mesoHot.m for selected biogeochemical variables and
% computes regression significance using the Benjiamini and Hochberg
% procedure.
% Results are collected in two tables that are printed as a result of
% calling this script.
%
% Benedetto Barone - Oct 2017

% Biogeochemical variables to be extracted
varnames = {'chl','l12','dic','pc','doc','nit','pn','don','phos','pp','dop','pro','hbact','euk','ph','sil'};
% Initialize variables
p_lin = NaN(8,16); 
p_AR = NaN(8,16);
r_lin = NaN(8,16); 
r_AR = NaN(8,16);
% Extract correlation coefficients and p values
for i = 1:length(varnames)
    [rs,pvals,rs_AR,pvals_AR] =  mesoHot(varnames{i});
    p_lin(:,i) = pvals(:,2);
    r_lin(:,i) = rs(:,2);
    p_AR(:,i) = pvals_AR;
    r_AR(:,i) = rs_AR;
    clear rs pvals rs_AR pval_AR
    close
end
% Significance from BENJAMINI & HOCHBERG PROCEDURE
[h_lin, ~, ~, ~]=fdr_bh(p_lin,0.05);
[h_AR, ~, ~, ~]=fdr_bh(p_AR,0.05);
% Plot tables
f1 = plot_REG_table_h([5 25 45 75 100 125 150 175],varnames,r_lin,h_lin);
print('./results/H&B_table_lin','-dpng');
f2 = plot_REG_table_h([5 25 45 75 100 125 150 175],varnames,r_AR,h_AR);
print('./results/H&B_table_AR','-dpng');