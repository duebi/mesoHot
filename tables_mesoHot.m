% script tables_mesoHot.m
%
% Calls the function mesoHot.m for selected biogeochemical variables and
% computes regression significance using the Benjiamini and Hochberg
% procedure.
% Results are collected in two tables that are printed as a result of
% calling this script.
%
% Benedetto Barone - Oct 2017

%% Tables with results from biogeochemical variables correlations with SLA
varnames = {'nit','phos','sil','dic','boxy','chl','pc','pn','pp','doc','don','dop','hbact','pro','euk','l12'};
%varnames = {'chlb','fuco','hex','but','peri','pras','diad','zeax','chlda','viol','mvchl','dvchl','mvdv'};
% Initialize variables
lv = length(varnames);

p_lin = NaN(8,lv); 
p_AR = NaN(8,lv);
p_rks = NaN(8,lv);
r_lin = NaN(8,lv); 
r_AR = NaN(8,lv);
% Extract correlation coefficients and p values
for i = 1:length(varnames)
    [rs,pvals,rs_AR,pvals_AR,pvals_rks] =  mesoHot(varnames{i},'Spearman');
    p_lin(:,i) = pvals(:,2);
    r_lin(:,i) = rs(:,2);
    p_AR(:,i) = pvals_AR;
    r_AR(:,i) = rs_AR;
    p_rks(:,i) = pvals_rks;
    clear rs pvals rs_AR pvals_AR pvals_rks
    close
end
% Significance from BENJAMINI & HOCHBERG PROCEDURE
[h_lin, ~, ~, ~]=fdr_bh(p_lin,0.05);
[h_AR, ~, ~, ~]=fdr_bh(p_AR,0.05);
[h_rks, ~, ~, ~]=fdr_bh(p_rks,0.05);
% Plot tables
f1 = plot_REG_table_h([5 25 45 75 100 125 150 175],varnames,r_lin,h_lin);
%print('./results/H&B_table_lin','-dpng');
f2 = plot_REG_table_h([5 25 45 75 100 125 150 175],varnames,r_AR,h_AR);
%print('./results/H&B_table_AR','-dpng');

f3 = plot_PVAL_table([5 25 45 75 100 125 150 175],varnames,p_lin);

%% Tables with results from pigment correlations with SLA
varnames = {'hplc','chlb','chlc','fuco','hex','but','peri','diad','zeax','chlda','viol','mvchl','dvchl','mvdv'};
% Initialize variables
p_lin = NaN(8,14); 
p_AR = NaN(8,14);
p_rks = NaN(8,14);
r_lin = NaN(8,14); 
r_AR = NaN(8,14);
% Extract correlation coefficients and p values
for i = 1:length(varnames)
    [rs,pvals,rs_AR,pvals_AR,pvals_rks] =  mesoHot(varnames{i},'Spearman');
    p_lin(:,i) = pvals(:,2);
    r_lin(:,i) = rs(:,2);
    p_AR(:,i) = pvals_AR;
    r_AR(:,i) = rs_AR;
    p_rks(:,i) = pvals_rks;
    clear rs pvals rs_AR pvals_AR pvals_rks
    close
end
% Significance from BENJAMINI & HOCHBERG PROCEDURE
[h_lin, ~, ~, ~]=fdr_bh(p_lin,0.05);
[h_AR, ~, ~, ~]=fdr_bh(p_AR,0.05);
[h_rks, ~, ~, ~]=fdr_bh(p_rks,0.05);
% Plot tables
f1 = plot_REG_table_h([5 25 45 75 100 125 150 175],varnames,r_lin,h_lin);
%print('./results/H&B_table_lin','-dpng');
f2 = plot_REG_table_h([5 25 45 75 100 125 150 175],varnames,r_AR,h_AR);
%print('./results/H&B_table_AR','-dpng');

%% Monthly analysis on surface chlorophyll
[rs,pvals,rs_AR,pvals_AR,pvals_rks,datad] =  mesoHot('chl','Spearman');
for i = 1:12
    pnew = datad{1}(datad{1}.c_sla>=5,:);
    A{i} = pnew.chl(pnew.month==i);
    pnew = datad{1}(datad{1}.c_sla<=-5,:);
    B{i} = pnew.chl(pnew.month==i);
    pp(1:3,i) = prctile(A{i},[2.5 50 97.5]);
    pm(1:3,i) = prctile(B{i},[2.5 50 97.5]);
    
    pmonth(i) = ranksum(A{i},B{i});
end
mnth = 1:12;

clf
patch([1:12 12:-1:1],[pp(1,:) pp(3,12:-1:1)],[0.7 0.5 0.5],'Edgecolor','none','Facealpha',0.4)
patch([1:12 12:-1:1],[pm(1,:) pm(3,12:-1:1)],[0.5 0.5 0.7],'Edgecolor','none','Facealpha',0.4)
hold on, plot(1:12,pp(2,:),'r',1:12,pm(2,:),'b'), hold off
hold on, plot(mnth(pmonth<0.05),0*mnth(pmonth<0.05),'ko','MarkerFaceColor','k','Markersize',12), hold off
xlim([1 12])
set(gca,'box','on')
legend('>5cm 95%','<-5cm 95%','median >5cm','median <5cm')
xlabel('month')
ylabel('chl a (mg m-3)')

%% Single value per HOT cruise - tables with results from biogeochemical variables correlations with SLA
% Biogeochemical variables to be extracted
varnames = {'nit','phos','sil','dic','boxy','chl','pc','pn','pp','doc','don','dop','hbact','pro','euk','l12'};
% Initialize variables
p_cr = NaN(8,16);
rho_cr = NaN(8,16);
for i = 1:length(varnames)
    [rs,pvals,rs_AR,pvals_AR,pvals_rks,datad] =  mesoHot(varnames{i},'Spearman');
    for j = 1:8
        pnew = grpstats(datad{j},'cruise');
        [rho_temp,p_temp] = corr(pnew.mean_c_sla,eval(['pnew.mean_' datad{1}.Properties.VariableNames{7}]),'type','Spearman');
        rho_cr(j,i) = rho_temp; 
        p_cr(j,i) = p_temp;
    end    
end

[h_cr, ~, ~, ~]=fdr_bh(p_cr,0.05);
% Plot tables
f1 = plot_REG_table_h([5 25 45 75 100 125 150 175],varnames,rho_cr,h_cr);