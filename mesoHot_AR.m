% Adapted from mesoHot.m.

% Regression between biogeochemical variable autoregression residuals and
% SLA (autoregression residuals, case 1; normalized timeseries, case 2).

% Dependencies: plot_REG_phase.m, plot_REG_results.m, plot_REG_table.m,
% redbluecmap.m, subplot_labels.m, zerocmap.m, as well as the dependencies
% of mesoHot.m (readHotDogs.m, data directory, gsw and lsq scripts).

% Ashley Coenen, 2017


clear;
close all;

% grab all the variables from the data folder
%varnames = variables_in_folder('./data');
varnames = {'chl','l12','dic','pc','doc','nit','pn','don','phos','pp','dop','pro','hbact','euk','ph','sil'};

load ./data/allHawaii
% loop over all the variables
for k = 1:length(varnames)
      
% Extract data for selected variable
data = readHotDogs(['./data/' varnames{k} '_0_200m.txt'],1);
% Compute pigment ratio in the case of mvdv
if length(data.Properties.VariableNames)>7 && strcmp(data.Properties.VariableNames(7),'mvchla') && strcmp(data.Properties.VariableNames(8),'dvchla')
   data.mvchla = data.mvchla./data.dvchla;
   data(:,8) = [];
   data.Properties.VariableNames(7) = {'mvdv'};
   data.Properties.VariableUnits(7) = {'g/g'};
end
% Identify variable name and units from data table
varname = data.Properties.VariableNames{end};
varunits = data.Properties.VariableUnits{end};
% Remove rows with missing values
data(isnan(data.date) | eval(['isnan(data.' varname ') | isinf(data.' varname ')']),:) = [];
% Interpolate SLA on bottle data dates
data.r_sla = interp1(allHawaii.date -10/24,allHawaii.aloha,data.date); % uncorrected SLA
data.c_sla = interp1(allHawaii.date -10/24,allHawaii.deseas,data.date); % de-trended, no season, SLA

depths = [5 25 45 75 100 125 150 175];

% Extract data at each depth, compute regression, compute statistics for
    for j = 1:length(depths)
        ind_dpt = data.depth >= depths(j)-3 & data.depth <= depths(j)+3;
        % increase depth layer width for 100m for tdn, tdp, don, and dop (Explained in the manuscript)
        if (strcmp(varname,'tdn') || strcmp(varname,'tdp') || strcmp(varname,'don') || strcmp(varname,'dop')) && depths(j)==100
            ind_dpt = data.depth >= depths(j)-7 & data.depth <= depths(j)+7;
        end
        % Isolate data from the depth layer & count data points
        ns(j) = sum(ind_dpt);
        ddata = data(ind_dpt,:);

        
        % FOR SAMPLING TIME HISTORGRAMS
        dt = ddata.date(2:end)-ddata.date(1:end-1);
        dt_bins = 0:15:120;
        [dt_counts, ~] = hist(dt,dt_bins);
        %dt_frac(:,j,k) = dt_counts/sum(dt_counts);
        dt_frac(:,j,k) = dt_counts;
        
        
        %%%%%%%%%%%%%%%%%%%%%
        % BEGIN AR ANALYSIS %
        %%%%%%%%%%%%%%%%%%%%%

        % autoregression of variable k at depth j
        x1 = ddata.(varname)(1:end-1); % 'x'
        x2 = ddata.(varname)(2:end);   % 'y'

        % normalize
        z1 = (x1-mean(x1))/std(x1);
        z2 = (x2-mean(x2))/std(x2);

        % autoregression, record residuals
        rho = sum(z2.*z1)/sum(z1.^2);
        eps_var{j,k} = z2 - rho*z1;
        
        clear x1 x2 z1 z2 rho;
        

        % autoregression of SLA for variable k at depth j
        x1 = ddata.c_sla(1:end-1);
        x2 = ddata.c_sla(2:end);

        % normalize
        z1 = (x1-mean(x1))/std(x1);
        z2 = (x2-mean(x2))/std(x2);

        % autoregression, record residuals
        rho = sum(z2.*z1)/sum(z1.^2);
        eps_sla{j,k} = z2 - rho*z1;

        % record normalized SLA
        z1_sla{j,k} = z1;
        
        clear x1 x2 z1 z2 rho;
        


        % correlation between variable residual and SLA residual
        [rho_eps(j,k),pval_eps(j,k)] = corr(eps_var{j,k},eps_sla{j,k});

        % correlation between variable residual and normalized SLA
        [rho_sla(j,k),pval_sla(j,k)] = corr(eps_var{j,k},z1_sla{j,k});


    end

end

% Significance from BENJAMINI & HOCHBERG PROCEDURE
[h_eps, ~, ~, ~]=fdr_bh(pval_eps,0.05);
[h_sla, ~, ~, ~]=fdr_bh(pval_sla,0.05);

for k = 1:length(varnames)
    
    plot_REG_phase_h(depths,eps_var(:,k),eps_sla(:,k),rho_eps(:,k),pval_eps(:,k),h_eps(:,k));
    subplot_labels(sprintf('%s vs SLA (autoregressed)',varnames{k}),'SLA residual',sprintf('%s residual',varnames{k}));
    print(sprintf('./results/SLAeps_%s',varnames{k}),'-dpng');
    close;

    plot_REG_phase_h(depths,eps_var(:,k),z1_sla(:,k),rho_sla(:,k),pval_sla(:,k),h_sla(:,k));
    subplot_labels(sprintf('%s vs SLA (normalized)',varnames{k}),'normalized SLA',sprintf('%s residual',varnames{k}));
    print(sprintf('./results/SLAnorm_%s',varnames{k}),'-dpng');
    close;
    
end

%% SUMMARY FIGURES 

results_rho = {rho_eps, rho_sla};
results_pval = {pval_eps, pval_sla};
results_h = {h_eps, h_sla};
results_title = strcat('regression coefficients ',{'(autoregressed SLA)','(normalized SLA)'});
results_save = strcat('./results/',{'SLAeps','SLAnorm'});

%{
% (0.05 & BONFERRONI)
% run for two different pvalue thresholds - uncorrected (U) and bonferroni
% corrected (B)
pval_thresh = [0.05 0.05/(length(depths)*length(varnames))];
pval_str = {'U','B'};

for p = 1:length(pval_thresh)
    for m = 1:length(results_rho)

        % heatmap summaries (H)
        plot_REG_results_h(depths,varnames,results_rho{m},results_pval{m}<pval_thresh(p));
        title(results_title{m});
        print(sprintf('%sH%s',results_save{m},pval_str{p}),'-dpng');

        % table summaries (T)
        plot_REG_table_h(depths,varnames,results_rho{m},results_pval{m}<pval_thresh(p));
        print(sprintf('%sT%s',results_save{m},pval_str{p}),'-dpng');

    end
end
%}

% BENJAMINI & HOCHBERG PROCEDURE
for m = 1:length(results_rho)
    
    % heatmap summaries (H)
    plot_REG_results_h(depths,varnames,results_rho{m},results_h{m});
    title(results_title{m});
    colormap(zerocmap(redbluecmap))
    print(sprintf('%sH%s',results_save{m},'_B&H'),'-dpng');
    
    % table summaries (T)
    plot_REG_table_h(depths,varnames,results_rho{m},results_h{m});
    print(sprintf('%sT%s',results_save{m},'_B&H'),'-dpng');
    
end

% sampling time histograms
fig = figure('Units','inches');
fig.Position(3:4) = [6.5 8];
for k = 1:length(varnames)
    subplot(4,4,k);
    plot(dt_bins,dt_frac(:,:,k),'.-')
    title(varnames{k});
    xlim([0 max(dt_bins)]);
    %ylim([0 0.7]);
    set(gca,'XTick',0:30:max(dt_bins));
end
subplot_labels('','days between consecutive samples','number of sample points');
print('./results/sampletimes','-dpng');






