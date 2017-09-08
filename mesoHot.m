% script mesoHot.m
%
% Regression analyses of biogeochemical measurements from the Hawaii Ocean
% Time-series (HOT) versus Sea Level Anomaly (SLA) from CMEMS.
% Statistical analyses are done on 8 depth levels in the upper 200 m
% coinciding with the standard sampling depths of HOT.
%
% LINEAR REGRESSIONS:
% Biogeochemical variables are regressed against both uncorrected SLA and 
% corrected SLA (i.e. the one with trend and seasonal cycle removed).
% A third regression is computed between corrected SLA and the monthly 
% anomaly of the biogeochemical variable.
% Regressions are of the Model II geometric mean type.
%
% COMPARISON OF DISTRIBUTION WHEN SLA > 5cm & SLA < -5cm
% Median, 25th and 75th percentiles are computed for the cases of high and
% low SLA (the corrected one is used in this analysis). The distributions
% of the variable in the two different cases are then compared using a rank
% sum test to check for significant differences (p<0.05).
%
% REQUIREMENTS:
% - Matlab ver 2013b or above
% - Gibbs SeaWater (GSW) Oceanographic Toolbox of TEOS-10 available rom http://www.TEOS-10.org
% - lsqfitgm.m, lsqfitx.m and lsqfity.m available from http://www.mbari.org/index-of-downloadable-files/ (by E.T. Peltzer) 
% - readHotDogs.m
% 
% The parameter for the analysis is chosen by renaming the variable <biovar>
% to be one of the following:
% acar, alk, atp, bcar, boxy, but, chl, chlb, chlc, chlda, diad, dic, doc,
% don, dop, dvchl, euk, fuco, hbact, hex, hplc, l12, lln, llp, lut, mvchl,
% mvdv, nit, pc, pe1-, pe4, pe5, peri, ph, phos, pn, pp, peas, pro, sil,
% syn, tdn, tdp, theta, viol, zeax
%
% Benedetto Barone - Sep 7, 2017 - Version 0.1

biovar = 'euk'; % biogeochemical variable name (INPUT)

% Load altimetry data
load ./data/allHawaii 
% Extract data for selected variable
data = readHotDogs(['./data/' biovar '_0_200m.txt'],1);
% Compute pigment ratio in the case of mvdv
if strcmp(data.Properties.VariableNames(7),'mvchla') & strcmp(data.Properties.VariableNames(8),'dvchla')
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
% Compute month (to be used for monthly anomalies)
month = datevec(data.date);
data.month = month(:,2);
clear month

% Define depths for regression analyses 
depths = [5 25 45 75 100 125 150 175];

% Initialize variables for the results of the regressions
ms = NaN(8,3); % slopes
r2s = NaN(8,3); % R2s
hs = NaN(8,3); % significance
ns = NaN(8,1); % number of data points

% Initialize variables for the results of the regressions
var5m = NaN(8,1); % median SLAcorr <-5cm
p75_5m = NaN(8,1); % 75th percentile SLAcorr <-5cm
p25_5m = NaN(8,1); % 25th percentile SLAcorr <-5cm
var5p = NaN(8,1); % median SLAcorr >5cm
p75_5p = NaN(8,1);
p25_5p = NaN(8,1);
rks = NaN(8,2); % ranksum test results

% Extract data at each depth, compute regression, compute statistics for
% SLA > 5cm & SLA < -5cm
for i = 1:8
    ind_dpt = data.depth >= depths(i)-3 & data.depth <= depths(i)+3;
    % increase depth layer width for 100m for tdn, tdp, don, and dop
    if (strcmp(varname,'tdn') | strcmp(varname,'tdp') | strcmp(varname,'don') | strcmp(varname,'dop')) & depths(i)==100
        ind_dpt = data.depth >= depths(i)-7 & data.depth <= depths(i)+7;
    end
    % Isolate data from the depth layer & count data points
    ns(i) = sum(ind_dpt);
    ddata = data(ind_dpt,:);
    % Transform data into monthly anomalies
    month_m = NaN(12,1);
    um = unique(ddata.month);
    month_m(um) = grpstats(eval(['ddata.' varname]),ddata.month,'mean');
    data_anom = eval(['ddata.' varname]) - month_m(ddata.month);
    % Model II linear regressions
    [m_r,b_r,r_r,sm_r,sb_r] = lsqfitgm(ddata.r_sla,eval(['ddata.' varname])); % uncorrected vs. uncorrected
    [m_c,b_c,r_c,sm_c,sb_c] = lsqfitgm(ddata.c_sla,eval(['ddata.' varname])); % corrected vs. uncorrected
    [m_cc,b_cc,r_cc,sm_cc,sb_cc] = lsqfitgm(ddata.c_sla,data_anom); % corrected vs. corrected
    % Coefficients of determination and significance
    [r_r2,p_r2] = corrcoef(ddata.r_sla*m_r+b_r,eval(['ddata.' varname]));
    [r_c2,p_c2] = corrcoef(ddata.c_sla*m_c+b_c,eval(['ddata.' varname]));
    [r_cc2,p_cc2] = corrcoef(ddata.c_sla*m_cc+b_cc,data_anom);
    ms(i,:) = [m_r m_c m_cc];
    r2s(i,:) = [r_r2(2)^2 r_c2(2)^2 r_cc2(2)^2];
    hs(i,:) = [p_r2(2)<0.05 p_c2(2)<0.05 p_cc2(2)<0.05];    
    % Find median values and percentiles for SLAcorr >5cm and <-5cm
    ind5m = ddata.c_sla < -5;
    ind5p = ddata.c_sla > 5;
    var5m(i) = nanmedian(eval(['ddata.' varname '(ind5m)']));
    var5p(i) = nanmedian(eval(['ddata.' varname '(ind5p)']));
    p75_5m(i) = prctile(eval(['ddata.' varname '(ind5m)']),75);
    p75_5p(i) = prctile(eval(['ddata.' varname '(ind5p)']),75);
    p25_5m(i) = prctile(eval(['ddata.' varname '(ind5m)']),25);
    p25_5p(i) = prctile(eval(['ddata.' varname '(ind5p)']),25);
    % Rank sum test
    [P,H] = ranksum(eval(['ddata.' varname '(ind5m)']),eval(['ddata.' varname '(ind5p)']));
    rks(i,:) = [P H];
    % Clear temporary variables
    clear ind_dpt ddata month_m um data_anom P H
    clear m_r m_c m_cc r_r r_c r_cc b_r b_c b_cc sm_r sm_c sm_cc sb_r sb_c sb_cc
    clear r_r2 r_c2 r_cc2 p_r2 p_c2 p_cc2
end

% Plot median profiles for SLA >5cm and <-5cm 
figure
clf
subplot(1,3,1:2)
patch([p75_5p; p25_5p(end:-1:1)],[depths'; depths(end:-1:1)'],[0.6 0.4 0.4],'Edgecolor','none','Facealpha',0.2)
hold on, patch([p75_5m; p25_5m(end:-1:1)],[depths'; depths(end:-1:1)'],[0.4 0.4 0.6],'Edgecolor','none','Facealpha',0.2); hold off
hold on,ln12 = plot(var5p,depths,'r',var5m,depths,'b'); hold off
tmpl = xlim;
hold on, plot(tmpl(2)+0*depths(logical(rks(:,2))),depths(logical(rks(:,2))),'k^','Markersize',12,'MarkerFaceColor','k'), hold off
set(gca,'ydir','rev','ytick',[5 25 45 75 100 125 150 175])
lg = legend(ln12,'SLA_{corr} > 5 cm','SLA_{corr} < -5 cm'); set(lg,'Fontsize',18);
ylabel('Depth (m)')
box('on')
ylim([0 180])
xlabel([varname ' (' varunits ')'])
title('Comparison of extreme quartiles')
% Report results of linear regressions 
subplot(1,3,3)
plot(0*depths(logical(hs(:,2)))+0.2,depths(logical(hs(:,2))),'k.','Markersize',40)
for i = 1:8
    if hs(i,2) == 1
       if sign(ms(i,2)) > 0
           sgn = '+';
       else
           sgn = '-';
       end
       text(0.3,depths(i),[sgn '(' num2str(r2s(i,2),2) ')']) 
    end
end
xlim([0 1]),ylim([0 180])
set(gca,'ydir','rev','box','on')
title('Lin. Regressions')

clear tmpl sgn ln12 lg i