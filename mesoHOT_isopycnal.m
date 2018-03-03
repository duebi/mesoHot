% mesoHot_isopycnal
%
% Includes several isopycnal analysis in separate cells:
%
% 1. Isopycnal depth vs SLA
% 2. Second horizontal derivative of SLA at Station ALOHA
% 3. Average HOT potential density profile
% 4. Isopycnal concentration anomalies
%
% Benedetto Barone - Feb 2018


%% 1. Analysis on isopycnal depth change with SLA
load allHawaii
A = csvread('pden0_1000.txt',2,0);
A(:,5) = [];
hotd.crn = A(1:501:end,1);
hotd.date = datenum(1988,10,1) + A(1:501:end,2) -1;
hotd.press = 0:2:1000;
hotd.depth = -gsw_z_from_p(hotd.press,22.75);
hotd.r_sla = interp1(allHawaii.date -10/24,allHawaii.aloha,hotd.date); % uncorrected SLA
hotd.c_sla = interp1(allHawaii.date -10/24,allHawaii.deseas,hotd.date); % de-trended, no season, SLA
hotd.sig = reshape(A(:,4),501,232);
clear A

hotd.sig(243:end,18) = NaN; hotd.sig(243:end,25) = NaN;
hotd.sig(244:end,31) = NaN; hotd.sig(244:end,33) = NaN;
hotd.sig(104:end,182) = NaN;

hotd.msig = nanmean(hotd.sig,2);
hotd.sigz = NaN*hotd.sig;
for i = 1:232
    ind_gd = ~isnan(hotd.sig(:,i)) & hotd.press' >20;
    hotd.sigz(:,i) = interp1(sort(hotd.sig(ind_gd,i))+(1:sum(ind_gd))'*1e-8,hotd.depth(ind_gd),hotd.msig);
end

% Regressions on potential density vertical displacement

sig_m = NaN(501,1);
sig_sm = NaN(501,1);
sigR2 = NaN(501,1);
sigpval = NaN(501,1);
sig_mc = NaN(501,1);
sig_smc = NaN(501,1);
sigR2c = NaN(501,1);
sigpvalc = NaN(501,1);
sig_n = NaN(501,1);
sig_bc = NaN(501,1);
sig_sbc = NaN(501,1);
for i = 2:501;
    ntnan = ~isnan(hotd.sigz(i,:));
    sig_n(i) = sum(ntnan);
    [m_r,b_r,r_r,sm_r,sb_r] = lsqfitgm(hotd.r_sla(ntnan)',hotd.sigz(i,ntnan)); % uncorrected vs. uncorrected
    [m_c,b_c,r_c,sm_c,sb_c] = lsqfitgm(hotd.c_sla(ntnan)',hotd.sigz(i,ntnan)); % corrected vs. uncorrected
    sig_m(i) = m_r; sig_sm(i) = sm_r;
    [r_r,p_r] = corrcoef(hotd.r_sla(ntnan)'*m_r + b_r,hotd.sigz(i,ntnan));
    sigR2(i) = r_r(1,2)^2; sigpval(i) = p_r(1,2);
    sig_mc(i) = m_c; sig_smc(i) = sm_c;
    sig_bc(i) = b_c; sig_sbc(i) = sb_c;
    [r_c,p_c] = corrcoef(hotd.c_sla(ntnan)'*m_c + b_c,hotd.sigz(i,ntnan));
    sigR2c(i) = r_c(1,2)^2; sigpvalc(i) = p_c(1,2);
end

clf
subplot(1,4,1:2)
patch([sig_mc(2:end)-sig_smc(2:end); sig_mc(end:-1:2)+sig_smc(end:-1:2)],hotd.depth([2:end end:-1:2]),[0.85 0.85 0.85],'Edgecolor','none','Facealpha',1);
hold on, plot(sig_mc,hotd.depth,'k'), hold off
set(gca,'ydir','rev','box','on')
xlabel('Slope (m/cm)'), ylabel('Average isopycnal depth (m)')
subplot(1,4,3)
hold on, plot(sigR2c,hotd.depth,'k',sigR2,hotd.depth,'k--'), hold off
set(gca,'ydir','rev','box','on','YaxisLocation','right')
xlabel('R^2 of regression'), ylabel('Average isopycnal depth (m)')
lg = legend('corrected SLA','raw SLA'); set(lg,'Fontsize',18)
subplot(1,4,4)
hold on, plot(sig_n./length(hotd.sigz(1,:)),hotd.depth,'k.-'), hold off
set(gca,'ydir','rev','box','on','YaxisLocation','right')
xlabel('n'), ylabel('Average isopycnal depth (m)')

%% 2. Compute second derivative of SLA with finite difference from SLA maps
load allHawaii
% Compute d2SLA/dx2 and d2SLA/dy2
d2lat = allHawaii.lat(12:15);
d2lon = allHawaii.lon(13:16);
sla_tmp = allHawaii.sla(12:15,13:16,:);
xdist = vdist(22.75,-d2lon(3),22.75,-d2lon(2))/1000; % Distance between consecutive longitude points in km
ydist = vdist(d2lat(3),-158,d2lat(2),-158)/1000; % Distance between consecutive latitude points in km
d2x = (sla_tmp(:,3:end,:) - 2*sla_tmp(:,2:(end-1),:) + sla_tmp(:,1:(end-2),:))./(xdist^2); % Zonal second derivative
d2y = (sla_tmp(3:end,:,:) - 2*sla_tmp(2:(end-1),:,:) + sla_tmp(1:(end-2),:,:))./(ydist^2); % Meridional second derivative
% Interpolate value for Station ALOHA
lat_tmp = d2lat(2:3);
lon_tmp = d2lon(2:3);
d2sla_tmp  = d2x(2:3,:,:) + d2y(:,2:3,:);
d2sla = interp3(lon_tmp,lat_tmp,allHawaii.date,d2sla_tmp,158,22.75,allHawaii.date);
d2sla = squeeze(d2sla);
% Plot SLA vs d2SLA/dx2+d2SLA/dy2+
clear d2lat d2lon sla_tmp xdist ydist d2x d2y lat_tmp lon_tmp d2sla_tmp
scatter(allHawaii.deseas,d2sla,35,'Markerfacecolor','k','Markeredgecolor','none','Markerfacealpha',0.2)
xlabel('$\mathrm{SLA_{corr} (cm)}$','Interpreter','latex'), ylabel('$\mathrm{\partial^2 SLA/\partial x^2 + \partial^2 SLA/\partial y^2 \> (cm \> km^2)}$','Interpreter','latex')
set(gca,'box','on')

%% 3. Compute average HOT potential density profile (1993-2015)

load allHawaii
A = csvread('pden0_1000.txt',2,0);
A(:,5) = [];
hotd.crn = A(1:501:end,1);
hotd.date = datenum(1988,10,1) + A(1:501:end,2) -1;
hotd.press = 0:2:1000;
hotd.depth = -gsw_z_from_p(hotd.press,22.75);
hotd.r_sla = interp1(allHawaii.date -10/24,allHawaii.aloha,hotd.date); % uncorrected SLA
hotd.c_sla = interp1(allHawaii.date -10/24,allHawaii.deseas,hotd.date); % de-trended, no season, SLA
hotd.sig = reshape(A(:,4),501,232);
hotd.sig(hotd.sig==-9) = NaN;
clear A

hotd.sig(243:end,18) = NaN; hotd.sig(243:end,25) = NaN;
hotd.sig(244:end,31) = NaN; hotd.sig(244:end,33) = NaN;
hotd.sig(104:end,182) = NaN;
% Average pot density profile  
hotd.msig = nanmean(hotd.sig,2);
% Density levels for data interpolation (25 levels)
depth_grid = 10:10:200;
pden_grid = interp1(hotd.depth,hotd.msig,depth_grid);

plot(hotd.msig,hotd.depth,'k-',pden_grid,depth_grid,'ko')
set(gca,'ydir','rev')
ylim([0 200])
pdengrid = table(depth_grid',pden_grid');
pdengrid.Properties.VariableNames = {'depth','sigma'};

%% 4. Isopycnal anomalies of biogeochemical variables

clear
load allHawaii
load pdengrid
data = readHotDogs('nit_0_300m_sig.txt',1);
varname = data.Properties.VariableNames{end};
varunits = data.Properties.VariableUnits{end};
data(isnan(data.date) | eval(['isnan(data.' varname ')']),:) = [];
% Interpolate SLA on bottle data dates
data.r_sla = interp1(allHawaii.date -10/24,allHawaii.aloha,data.date); % uncorrected SLA
data.c_sla = interp1(allHawaii.date -10/24,allHawaii.deseas,data.date); % de-trended, no season, SLA
% Compute month
month = datevec(data.date);
data.month = month(:,2);
clear month
% Change obs <0 into 0
dtmp = data{:,varname};
dtmp(dtmp<0) = 0;
eval(['data.' varname ' = dtmp;']);
% From umol/kg to umol L
%eval(['data.' varname ' = data.' varname '.*(1000+data.sigma)/1000;']);

% isolate profiles and interpolate on density grid
id = 1;
datasig = [];
ls = length(pdengrid.sigma);
while id < height(data)
    ipr = data.cruise == data.cruise(id) & data.cast == data.cast(id);
    ipr = find(ipr);
    if length(ipr)>1
        datapr = data{ipr,{'sigma';varname;'depth'}};
        if sum(datapr(:,3)<200) >=5
            datapr=grpstats(datapr,datapr(:,1),'mean'); %take the average of duplicate sigma values
            tmp = [data.cruise(id)*ones(ls,1) data.cast(id)*ones(ls,1) data.date(id)*ones(ls,1) data.r_sla(id)*ones(ls,1) data.c_sla(id)*ones(ls,1) data.month(id)*ones(ls,1) pdengrid.depth pdengrid.sigma];
            varsig = interp1(datapr(:,1),datapr(:,2),pdengrid.sigma,'linear');
            depthsig = interp1(datapr(:,1),datapr(:,3),pdengrid.sigma,'linear'); % keep original measurement depth
            datasig = [datasig; [tmp varsig depthsig]];
            %plot(varsig,pdengrid.sigma,'rx',datapr(:,2),datapr(:,1),'k-o',pdengrid.sigma*0,pdengrid.sigma,'rs')
            %pause(0.1)
            clear varsig depthsig tmp datapr
        end
    end
    id = ipr(end)+1;
end
datasig = array2table(datasig);
datasig.Properties.VariableNames = {'cruise','cast','date','r_sla','c_sla','month','depth','sigma',varname,'odepth'};
datasig.anom = NaN(height(datasig),1);
% Compute isopycnal anomaly
for i = 1:20
    ind_dpt = datasig.depth == pdengrid.depth(i);
    datasig.anom(ind_dpt) = eval(['datasig.' varname '(ind_dpt) - nanmedian(datasig.' varname '(ind_dpt));']);
end


% Initialize final variables (USING ISOPYCNAL DEPTH)
    depth = pdengrid.depth; sigma = pdengrid.sigma; n = NaN(20,1); m = NaN(20,1); b = NaN(20,1); r = NaN(20,1); pval = NaN(20,1); h = NaN(20,1);   
    m_all = NaN(20,1); m_5m = NaN(20,1); m_5p = NaN(20,1); p75_5m = NaN(20,1); p75_5p = NaN(20,1); p25_5m = NaN(20,1); p25_5p = NaN(20,1);
    results = table(depth, sigma, n, m, b, r, pval, h, m_all, m_5m, m_5p, p25_5m, p75_5m, p25_5p, p75_5p);
    alldata = {};
    clear depth sigma ns m b r pval h m_all m_5m m_5p p25_5m p75_5m p25_5p p75_5p
for i = 1:20
    ind_dpt = datasig.depth == pdengrid.depth(i) & eval(['~isnan(datasig.' varname ')']);
    % Isolate data from the depth layer & count data points
    results.n(i) = sum(ind_dpt);
    ddata = datasig(ind_dpt,:);
    alldata{i} = ddata;
    % Model I linear regressions
    [m_c,b_c,r_c,sm_c,sb_c] = lsqfity(ddata.c_sla,eval(['ddata.' varname])); % corrected vs. uncorrected
    % Coefficients of determination and significance
    [r_c2,p_c2] = corrcoef(ddata.c_sla*m_c+b_c,eval(['ddata.' varname]));
    results.m(i) = m_c;
    results.b(i) = b_c;
    results.r(i) = r_c;
    results.pval(i,:) = p_c2(2);
    results.h(i) = p_c2(2)<0.05;    
    % Find median values and percentiles for SLAcorr >5cm and <-5cm
    ind5m = ddata.c_sla < -5;
    ind5p = ddata.c_sla > 5;
    results.m_all(i) = nanmedian(eval(['ddata.' varname]));
    results.m_5m(i) = nanmedian(eval(['ddata.' varname '(ind5m)']));
    results.m_5p(i) = nanmedian(eval(['ddata.' varname '(ind5p)']));
    results.p75_5m(i) = prctile(eval(['ddata.' varname '(ind5m)']),75);
    results.p75_5p(i) = prctile(eval(['ddata.' varname '(ind5p)']),75);
    results.p25_5m(i) = prctile(eval(['ddata.' varname '(ind5m)']),25);
    results.p25_5p(i) = prctile(eval(['ddata.' varname '(ind5p)']),25);
    % Rank sum test
    [P,H] = ranksum(eval(['ddata.' varname '(ind5m)']),eval(['ddata.' varname '(ind5p)']));
    rks(i,:) = [P H];
    
    % Clear temporary variables
    clear ind_dpt P H
    clear m_c r_c b_c sm_c sb_c
    clear r_c2 p_c2 
    clear ind5m ind5p
end

% Initialize final variables (USING ORIGINAL DEPTH)
    depth = pdengrid.depth; sigma = pdengrid.sigma; n = NaN(20,1); m = NaN(20,1); b = NaN(20,1); r = NaN(20,1); pval = NaN(20,1); h = NaN(20,1);   
    m_all = NaN(20,1); m_5m = NaN(20,1); m_5p = NaN(20,1); p75_5m = NaN(20,1); p75_5p = NaN(20,1); p25_5m = NaN(20,1); p25_5p = NaN(20,1);
    results2 = table(depth, sigma, n, m, b, r, pval, h, m_all, m_5m, m_5p, p25_5m, p75_5m, p25_5p, p75_5p);
    alldata2 = {};
    clear depth sigma ns m r pval h m_all m_5m m_5p p25_5m p75_5m p25_5p p75_5p
for i = 1:20
    ind_dpt = datasig.odepth <= pdengrid.depth(i)+5 & datasig.odepth >= pdengrid.depth(i)-5 & ~isnan(datasig.anom);
    % Isolate data from the depth layer & count data points
    results2.n(i) = sum(ind_dpt);
    ddata = datasig(ind_dpt,:);
    alldata2{i} = ddata;
    % Model I linear regressions
    [m_c,b_c,r_c,sm_c,sb_c] = lsqfity(ddata.c_sla,ddata.anom); % corrected vs. uncorrected
    % Coefficients of determination and significance
    [r_c2,p_c2] = corrcoef(ddata.c_sla*m_c+b_c,ddata.anom);
    results2.m(i) = m_c;
    results2.b(i) = b_c;
    results2.r(i) = r_c;
    results2.pval(i,:) = p_c2(2);
    results2.h(i) = p_c2(2)<0.05;    
    % Find median values and percentiles for SLAcorr >5cm and <-5cm
    ind5m = ddata.c_sla < -5;
    ind5p = ddata.c_sla > 5;
    results2.m_all(i) = nanmedian(ddata.anom);
    results2.m_5m(i) = nanmedian(ddata.anom(ind5m));
    results2.m_5p(i) = nanmedian(ddata.anom(ind5p));
    results2.p75_5m(i) = prctile(ddata.anom(ind5m),75);
    results2.p75_5p(i) = prctile(ddata.anom(ind5p),75);
    results2.p25_5m(i) = prctile(ddata.anom(ind5m),25);
    results2.p25_5p(i) = prctile(ddata.anom(ind5p),25);
    % Rank sum test
    %[P,H] = ranksum(ddata.anom(ind5m),ddata.anom(ind5p));
    %rks2(i,:) = [P H];
    
    % Clear temporary variables
    clear ind_dpt P H
    clear m_c r_c b_c sm_c sb_c
    clear r_c2 p_c2 
    clear ind5m ind5p
end

clf
subplot(1,3,1)
patch([results.p75_5p - results.m_all; results.p25_5p(end:-1:1)- results.m_all(end:-1:1)],[results.depth; results.depth(end:-1:1)],[0.6 0.4 0.4],'Edgecolor','none','Facealpha',0.2)
hold on, patch([results.p75_5m - results.m_all; results.p25_5m(end:-1:1)- results.m_all(end:-1:1)],[results.depth; results.depth(end:-1:1)],[0.4 0.4 0.6],'Edgecolor','none','Facealpha',0.2); hold off
hold on,ln12 = plot(results.m_5p - results.m_all,results.depth,'r',results.m_5m - results.m_all,results.depth,'b','Linewidth',2); hold off
hold on, plot([0 0],ylim,'k--'),hold off
tmpl = xlim;
%hold on, plot(tmpl(2)+0*pdengrid.depth(logical(rks(:,2))),pdengrid.depth(logical(rks(:,2))),'k^','Markersize',12,'MarkerFaceColor','k'), hold off
set(gca,'ydir','rev','ytick',pdengrid.depth(2:2:end))
set(gca,'YTickLabel',pdengrid.sigma(2:2:end),'YAxisLocation','right')
lg = legend(ln12,'SLA_{corr} > 5 cm','SLA_{corr} < -5 cm'); set(lg,'Fontsize',18);
ylabel('Average isopycnal depth (m)')
box('on')
ylim([0 200])
xlabel([varname ' anomaly (' varunits ')'])
subplot(1,3,2)
patch([results2.p75_5p; results2.p25_5p(end:-1:1)],[results2.depth; results2.depth(end:-1:1)],[0.6 0.4 0.4],'Edgecolor','none','Facealpha',0.2)
hold on, patch([results2.p75_5m; results2.p25_5m(end:-1:1)],[results2.depth; results2.depth(end:-1:1)],[0.4 0.4 0.6],'Edgecolor','none','Facealpha',0.2); hold off
hold on,ln12 = plot(results2.m_5p,results2.depth,'r',results2.m_5m,results2.depth,'b',results2.m_all,results2.depth,'k--','Linewidth',2); hold off
hold on, plot([0 0],ylim,'k--'),hold off
tmpl = xlim;
%hold on, plot(tmpl(2)+0*pdengrid.depth(logical(rks2(:,2))),pdengrid.depth(logical(rks2(:,2))),'k^','Markersize',12,'MarkerFaceColor','k'), hold off
set(gca,'ydir','rev','ytick',pdengrid.depth(2:2:end))
lg = legend(ln12,'SLA_{corr} > 5 cm','SLA_{corr} < -5 cm','median depth anomaly'); set(lg,'Fontsize',18);
ylabel('Depth (m)')
box('on')
ylim([0 200])
xlabel([varname ' isopycnal anomaly (' varunits ')'])
subplot(1,3,3)
ind100 = find(pdengrid.depth == 100);
datareg = alldata{ind100};
[m_c,b_c,r_c,sm_c,sb_c] = lsqfity(datareg.c_sla,eval(['datareg.' varname])); 
scatter(datareg.c_sla,eval(['datareg.' varname]),50,'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFacealpha',0.3)
hold on, plot(datareg.c_sla,datareg.c_sla*m_c+b_c,'r-','linewidth',2), hold off
set(gca,'box','on')
ylabel([varname '(' varunits ')'])