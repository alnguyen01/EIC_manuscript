%% read in time series array 
coh=h5read('../geo/geo_avgSpatialCoh.h5','/coherence');
ts=h5read('../geo/geo_timeseries_ERA5_ramp_demErr.h5','/timeseries');
dates=h5read('../geo/geo_timeseries_ERA5_ramp_demErr.h5','/date');

doys=zeros(length(dates),1);

for i=1:length(dates)
    t = datetime(str2double(dates(i)),'ConvertFrom','yyyymmdd');
    doys(i,1) = day(t,'dayofyear');
end

year = t.Year;
maskstring = strcat("../../mintpy2018_2025_NW_drew/ReSALT/otsu_",num2str(year),".tif");

[mask,ref] = readgeoraster(maskstring);

coh = coh.*mask';
ts = ts.*mask';
nr = size(ts,1);
naz = size(ts,2);

[theta_inc,ff]=readgeoraster("../../warped_files/NW_drew/warped_S1AA_20180907T165937_20180919T165937_VVP012_INT40_G_weF_A162_inc_map.tif");

%% locate pixels of interest
pois = readtable("../../warped_files/NW_drew/test_pixels.xls",'VariableNamingRule','preserve'); %read in table of pixels of interest
longitude = pois.lon;
latitude = pois.lat;
entries = length(longitude);
names = pois.POI;

xfirst = str2double(h5readatt('../geo/geo_temporalCoherence.h5','/','X_FIRST'));
yfirst = str2double(h5readatt('../geo/geo_temporalCoherence.h5','/','Y_FIRST'));

xstep = str2double(h5readatt('../geo/geo_temporalCoherence.h5','/','X_STEP'));
ystep = str2double(h5readatt('../geo/geo_temporalCoherence.h5','/','Y_STEP'));

pixelX = zeros(entries,1);
pixelY = zeros(entries,1);

for i = 1:entries
    pixelX(i,1) = ceil(abs((xfirst+(-longitude(i,1)))/xstep));
    pixelY(i,1) = ceil(abs((yfirst-latitude(i,1))/ystep));
end

%%
names_ts = strcat(names,' time series');
names_fits = strcat(names,' fit');

names_extrap = strcat(names, ' extrapolation');
naddt_string = "NADDT";

vec = [names_ts names_fits]';
originallegendnames = vec(:);
originallegendnames{end+1} = 'NADDT';

vecextrap = [names_ts names_extrap]';
extraplegendnames = vecextrap(:);
extraplegendnames{end+1} = 'NADDT';

colors_full = {'red','blue','green','magenta','black','yellow'}; %specific for 6 pixels
colors = {'r--','b--','g--','m--','k--','y--'}; %specific for 6 pixels

startofyear = datetime(year, 1,1);
dates_datetime = startofyear + doys;
%%
[addt] = addt_curve(year);

freezeidx = find(addt==1);
beforefreeze = (min(freezeidx))-1; %1st day above max ADDT; assuming freezing
freezedate = startofyear + caldays(beforefreeze); %check date relative to timeseries dates

thawidx = find(addt>0);
thawfirst = thawidx(1);
thawdate = startofyear + caldays(thawfirst);

thawlength = days(freezedate-thawdate)+1;

num_dates = length(dates);
dates_datetime_full = startofyear + day(thawdate-1:freezedate-1,'dayofyear');
dates_full_doy = (thawfirst:beforefreeze)';

%% copy over from other if messed up 
if freezedate > dates_datetime(num_dates)
    disp("Maximum ADDT date greater than last time series date. Use all timestamps.")
    addtseries=addt(doys); %not using square root b/c this is done in addt_curve.m
    addtfull = addt(thawfirst:beforefreeze);
    ts_thaw=squeeze(ts(:,:,:))-squeeze(ts(:,:,1));

    yfitaddt = zeros(nr,naz,num_dates);
    Fit_addt = zeros(nr,naz,2);
    yfit = zeros(nr,naz,thawlength);
    ts_thaw(isnan(ts_thaw))=0;
    for j = 1:nr
        disp(j)
        for k = 1:naz
            Fit_addt(j,k,:) = polyfit(addtseries,ts_thaw(j,k,:),1);
            yfitaddt(j,k,:) = (Fit_addt(j,k,1)*addtseries)+Fit_addt(j,k,2);
            yfit(j,k,:) = (Fit_addt(j,k,1)*addtfull)+Fit_addt(j,k,2);
        end
    end
    yfitaddt =(yfitaddt*100)./(cos(theta_inc')); %cm
    yfit = (yfit*100)./(cos(theta_inc')); %cm
    ts_thaw =(ts_thaw*100)./(cos(theta_inc')); %cm

    ts_thaw_extrap = ts_thaw(:,:,:) - yfit(:,:,1);
    yfit(:,:,:) = yfit(:,:,:) - yfit(:,:,1);

    E = yfit(:,:,size(yfit,3)).*-1;

    NegativeE = E;

    E(E<=0) = NaN; %masking out waterMask pixels and negative E pixels (heaving) to NaN; unreliable

    yfit_InSAR_dates = zeros(size(doys,1),1);
    
    for i = 1:size(dates_full_doy,1)
        for j = 1:size(doys,1)
            if dates_full_doy(i) == doys(j)
                yfit_InSAR_dates(j) = i;
            else
            end
        end
    end

    figure
    tiledlayout(2,2,'TileSpacing','tight')
    ax1 = nexttile;
    imagesc((coh)');
    colormap(ax1,'gray');
    axis on
    hold on
    for i = 1:entries
        h = scatter(pixelX(i),pixelY(i),100,'white','filled','o','LineWidth',4);
        h = plot(pixelX(i),pixelY(i),'Color',colors_full{i},'Linewidth', 2,'MarkerSize',10);
        h.Marker = 'x';
%         text(pixelX(i)+15,pixelY(i)+15,names(i),'Color',colors_full{i},'FontSize',20)
    end
    title('Masked Average Spatial Coherence')
    set(gca,'FontSize',20)
    colorbar

    ax2 = nexttile;
    h =imagesc(E');
    colormap(ax2, 'parula');
    set(gca,'color',0.5*[1 1 1])
    set(h, 'AlphaData', ~isnan(E'))
    axis on
    hold on
    for i = 1:entries
        h = scatter(pixelX(i),pixelY(i),100,'white','filled','o','LineWidth',4);
        h = plot(pixelX(i),pixelY(i),'Color',colors_full{i},'Linewidth', 2,'MarkerSize',10);
        h.Marker = 'x';
%         text(pixelX(i)+15,pixelY(i)+15,names(i),'Color',colors_full{i},'FontSize',20)
    end
    title([num2str(year),' Seasonal Subsidence from Extrapolated Fit'])
    set(gca,'FontSize',20)
    c = colorbar;
    xlabel(c, 'E [cm]')

    stats_array = squeeze(E);
    row_stats_array = (stats_array(:));
    percent5 = prctile(row_stats_array,5); %colorbar lims
    percent95 = prctile(row_stats_array,95);
    clim([percent5 percent95])

    nexttile
    for i = 1:entries
        yyaxis left
        scatter(dates_datetime,squeeze(ts_thaw(pixelX(i),pixelY(i),:)),100,colors_full{i},'filled')
        hold on
        plot(dates_datetime,squeeze(yfitaddt(pixelX(i),pixelY(i),:)),colors{i},'LineWidth',3)
        yyaxis right
        set(gca,'YDir','reverse')
        plot(dates_datetime_full,addt(dates_full_doy))
    end
    title([num2str(year),' Seasonal Thaw Subsidence Original Fit'])
    xlabel('Date')
    yyaxis left
    ylabel('Elevation (cm)')
    ylim_min = min(min(min(min(yfitaddt(pixelX,pixelY,:)))),min(min(min(yfit(pixelX,pixelY,:)))));
    ylim_max = max(max(max(max(yfitaddt(pixelX,pixelY,:)))),max(max(max(yfit(pixelX,pixelY,:)))));
    ylim([ylim_min,ylim_max])
    yyaxis right
    ylabel('NADDT (-)')
%     legend(originallegendnames,'Location','bestoutside')
    legend('','','','','','','','','InSAR time series','Fit','','','NADDT',"Location","southwest")
    set(gca,'FontSize',20)
    set(gcf,'Position',[100 100 1500 1000])
    
    nexttile

    for i = 1:entries
        yyaxis left
        scatter(dates_datetime,squeeze(ts_thaw_extrap(pixelX(i),pixelY(i),:)),100,colors_full{i},'filled')
        hold on
        plot(dates_datetime_full,squeeze(yfit(pixelX(i),pixelY(i),:)),colors{i},'LineWidth',3,'LineStyle',':')
        yyaxis right
        set(gca,'YDir','reverse')
        plot(dates_datetime_full,addt(dates_full_doy))
    end
    title([num2str(year),' Seasonal Thaw Subsidence Extrapolated Fit'])
    xlabel('Date')
    yyaxis left
    ylabel('Elevation (cm)')
    ylim([ylim_min,ylim_max])
    yyaxis right
    ylabel('NADDT (-)')
%     legend(extraplegendnames,'Location','bestoutside')
    legend('','','','','','','','','InSAR time series','Extrapolation','','','NADDT',"Location","southwest")
    set(gca,'FontSize',20)
    set(gcf,'Position',[100 100 1500 1000])

    %AGU ts plot
    figure
    tiledlayout(1,2)
    nexttile
    for i = 1:entries
        yyaxis left
        scatter(dates_datetime,squeeze(ts_thaw(pixelX(i),pixelY(i),:)),100,colors_full{i},'filled')
        hold on
        plot(dates_datetime,squeeze(yfitaddt(pixelX(i),pixelY(i),:)),colors{i},'LineWidth',3)
        yyaxis right
        set(gca,'YDir','reverse')
        plot(dates_datetime_full,addt(dates_full_doy))
    end
    title([num2str(year),' Seasonal Thaw Subsidence Original Fit'])
    xlabel('Date')
    yyaxis left
    ylabel('Elevation (cm)')
    ylim_min = min(min(min(min(yfitaddt(pixelX,pixelY,:)))),min(min(min(yfit(pixelX,pixelY,:)))));
    ylim_max = max(max(max(max(yfitaddt(pixelX,pixelY,:)))),max(max(max(yfit(pixelX,pixelY,:)))));
    ylim([ylim_min,ylim_max])
    yyaxis right
    ylabel('NADDT (-)')
%     legend(originallegendnames,'Location','bestoutside')
    legend('','','','','','','','','InSAR time series','Fit','','','NADDT',"Location","southwest")
    set(gca,'FontSize',20)
    set(gcf,'Position',[100 100 1500 1000])
    
    nexttile

    for i = 1:entries
        yyaxis left
        scatter(dates_datetime,squeeze(ts_thaw_extrap(pixelX(i),pixelY(i),:)),100,colors_full{i},'filled')
        hold on
        plot(dates_datetime_full,squeeze(yfit(pixelX(i),pixelY(i),:)),colors{i},'LineWidth',3,'LineStyle',':')
        yyaxis right
        set(gca,'YDir','reverse')
        plot(dates_datetime_full,addt(dates_full_doy))
    end
    title([num2str(year),' Seasonal Thaw Subsidence Extrapolated Fit'])
    xlabel('Date')
    yyaxis left
    ylabel('Elevation (cm)')
    ylim([ylim_min,ylim_max])
    yyaxis right
    ylabel('NADDT (-)')
%     legend(extraplegendnames,'Location','bestoutside')
    legend('','','','','','','','','InSAR time series','Extrapolation','','','NADDT',"Location","southwest")
    set(gca,'FontSize',20)
    set(gcf,'Position',[100 100 1500 1000])
    f = gcf;
    f.WindowState = 'maximized';
    aguplot = strcat(num2str(year),'_ts.png');
    exportgraphics(f,aguplot)
else
    disp("Maximum ADDT date less than last time series date. Grabbing date fit prior as E.")
    finaldate = max(find(dates_datetime<freezedate));
    addtseries=addt(doys(1:finaldate));
    addtfull = addt(thawfirst:beforefreeze);

    ts_thaw=squeeze(ts(:,:,1:finaldate))-squeeze(ts(:,:,1));

    yfitaddt = zeros(nr,naz,finaldate);
    Fit_addt = zeros(nr,naz,2);
    yfit = zeros(nr,naz,thawlength);
    ts_thaw(isnan(ts_thaw))=0;
    for j = 1:nr
        disp(j)
        for k = 1:naz
            Fit_addt(j,k,:) = polyfit(addtseries,ts_thaw(j,k,:),1);
            yfitaddt(j,k,:) = (Fit_addt(j,k,1)*addtseries)+Fit_addt(j,k,2);
            yfit(j,k,:) = (Fit_addt(j,k,1)*addtfull)+Fit_addt(j,k,2);
        end
    end
    yfitaddt =(yfitaddt*100)./(cos(theta_inc')); %cm
    yfit = (yfit*100)./(cos(theta_inc')); %cm
    ts_thaw =(ts_thaw*100)./(cos(theta_inc')); %cm

    yfit(:,:,:) = yfit(:,:,:) - yfit(:,:,1);

    E = yfit(:,:,size(yfit,3)).*-1;

    NegativeE = E;

    E(E<=0) = NaN; %masking out waterMask pixels and negative E pixels (heaving) to NaN; unreliable

    yfit_InSAR_dates = zeros(size(doys(1:finaldate),1));

    for i = 1:size(dates_full_doy,1)
        for j = 1:size(doys,1)
            if dates_full_doy(i) == doys(j)
                yfit_InSAR_dates(j) = i;
            else
            end
        end
    end

    figure
    tiledlayout(2,2,'TileSpacing','tight')
    ax1 = nexttile;
    imagesc((coh.*mask')');
    colormap(ax1,'gray');
    axis on
    hold on
    for i = 1:entries
        h = scatter(pixelX(i),pixelY(i),100,'white','filled','o','LineWidth',4);
        h = plot(pixelX(i),pixelY(i),'Color',colors_full{i},'Linewidth', 2,'MarkerSize',10);
        h.Marker = 'x';
%         text(pixelX(i)+15,pixelY(i)+15,names(i),'Color',colors_full{i},'FontSize',20)
    end
    title('Masked Average Spatial Coherence')
    set(gca,'FontSize',20)
    colorbar

    ax2 = nexttile;
    h =imagesc(E');
    colormap(ax2, 'parula');
    set(gca,'color',0.5*[1 1 1])
    set(h, 'AlphaData', ~isnan(E'))
    axis on
    hold on
    for i = 1:entries
        h = scatter(pixelX(i),pixelY(i),100,'white','filled','o','LineWidth',4);
        h = plot(pixelX(i),pixelY(i),'Color',colors_full{i},'Linewidth', 2,'MarkerSize',10);
        h.Marker = 'x';
%         text(pixelX(i)+15,pixelY(i)+15,names(i),'Color',colors_full{i},'FontSize',20)
    end
    title([num2str(year),' Seasonal Subsidence from Extrapolated Fit'])
    set(gca,'FontSize',20)
    c = colorbar;
    xlabel(c, 'E [cm]')

    stats_array = squeeze(E);
    row_stats_array = (stats_array(:));
    percent5 = prctile(row_stats_array,5); %colorbar lims
    percent95 = prctile(row_stats_array,95);
    clim([percent5 percent95])

    nexttile
    for i = 1:entries
        yyaxis left
        scatter(dates_datetime,squeeze(ts_thaw(pixelX(i),pixelY(i),:)),100,colors_full{i},'filled')
        hold on
        plot(dates_datetime,squeeze(yfitaddt(pixelX(i),pixelY(i),:)),colors{i},'LineWidth',3)
        yyaxis right
        set(gca,'YDir','reverse')
        plot(dates_datetime_full,addt(dates_full_doy))
    end
    title([num2str(year),' Seasonal Thaw Subsidence Original Fit'])
    xlabel('Date')
    yyaxis left
    ylabel('Elevation (cm)')
    ylim_min = min(min(min(min(yfitaddt(pixelX,pixelY,:)))),min(min(min(yfit(pixelX,pixelY,:)))));
    ylim_max = max(max(max(max(yfitaddt(pixelX,pixelY,:)))),max(max(max(yfit(pixelX,pixelY,:)))));
    ylim([ylim_min,ylim_max])
    yyaxis right
    ylabel('NADDT (-)')
%     legend(originallegendnames,'Location','bestoutside')
    legend('','','','','','','','','InSAR time series','Fit','','','NADDT',"Location","southwest")
    set(gca,'FontSize',20)
    set(gcf,'Position',[100 100 1500 1000])
    
    nexttile

    for i = 1:entries
        yyaxis left
        scatter(dates_datetime,squeeze(ts_thaw_extrap(pixelX(i),pixelY(i),:)),100,colors_full{i},'filled')
        hold on
        plot(dates_datetime_full,squeeze(yfit(pixelX(i),pixelY(i),:)),colors{i},'LineWidth',3,'LineStyle',':')
        yyaxis right
        set(gca,'YDir','reverse')
        plot(dates_datetime_full,addt(dates_full_doy))
    end
    title([num2str(year),' Seasonal Thaw Subsidence Extrapolated Fit'])
    xlabel('Date')
    yyaxis left
    ylabel('Elevation (cm)')
    ylim([ylim_min,ylim_max])
    yyaxis right
    ylabel('NADDT (-)')
%     legend(extraplegendnames,'Location','bestoutside')
    legend('','','','','','','','','InSAR time series','Extrapolation','','','NADDT',"Location","southwest")
    set(gca,'FontSize',20)
    set(gcf,'Position',[100 100 1500 1000])
end
%% error
RMSE = zeros(nr,naz);
RMSE_extrap = zeros(nr,naz);

for j = 1:nr
    disp(j)
    for k= 1:naz
        RMSE(j,k) = sqrt(mean((yfitaddt(j,k,:)-ts_thaw(j,k,:)).^2));
        RMSE_extrap(j,k) = sqrt(mean((yfit(j,k,yfit_InSAR_dates)-ts_thaw(j,k,:)).^2));
    end
end

RMSE(RMSE==0) = NaN;
RMSE_extrap(RMSE_extrap==0)=NaN;
%% RMSE figure
figure
tiledlayout(1,2)
nexttile
h = imagesc(RMSE');
set(gca,'color',0.5*[1 1 1])
set(h, 'AlphaData', ~isnan(RMSE'))
title('Original Fit RMSE')
c = colorbar;
xlabel(c, 'RMSE [cm]')

stats_array = squeeze(RMSE);
row_stats_array = (stats_array(:));
percent5 = prctile(row_stats_array,5); %colorbar lims
percent95 = prctile(row_stats_array,95);
clim([percent5 percent95])

nexttile
h = imagesc(RMSE_extrap');
set(gca,'color',0.5*[1 1 1])
set(h, 'AlphaData', ~isnan(RMSE_extrap'))
title('Extrapolated Fit RMSE')
c = colorbar;
xlabel(c, 'RMSE [cm]')

stats_array = squeeze(RMSE_extrap);
row_stats_array = (stats_array(:));
percent5 = prctile(row_stats_array,5); %colorbar lims
percent95 = prctile(row_stats_array,95);
clim([percent5 percent95])

%% BELOW IS ALT + EIC CALC
%% calculate ALT
[deform] = E';
deform(deform==0)=NaN; %setting pixels w/ 0 equal to NaN

%NOTE!!!~~~ Keeping values for parameters default from Liu et al. 2012...
%May potentially be different depending on datasets being used and area
%covered
err_SaturationFraction = 0.01;
%SandFraction = 36.6; %old
SandFraction = 43.5; %GLDAS soil fraction
%err_SandFraction = 14.7; %old
err_SandFraction = 4.89; %std. dev of GLDAS soil fraction InSAR extent
kRoot = 5.5;
err_kRoot = 0.5;
% OrganicMass = 30.0; %old
% err_OrganicMass = 10.0; %old
OrganicMass = 51.9; %Johnson et al. 2011 and NSSI Land Cover for wet sedge
err_OrganicMass = 19.1; %Johnson et al. 2011 std. dev of sandy lowland
% RootDepth = 0.7; %old
% err_RootDepth = 0.2; %old
RootDepth = 0.75; %maximum observed ALT value (found in CALM TLO)
err_RootDepth = 0.14; %stdev of max alt for each year
MaxOrganicMatter = 140.0;
err_MaxOrganicMatter = 10.0;
OrganicSoilPorosity = 0.9;
err_OrganicSoilPorosity = 0.05;
org_depth = 0.18; %! input for alt_from_deform subroutine, but not used in the uncertainty calculation
%SaturationFraction = sat(depth_top,depth_bot,vwc,SandFraction,kRoot,OrganicMass,RootDepth,org_depth,MaxOrganicMatter,OrganicSoilPorosity) % ! 1.0
SaturationFraction=ones(size(deform));      %assume fully saturated for all pixels

dev_ratio = 0.99; %! keep the ratio as a constant

[ABoVE_ALT,reference] = readgeoraster("../../warped_files/NW_drew/averageALTproduct.tif","OutputType","double","CoordinateSystemType","geographic");
squeezeALT = squeeze(ABoVE_ALT(:,:));
rowALT = squeezeALT(:);
avgALT = mean(rowALT,'omitnan'); %avg ALT from ABoVE map in cm

[dummy_alt1,ref_alt,dummy_alt2,dummy_alt3,defo]=alt_from_deform_V3(deform,SaturationFraction, SandFraction, kRoot, OrganicMass, RootDepth,org_depth, MaxOrganicMatter, OrganicSoilPorosity, ABoVE_ALT);
size1=size(deform);
err_alt = zeros(size1(1),size1(2));
nominal_alt=ref_alt;
err_defo = zeros(size1(1),size1(2));

%% calculate ALT uncertainty
tic
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% ! deformation uncertainty
dev_param = dev_ratio*deform;
if dev_param<0.1 
    dev_param=0.1; %!!cautious, don't undersdand why keep deform > 0.1 cm
end
%! calculate new alt
[dummy_alt1, dev_alt, dummy_alt2, dummy_alt3]=alt_from_deform_V3(dev_param, SaturationFraction, SandFraction, kRoot, OrganicMass, RootDepth,org_depth, MaxOrganicMatter, OrganicSoilPorosity, ABoVE_ALT);

%! calculate the gradient
slope=(dev_alt-ref_alt)./(dev_param-deform);

err_deform = RMSE_extrap.'; 
% err_deform_mask = err_deform.*mask_51;
% avg_err_deform_mask = mean(err_deform_mask(:),'omitmissing');
%! calculate the uncertainty
single_err_alt1=slope.*err_deform;
% single_err_alt1_1=slope.*err_deform.*Armask; %not sure about this?
% single_err_alt1_2=slope.*err_deform.*Armask./2;
% single_err_alt1_4=slope.*err_deform.*Armask./4;


err_alt=err_alt+single_err_alt1.^2;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%! satfrac uncertainty
dev_param = dev_ratio*SaturationFraction;
%! calculate new alt
[dummy_alt1, dev_alt, dummy_alt2, dummy_alt3]=alt_from_deform_V3(deform, dev_param, SandFraction, kRoot, OrganicMass, RootDepth,org_depth, MaxOrganicMatter, OrganicSoilPorosity, ABoVE_ALT);
%! calculate the gradient
slope = (dev_alt-ref_alt) ./ (dev_param-SaturationFraction);
	      
%! calculate the uncertainty
single_err_alt2 = slope .* err_SaturationFraction;
err_alt=err_alt+single_err_alt2.^2;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%! sand fraction uncertainty
dev_param = dev_ratio*SandFraction;
%! calculate new alt
[dummy_alt1, dev_alt, dummy_alt2, dummy_alt3]=alt_from_deform_V3(deform, SaturationFraction, dev_param, kRoot, OrganicMass, RootDepth,org_depth, MaxOrganicMatter, OrganicSoilPorosity, ABoVE_ALT);
%! calculate the gradient
slope = (dev_alt-ref_alt) ./ (dev_param-SandFraction);
	      
%! calculate the uncertainty
single_err_alt3 = slope .* err_SandFraction;
err_alt = err_alt + single_err_alt3.^2;


%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%! kroot uncertainty
dev_param = dev_ratio*kRoot;
% calculate new alt
[dummy_alt1, dev_alt, dummy_alt2, dummy_alt3]=alt_from_deform_V3(deform, SaturationFraction, SandFraction, dev_param, OrganicMass, RootDepth,org_depth, MaxOrganicMatter, OrganicSoilPorosity, ABoVE_ALT);
% calculate the gradient
slope = (dev_alt-ref_alt) ./ (dev_param-kRoot);
	      
%! %calculate the uncertainty
single_err_alt4 = slope .* err_kRoot;
err_alt = err_alt + single_err_alt4.^2;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%! total organic mass uncertainty
dev_param = dev_ratio*OrganicMass;
%! calculate new alt
[dummy_alt1, dev_alt, dummy_alt2, dummy_alt3]=alt_from_deform_V3(deform, SaturationFraction, SandFraction, kRoot, dev_param, RootDepth,org_depth, MaxOrganicMatter, OrganicSoilPorosity, ABoVE_ALT);
%! calculate the gradient
slope = (dev_alt-ref_alt) ./ (dev_param-OrganicMass);
	      
%! calculate the uncertainty
single_err_alt5 = slope .* err_OrganicMass;
err_alt = err_alt + single_err_alt5.^2;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%! rooting depth uncertainty
dev_param = dev_ratio*RootDepth;
%! calculate new alt
[dummy_alt1, dev_alt, dummy_alt2, dummy_alt3]=alt_from_deform_V3(deform, SaturationFraction, SandFraction, kRoot, OrganicMass, dev_param,org_depth, MaxOrganicMatter, OrganicSoilPorosity, ABoVE_ALT);
%! calculate the gradient
slope = (dev_alt-ref_alt) ./ (dev_param-RootDepth);
	      
%! calculate the uncertainty
single_err_alt6 = slope .* err_RootDepth;
err_alt = err_alt + single_err_alt6.^2;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%! density organic matter max uncertainty
dev_param = dev_ratio*MaxOrganicMatter;
%! calculate new alt
[dummy_alt1, dev_alt, dummy_alt2, dummy_alt3]=alt_from_deform_V3(deform, SaturationFraction, SandFraction, kRoot, OrganicMass, RootDepth,org_depth, dev_param, OrganicSoilPorosity, ABoVE_ALT);
%! calculate the gradient
slope = (dev_alt-ref_alt) ./ (dev_param-MaxOrganicMatter);
	      
%! calculate the uncertainty
single_err_alt7 = slope .* err_MaxOrganicMatter;
err_alt = err_alt + single_err_alt7.^2;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%! porosity organic matter uncertainty
dev_param = dev_ratio*OrganicSoilPorosity;
%! calculate new alt
[dummy_alt1, dev_alt, dummy_alt2, dummy_alt3]=alt_from_deform_V3(deform, SaturationFraction, SandFraction, kRoot, OrganicMass, RootDepth,org_depth, MaxOrganicMatter, dev_param, ABoVE_ALT);
%! calculate the gradient
slope = (dev_alt-ref_alt) ./ (dev_param-OrganicSoilPorosity);
	      
%! calculate the uncertainty
single_err_alt8 = slope .* err_OrganicSoilPorosity;
err_alt = err_alt + single_err_alt8.^2;


joint_error=sqrt(err_alt);
alt=nominal_alt;
toc

%% uncertainty in ALT -> deform;
tic
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% ! ALT uncertainty
dev_param = dev_ratio*ABoVE_ALT; %cm
if dev_param<0.1 
    dev_param=0.1; %!!cautious, don't undersdand why keep deform > 0.1 cm
end
%! calculate new defo
[dummy_alt1, ref_alt, dummy_alt2, dummy_alt3,dev_defo]=alt_from_deform_V3(dev_param, SaturationFraction, SandFraction, kRoot, OrganicMass, RootDepth,org_depth, MaxOrganicMatter, OrganicSoilPorosity,ABoVE_ALT);

%! calculate the gradient
slope=(dev_defo-defo)./(dev_param-ABoVE_ALT);
[err_ALT_cm,reference] = readgeoraster("../../warped_files/NW_drew/combinedALTproductUncertainty_v2.tif"); %combination of Yi and Whitcomb Uncertainties


%! calculate the uncertainty
single_err_defo1=slope.*err_ALT_cm;
% single_err_alt1_1=slope.*err_deform.*Armask; %not sure about this?
% single_err_alt1_2=slope.*err_deform.*Armask./2;
% single_err_alt1_4=slope.*err_deform.*Armask./4;


err_defo=err_defo+single_err_defo1.^2;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%! satfrac uncertainty
dev_param = dev_ratio*SaturationFraction;
%! calculate new defo
[dummy_alt1, ref_alt, dummy_alt2, dummy_alt3,dev_defo]=alt_from_deform_V3(deform, dev_param, SandFraction, kRoot, OrganicMass, RootDepth,org_depth, MaxOrganicMatter, OrganicSoilPorosity, ABoVE_ALT);
%! calculate the gradient
slope = (dev_defo-defo) ./ (dev_param-SaturationFraction);
	      
%! calculate the uncertainty
single_err_defo2 = slope .* err_SaturationFraction;
err_defo=err_defo+single_err_defo2.^2;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%! sand fraction uncertainty
dev_param = dev_ratio*SandFraction;
%! calculate new defo
[dummy_alt1, ref_alt, dummy_alt2, dummy_alt3,dev_defo]=alt_from_deform_V3(deform, SaturationFraction, dev_param, kRoot, OrganicMass, RootDepth,org_depth, MaxOrganicMatter, OrganicSoilPorosity, ABoVE_ALT);
%! calculate the gradient
slope = (dev_defo-defo) ./ (dev_param-SandFraction);
	      
%! calculate the uncertainty
single_err_defo3 = slope .* err_SandFraction;
err_defo = err_defo + single_err_defo3.^2;


%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%! kroot uncertainty
dev_param = dev_ratio*kRoot;
% calculate new defo
[dummy_alt1, ref_alt, dummy_alt2, dummy_alt3,dev_defo]=alt_from_deform_V3(deform, SaturationFraction, SandFraction, dev_param, OrganicMass, RootDepth,org_depth, MaxOrganicMatter, OrganicSoilPorosity, ABoVE_ALT);
% calculate the gradient
slope = (dev_defo-defo) ./ (dev_param-kRoot);
	      
%! %calculate the uncertainty
single_err_defo4 = slope .* err_kRoot;
err_defo = err_defo + single_err_defo4.^2;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%! total organic mass uncertainty
dev_param = dev_ratio*OrganicMass;
%! calculate new defo
[dummy_alt1, ref_alt, dummy_alt2, dummy_alt3,dev_defo]=alt_from_deform_V3(deform, SaturationFraction, SandFraction, kRoot, dev_param, RootDepth,org_depth, MaxOrganicMatter, OrganicSoilPorosity, ABoVE_ALT);
%! calculate the gradient
slope = (dev_defo-defo) ./ (dev_param-OrganicMass);
	      
%! calculate the uncertainty
single_err_defo5 = slope .* err_OrganicMass;
err_defo = err_defo + single_err_defo5.^2;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%! rooting depth uncertainty
dev_param = dev_ratio*RootDepth;
%! calculate new defo
[dummy_alt1, ref_alt, dummy_alt2, dummy_alt3,dev_defo]=alt_from_deform_V3(deform, SaturationFraction, SandFraction, kRoot, OrganicMass, dev_param,org_depth, MaxOrganicMatter, OrganicSoilPorosity, ABoVE_ALT);
%! calculate the gradient
slope = (dev_defo-defo) ./ (dev_param-RootDepth);
	      
%! calculate the uncertainty
single_err_defo6 = slope .* err_RootDepth;
err_defo = err_defo + single_err_defo6.^2;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%! density organic matter max uncertainty
dev_param = dev_ratio*MaxOrganicMatter;
%! calculate new defo
[dummy_alt1, ref_alt, dummy_alt2, dummy_alt3,dev_defo]=alt_from_deform_V3(deform, SaturationFraction, SandFraction, kRoot, OrganicMass, RootDepth,org_depth, dev_param, OrganicSoilPorosity, ABoVE_ALT);
%! calculate the gradient
slope = (dev_defo-defo) ./ (dev_param-MaxOrganicMatter);
	      
%! calculate the uncertainty
single_err_defo7 = slope .* err_MaxOrganicMatter;
err_defo = err_defo + single_err_defo7.^2;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%! porosity organic matter uncertainty
dev_param = dev_ratio*OrganicSoilPorosity;
%! calculate new defo
[dummy_alt1, ref_alt, dummy_alt2, dummy_alt3,dev_defo]=alt_from_deform_V3(deform, SaturationFraction, SandFraction, kRoot, OrganicMass, RootDepth,org_depth, MaxOrganicMatter, dev_param, ABoVE_ALT);
%! calculate the gradient
slope = (dev_defo-defo) ./ (dev_param-OrganicSoilPorosity);
	      
%! calculate the uncertainty
single_err_defo8 = slope .* err_OrganicSoilPorosity;
err_defo = err_defo + single_err_defo8.^2;



joint_error_defo=sqrt(err_defo); %add this in quadrature with deform uncertainty; then multiply or divide w/ error in ReSALT ALT for EIC%
alt=nominal_alt;
toc
%% calculate %EIC uncert.
defo_FALT = defo; %defo, as expected from organic layer model
Delta_deform = deform-defo_FALT; %functionally is thickness of GIC that is seasonally melted
Delta_deform( Delta_deform < 0) = 0;

Percent_GIC = (Delta_deform./ABoVE_ALT).*100; %percent volume of GIC in the active layer

error_EIC_pt1 = sqrt((err_deform.^2) + (joint_error_defo.^2));
error_EICpt2 = (Percent_GIC./100).*sqrt(((error_EIC_pt1./Delta_deform).^2)+((err_ALT_cm./ABoVE_ALT).^2));
error_EIC = (error_EICpt2.*100);

%% FIGURES BELOW
%% E & ALT fig
ref_alt_nan = ref_alt;
ref_alt_nan(ref_alt==0) = NaN;
figure
tiledlayout(2,2)
nexttile
h =imagesc(E');
set(gca,'color',0.5*[1 1 1])
set(h, 'AlphaData', ~isnan(E'))
title([num2str(year),' Seasonal Subsidence'])
c = colorbar;
xlabel(c, 'E [cm]')

stats_array = squeeze(E);
row_stats_array = (stats_array(:));
percent5 = prctile(row_stats_array,5); %colorbar lims
percent95 = prctile(row_stats_array,95);
clim([percent5 percent95])

nexttile
RMSE(RMSE==0)=NaN;
h =imagesc(RMSE');
set(gca,'color',0.5*[1 1 1])
set(h, 'AlphaData', ~isnan(RMSE'))
title([num2str(year),' Uncertainty in Seasonal Subsidence'])
c = colorbar;
xlabel(c, 'E [cm]')

stats_array = squeeze(RMSE);
row_stats_array = (stats_array(:));
percent5 = prctile(row_stats_array,5); %colorbar lims
percent95 = prctile(row_stats_array,95);
clim([percent5 percent95])



nexttile
h =imagesc(ref_alt_nan);
set(gca,'color',0.5*[1 1 1])
set(h, 'AlphaData', ~isnan(ref_alt_nan))
title([num2str(year),' Active Layer Thickness'])
c = colorbar;
xlabel(c, 'ALT [cm]')

stats_array = squeeze(ref_alt_nan);
row_stats_array = (stats_array(:));
percent5 = prctile(row_stats_array,5); %colorbar lims
percent95 = prctile(row_stats_array,95);
clim([percent5 percent95])

nexttile
h = imagesc(joint_error);
set(gca,'color',0.5*[1 1 1])
set(h, 'AlphaData', ~isnan(joint_error))
title([num2str(year),' Active Layer Thickness Uncertainty'])
c = colorbar;
xlabel(c, 'ALT [cm]')

stats_array = squeeze(joint_error);
row_stats_array = (stats_array(:));
percent5 = prctile(row_stats_array,5); %colorbar lims
percent95 = prctile(row_stats_array,95);
clim([percent5 percent95])
%% EIC plotting
deform_mask = deform;
defo_FALT_mask = defo_FALT;
Delta_deform_mask = Delta_deform;
Percent_GIC_mask = Percent_GIC;

imAlpha=ones(size(defo_FALT_mask));
imAlpha(isnan(defo_FALT_mask))=0; %maybe this is a workaround? all -1 pixels are masked out?

figure
tiledlayout(2,2)
nexttile
h=imagesc(deform_mask);
set(gca,'color',0.5*[1 1 1]);
set(h, 'AlphaData', ~isnan(deform_mask))
title('InSAR-derived Seasonal Subsidence')
c = colorbar;
xlabel(c, 'Subsidence [cm]')
set(gca,'XTick',[378 756 1134 1512 1890])      %~pixel values corresponding to long.
set(gca,'XTickLabel',[-153.79 -153.37 -152.95 -152.53 -152.11])
set(gca,'YTick',[254 508 762 1016])      %~pixel values corresponding to lat.
set(gca,'YTickLabel',[70.81 70.67 70.53 70.38])
xlabel('Longitude')
ylabel('Latitude')
stats_array = squeeze(deform);
row_stats_array = (stats_array(:));
percent5 = prctile(row_stats_array,5); %colorbar lims
percent95 = prctile(row_stats_array,95);
clim([percent5 percent95])

nexttile
defo_FALT_mask(defo_FALT_mask==0) = NaN;
h=imagesc(defo_FALT_mask);
set(gca,'color',0.5*[1 1 1]);
set(h, 'AlphaData', ~isnan(defo_FALT_mask))
title(['Field ALT-inverted Seasonal Subsidence'])
c = colorbar;
xlabel(c, 'Subsidence [cm]')
set(gca,'XTick',[378 756 1134 1512 1890])      %~pixel values corresponding to long.
set(gca,'XTickLabel',[-153.79 -153.37 -152.95 -152.53 -152.11])
set(gca,'YTick',[254 508 762 1016])      %~pixel values corresponding to lat.
set(gca,'YTickLabel',[70.81 70.67 70.53 70.38])
xlabel('Longitude')
ylabel('Latitude')

stats_array = squeeze(defo_FALT_mask);
row_stats_array = (stats_array(:));
percent5 = prctile(row_stats_array,5); %colorbar lims
percent95 = prctile(row_stats_array,95);
clim([percent5 percent95])

nexttile
h=imagesc(Delta_deform_mask);
set(gca,'color',0.5*[1 1 1]);
set(h, 'AlphaData', ~isnan(Delta_deform_mask))
title('Ground Ice Content Thickness')
c = colorbar;
xlabel(c, 'Thickness [cm]')
set(gca,'XTick',[378 756 1134 1512 1890])      %~pixel values corresponding to long.
set(gca,'XTickLabel',[-153.79 -153.37 -152.95 -152.53 -152.11])
set(gca,'YTick',[254 508 762 1016])      %~pixel values corresponding to lat.
set(gca,'YTickLabel',[70.81 70.67 70.53 70.38])
xlabel('Longitude')
ylabel('Latitude')

stats_array = squeeze(Delta_deform_mask);
row_stats_array = (stats_array(:));
percent5 = prctile(row_stats_array,5); %colorbar lims
percent95 = prctile(row_stats_array,95);
clim([percent5 percent95])

nexttile
h=imagesc(Percent_GIC_mask);
set(gca,'color',0.5*[1 1 1]);
set(h, 'AlphaData', ~isnan(Percent_GIC_mask))
title('Percent Ground Ice Content in the Active Layer')
c = colorbar;
xlabel(c, 'Percent')
set(gca,'XTick',[378 756 1134 1512 1890])      %~pixel values corresponding to long.
set(gca,'XTickLabel',[-153.79 -153.37 -152.95 -152.53 -152.11])
set(gca,'YTick',[254 508 762 1016])      %~pixel values corresponding to lat.
set(gca,'YTickLabel',[70.81 70.67 70.53 70.38])
xlabel('Longitude')
ylabel('Latitude')

stats_array = squeeze(Percent_GIC_mask);
row_stats_array = (stats_array(:));
percent5 = prctile(row_stats_array,5); %colorbar lims
percent95 = prctile(row_stats_array,95);
clim([percent5 percent95])

%% EIC 3x2
figure
tiledlayout(3,2)
nexttile
h = imagesc(defo_FALT);
set(gca,'color',0.5*[1 1 1]);
set(h, 'AlphaData', ~isnan(defo_FALT))
title('Independent Forward-modeled Seasonal Subsidence')
c = colorbar;
xlabel(c, 'E [cm]')

stats_array = squeeze(defo_FALT);
row_stats_array = (stats_array(:));
percent5 = prctile(row_stats_array,5); %colorbar lims
percent95 = prctile(row_stats_array,95);
clim([percent5 percent95])

nexttile
h = imagesc(joint_error_defo);
set(gca,'color',0.5*[1 1 1]);
set(h, 'AlphaData', ~isnan(joint_error_defo))
title('Independent Forward-modeled Seasonal Subsidence Uncertainty')
c = colorbar;
xlabel(c, 'E [cm]')

stats_array = squeeze(joint_error_defo);
row_stats_array = (stats_array(:));
percent5 = prctile(row_stats_array,5); %colorbar lims
percent95 = prctile(row_stats_array,95);
clim([percent5 percent95])

nexttile
h = imagesc(Delta_deform);
set(gca,'color',0.5*[1 1 1]);
set(h, 'AlphaData', ~isnan(Delta_deform))
title('Excess Ice Content Thickness')
c = colorbar;
xlabel(c, 'EIC [cm]')

stats_array = squeeze(Delta_deform);
row_stats_array = (stats_array(:));
percent5 = prctile(row_stats_array,5); %colorbar lims
percent95 = prctile(row_stats_array,95);
clim([percent5 percent95])

nexttile
h = imagesc(error_EIC_pt1);
set(gca,'color',0.5*[1 1 1]);
set(h, 'AlphaData', ~isnan(error_EIC_pt1))
title('Excess Ice Content Thickness Uncertainty')
c = colorbar;
xlabel(c, 'E [cm]')

stats_array = squeeze(error_EIC_pt1);
row_stats_array = (stats_array(:));
percent5 = prctile(row_stats_array,5); %colorbar lims
percent95 = prctile(row_stats_array,95);
clim([percent5 percent95])

nexttile
h = imagesc(Percent_GIC);
set(gca,'color',0.5*[1 1 1]);
set(h, 'AlphaData', ~isnan(Percent_GIC))
title('Percent Excess Ice Content in the Active Layer')
c = colorbar;
xlabel(c, '[%]')

stats_array = squeeze(Percent_GIC);
row_stats_array = (stats_array(:));
percent5 = prctile(row_stats_array,5); %colorbar lims
percent95 = prctile(row_stats_array,95);
clim([percent5 percent95])

nexttile
h = imagesc(error_EIC);
set(gca,'color',0.5*[1 1 1]);
set(h, 'AlphaData', ~isnan(error_EIC))
title('Percent Excess Ice Content in the Active Layer Uncertainty')
c = colorbar;
xlabel(c, '[%]')

stats_array = squeeze(error_EIC);
row_stats_array = (stats_array(:));
percent5 = prctile(row_stats_array,5); %colorbar lims
percent95 = prctile(row_stats_array,95);
clim([percent5 percent95])

%% ALT cores and probes

ALT_core_string = strcat("../../warped_files/NW_drew/ALT_cores_",num2str(year),".xls");
if isfile(ALT_core_string)
    disp('Core file exits. Continue to output ALT 1-to-1.')

    ALT_cores = readtable(ALT_core_string,'VariableNamingRule','preserve'); %read in table of pixels of interest
    longitude = ALT_cores.Longitude;
    latitude = ALT_cores.Latitude;
    entries_ALT = length(longitude);
    %names_ALT = ALT_cores.Borehole;
    
    insitu_ALT = ALT_cores.("ALT (cm)");

    insitu_ALT(isnan(insitu_ALT)) = [];

    pixelX_ALT = zeros(entries_ALT,1);
    pixelY_ALT = zeros(entries_ALT,1);

    for i = 1:entries_ALT
        pixelX_ALT(i,1) = ceil(abs((xfirst+(-longitude(i,1)))/xstep));
        pixelY_ALT(i,1) = ceil(abs((yfirst-latitude(i,1))/ystep));
    end
    
    %Look for cores within the same pixel
    ALT_radar_coords = [pixelX_ALT,pixelY_ALT];
    ALT_radar_coords = strcat(num2str(ALT_radar_coords(:,1)),",",num2str(ALT_radar_coords(:,2)));

    [uniqueStrings_ALT, ~, indices_ALT] = unique(ALT_radar_coords);

    counts_ALT = histcounts(indices_ALT, 1:length(uniqueStrings_ALT)+1);
    duplicateStrings_ALT = uniqueStrings_ALT(counts_ALT>1);

    if size(duplicateStrings_ALT)>0
        disp('Multiple ALT cores in the same pixel. Spatially average.')
        insitu_tbl = table(insitu_ALT,indices_ALT);
        insitu_tbl = sortrows(insitu_tbl,"indices_ALT","ascend");
        average_tbl = splitapply(@mean,insitu_tbl(:,"insitu_ALT"),findgroups(insitu_tbl(:,"indices_ALT")));
        stdev_tbl = splitapply(@std,insitu_tbl(:,"insitu_ALT"),findgroups(insitu_tbl(:,"indices_ALT")));
        stdev_tbl = sqrt((stdev_tbl.^2)+(5^2)); %spatial avg std + measurement uncertainty in quadrature

        splitstring = split(uniqueStrings_ALT,",");
        tblnames = ["ALT_cm","PixelX","PixelY","Uncert_cm"];
        combined_core_table = table(average_tbl,str2double(splitstring(:,1)),str2double(splitstring(:,2)),stdev_tbl,'VariableNames',tblnames);
    
        PixelX_ALT_avg =combined_core_table(:,2);
        PixelY_ALT_avg = combined_core_table(:,3);

        ReSALT_ALT = zeros(size(combined_core_table,1),1);
        Error_ALT = zeros(size(combined_core_table,1),1);

        for i = 1:size(combined_core_table,1)
            ReSALT_ALT(i,1) = ref_alt_nan(PixelY_ALT_avg(i,1),PixelX_ALT_avg(i,1));
            Error_ALT(i,1)= joint_error(PixelY_ALT_avg(i,1),PixelX_ALT_avg(i,1));
        end    
    
    else
        disp('No ALT cores in the same pixel, so do not spatially average cores.')
        stdev_tbl = zeros(entries_ALT,1);
        stdev_tbl(:) = sqrt((5^2)); %just measurement uncertainty; should be for some years w/ only cores
        tblnames = ["ALT_cm","PixelX","PixelY","Uncert_cm"];
        combined_core_table = table(insitu_ALT,pixelX_ALT,pixelY_ALT,stdev_tbl,'VariableNames',tblnames);
    
        ReSALT_ALT = zeros(size(combined_core_table,1),1);
        Error_ALT = zeros(size(combined_core_table,1),1);

        for i = 1:size(combined_core_table,1)
            ReSALT_ALT(i,1) = ref_alt_nan(pixelY_ALT(i,1),pixelX_ALT(i,1));
            Error_ALT(i,1)= joint_error(pixelY_ALT(i,1),pixelX_ALT(i,1));
        end 
    end

    figure
    scatter(combined_core_table.ALT_cm,ReSALT_ALT,'filled')
    hold on

    xlabel('In situ ALT [cm]')
    ylabel('InSAR-derived ALT [cm]')

    lim = max(max(combined_core_table.ALT_cm),max(ReSALT_ALT));
    lim = round((lim + max(Error_ALT) + 10)/5)*5;
    xlim([0, lim])
    ylim([0, lim])
    plot([0 xlim],[0 ylim],'--r')

    f = errorbar(combined_core_table.ALT_cm,ReSALT_ALT,Error_ALT,'vertical','Linestyle','none');
    f.Color = 'k';
    e = errorbar(combined_core_table.ALT_cm,ReSALT_ALT,combined_core_table.Uncert_cm,'horizontal','linestyle','none');
    e.Color = 'k';

    legend('ALT','Idealized 1-to-1 fit')
    ax = gca;
    ax.FontSize = 20;

else
    disp("There is no ALT core data for this year.")
end

ALT_probe_string = strcat("ALT_probes_",num2str(year),".xls");
if isfile(ALT_probe_string)
    disp('Probe file exits. Continue to output ALT 1-to-1.')

    ALT_probes = readtable(ALT_probe_string,'VariableNamingRule','preserve'); %read in table of pixels of interest
    longitude = ALT_probes.Longitude;
    latitude = ALT_probes.Latitude;
    entries_ALT = length(longitude);
    %names_ALT = ALT_cores.Borehole;
    
    insitu_ALT = ALT_probes.("ALT (cm)");

    insitu_ALT(isnan(insitu_ALT)) = [];

    pixelX_ALT = zeros(entries_ALT,1);
    pixelY_ALT = zeros(entries_ALT,1);

    for i = 1:entries_ALT
        pixelX_ALT(i,1) = ceil(abs((xfirst+(-longitude(i,1)))/xstep));
        pixelY_ALT(i,1) = ceil(abs((yfirst-latitude(i,1))/ystep));
    end
    
    %Look for cores within the same pixel
    ALT_radar_coords = [pixelX_ALT,pixelY_ALT];
    ALT_radar_coords = strcat(num2str(ALT_radar_coords(:,1)),",",num2str(ALT_radar_coords(:,2)));

    [uniqueStrings_ALT, ~, indices_ALT] = unique(ALT_radar_coords);

    counts_ALT = histcounts(indices_ALT, 1:length(uniqueStrings_ALT)+1);
    duplicateStrings_ALT = uniqueStrings_ALT(counts_ALT>1);

    if size(duplicateStrings_ALT)>0
        disp('Multiple ALT probes in the same pixel. Spatially average.')
        insitu_tbl = table(insitu_ALT,indices_ALT);
        insitu_tbl = sortrows(insitu_tbl,"indices_ALT","ascend");
        average_probe_tbl = splitapply(@mean,insitu_tbl(:,"insitu_ALT"),findgroups(insitu_tbl(:,"indices_ALT")));
        stdev_tbl = splitapply(@std,insitu_tbl(:,"insitu_ALT"),findgroups(insitu_tbl(:,"indices_ALT")));
        stdev_tbl = sqrt((stdev_tbl.^2)+(3^2)); %spatial avg std + measurement uncertainty in quadrature

        splitstring = split(uniqueStrings_ALT,",");
        tblnames = ["ALT_cm","PixelX","PixelY","Uncert_cm"];
        combined_probe_table = table(average_probe_tbl,str2double(splitstring(:,1)),str2double(splitstring(:,2)),stdev_tbl,'VariableNames',tblnames);
    end

    PixelX_ALT_avg =combined_probe_table(:,2);
    PixelY_ALT_avg = combined_probe_table(:,3);

    ReSALT_ALT = zeros(size(combined_probe_table,1),1);
    Error_ALT = zeros(size(combined_probe_table,1),1);

    for i = 1:size(combined_probe_table,1)
        ReSALT_ALT(i,1) = ref_alt_nan(PixelY_ALT_avg(i,1),PixelX_ALT_avg(i,1));
        Error_ALT(i,1)= joint_error(PixelY_ALT_avg(i,1),PixelX_ALT_avg(i,1));
    end

    figure
    scatter(combined_probe_table.ALT_cm,ReSALT_ALT,'filled')
    hold on

    xlabel('In situ ALT [cm]')
    ylabel('InSAR-derived ALT [cm]')

    lim = max(max(combined_probe_table.ATL_cm),max(ReSALT_ALT));
    lim = round((lim + max(Error_ALT) + 10)/5)*5;
    xlim([0, lim])
    ylim([0, lim])
    plot([0 xlim],[0 ylim],'--r')

    f = errorbar(combined_probe_table.ALT_cm,ReSALT_ALT,Error_ALT,'vertical','Linestyle','none');
    f.Color = 'k';
    e = errorbar(combined_probe_table.ALT_cm,ReSALT_ALT,combined_probe_table.Uncert_cm,'horizontal','linestyle','none');
    e.Color = 'k';

    legend('ALT','Idealized 1-to-1 fit')
    ax = gca;
    ax.FontSize = 20;

else
    disp("There is no probe data for this year.")
end

if isfile(ALT_core_string) %check if both core and probe data exists; if so, combine
    disp("Core data exits.")
    if isfile(ALT_probe_string)
        disp("Probe data also exists. Combine tables.")
        full_ALT_table = [combined_core_table;combined_probe_table];
        
        PixelX_ALT_com = full_ALT_table(:,2);
        PixelY_ALT_com = full_ALT_table(:,3);
        
        ReSALT_ALT = zeros(size(full_ALT_table,1),1);
        Error_ALT = zeros(size(full_ALT_table,1),1);

        for i = 1:size(full_ALT_table,1)
            ReSALT_ALT(i,1) = ref_alt_nan(PixelY_ALT_com(i,1),PixelX_ALT_com(i,1));
            Error_ALT(i,1)= joint_error(PixelY_ALT_com(i,1),PixelX_ALT_com(i,1));
        end
        
        figure
        scatter(full_ALT_table.ALT_cm,ReSALT_ALT,'filled')
        hold on

        xlabel('In situ ALT [cm]')
        ylabel('InSAR-derived ALT [cm]')

        lim = max(max(full_ALT_table.ALT_cm),max(ReSALT_ALT));
        lim = round((lim + max(Error_ALT) + 10)/5)*5;
        xlim([0, lim])
        ylim([0, lim])
        plot([0 xlim],[0 ylim],'--r')

        f = errorbar(full_ALT_table.ALT_cm,ReSALT_ALT,Error_ALT,'vertical','Linestyle','none');
        f.Color = 'k';
        e = errorbar(full_ALT_table.ALT_cm,ReSALT_ALT,full_ALT_table.Uncert_cm,'horizontal','linestyle','none');
        e.Color = 'k';

        legend('ALT','Idealized 1-to-1 fit')
        ax = gca;
        ax.FontSize = 20;

    else
        disp("Probe data does not exist. Only core data, so no need for combination of tables.")
    end
    
else
    disp("Core data does not exist, so no need for combination of tables.")
end
disp("In situ ALT plotting section complete!")

%% EIC cores
EIC_core_string = strcat("../../warped_files/NW_drew/EIC_cores_",num2str(year),".xls");
if isfile(EIC_core_string)
    disp('EIC file exits. Continue to output EIC 1-to-1.')

    EIC_cores = readtable(EIC_core_string,'VariableNamingRule','preserve'); %read in table of pixels of interest
    longitude = EIC_cores.Longitude;
    latitude = EIC_cores.Latitude;
    entries = length(longitude);
    %names = EIC_cores.Borehole;
    
    insitu_EIC = EIC_cores.EIC____;

    pixelX_EIC = zeros(entries,1);
    pixelY_EIC = zeros(entries,1);

    for i = 1:entries
        pixelX_EIC(i,1) = ceil(abs((xfirst+(-longitude(i,1)))/xstep));
        pixelY_EIC(i,1) = ceil(abs((yfirst-latitude(i,1))/ystep));
    end

    %Look for cores within the same pixel
    EIC_radar_coords = [pixelX_EIC,pixelY_EIC];
    EIC_radar_coords = strcat(num2str(EIC_radar_coords(:,1)),",",num2str(EIC_radar_coords(:,2)));

    [uniqueStrings, ~, indices_EIC] = unique(EIC_radar_coords);

    counts = histcounts(indices_EIC, 1:length(uniqueStrings)+1);
    duplicateStrings = uniqueStrings(counts>1);

    if size(duplicateStrings)>0
        disp('Multiple EIC cores in the same pixel. Spatially average.')

        insitu_tbl = table(insitu_EIC,indices_EIC);
        insitu_tbl = sortrows(insitu_tbl,"indices_EIC","ascend");
        average_EIC_tbl = splitapply(@mean,insitu_tbl(:,"insitu_EIC"),findgroups(insitu_tbl(:,"indices_EIC")));
        stdev_tbl = splitapply(@std,insitu_tbl(:,"insitu_EIC"),findgroups(insitu_tbl(:,"indices_EIC")));
        stdev_tbl = sqrt((stdev_tbl.^2)+(std(insitu_EIC).^2)); %spatial avg std + measurement uncertainty in quadrature

        splitstring = split(uniqueStrings,",");
        tblnames = ["EIC_%","PixelX","PixelY","Uncert_%"];
        combined_EIC_core_table = table(average_EIC_tbl,str2double(splitstring(:,1)),str2double(splitstring(:,2)),stdev_tbl,'VariableNames',tblnames);
    else
        disp('No EIC cores in the same pixel, so do not spatially average cores.')
        stdev_tbl = zeros(entries,1);
        stdev_tbl(:) = sqrt(std(insitu_EIC).^2); %just measurement uncertainty; should be for some years w/ only cores
        tblnames = ["EIC_%","PixelX","PixelY","Uncert_%"];
        combined_EIC_core_table = table(insitu_EIC,pixelX_EIC,pixelY_EIC,stdev_tbl,'VariableNames',tblnames);
    end

    PixelX_EIC = combined_EIC_core_table.PixelX;
    PixelY_EIC = combined_EIC_core_table.PixelY;
    
    ReSALT_EIC = zeros(entries,1);
    Error_EIC = zeros(entries,1);

    for i = 1:entries
        ReSALT_EIC(i,1) = Percent_GIC(PixelY_EIC(i,1),PixelX_EIC(i,1));
        Error_EIC(i,1) = error_EIC(PixelY_EIC(i,1),PixelX_EIC(i,1));
    end
    figure
    scatter(combined_EIC_core_table.("EIC_%"),ReSALT_EIC,'filled')
    hold on
    xlabel('In situ EIC [%]')
    ylabel('InSAR-derived EIC [%]')

    lim = max(max(insitu_EIC),max(ReSALT_EIC));
    xlim([0,round((lim+10)/5)*5])
    ylim([0,round((lim+10)/5)*5])
    plot(xlim,ylim,'--r')

    f = errorbar(combined_EIC_core_table.("EIC_%"),ReSALT_EIC,Error_EIC,'vertical','Linestyle','none');
    f.Color = 'k';
    e = errorbar(combined_EIC_core_table.("EIC_%"),ReSALT_EIC,combined_EIC_core_table.("Uncert_%"),'horizontal','linestyle','none');
    e.Color = 'k';

    legend('EIC','Idealized 1-to-1 fit')
    ax = gca;
    ax.FontSize = 20;
else
    disp("There is no EIC core data for this year.")
end
disp("In situ EIC plotting section complete!")
%% output geotiffs
outputvars_tiff = {'E','E_uncertainty','ALT','ALT_uncertainty','E_ABoVE','E_ABoVE_uncertainty',...
    'EIC_thickness','EIC_thickness_uncertainty','EIC_percent','EIC_percent_uncertainty','RMSE_extrap', 'Coherence'};
outputvars_tiff = cellfun(@(x) [strcat(num2str(year),'_'), x],outputvars_tiff,'UniformOutput',false);
testing_var = cat(3, deform, RMSE', ref_alt_nan, joint_error, defo_FALT, joint_error_defo, Delta_deform,...
    error_EIC_pt1, Percent_GIC, error_EIC, RMSE_extrap', coh');

for i = 1:length(outputvars_tiff)
    tempvar = testing_var(:,:,i);
    tempvar(isnan(tempvar))=-9999;
    geotiffwrite(outputvars_tiff{i},tempvar,ff)
end