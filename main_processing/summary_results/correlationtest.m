E = imread("../mintpy2019_v8_NW_drew_water/ReSALT/2019_E.tif");
EIC = imread("../mintpy2019_v8_NW_drew_water/ReSALT/2019_EIC_percent.tif");
DEM = imread("warped_elevation.dem.wgs84.tif");
slope = imread("slope.tif");

E(E==-9999) = NaN;
EIC(EIC==-9999) = NaN;
DEM(slope==-9999) = NaN;
slope(slope==-9999) = NaN;

E(isnan(E)) = 0;
EIC(isnan(EIC)) = 0;
DEM(isnan(DEM)) = 0;
slope(isnan(slope)) = 0;

figure
tiledlayout(2,2)
nexttile
imagesc(E)
title('E')
colorbar

nexttile
imagesc(EIC)
title('EIC')
colorbar

nexttile
imagesc(DEM)
title('DEM')
colorbar

nexttile
imagesc(slope)
title('Slope')
colorbar

E_DEM_xcorr = normxcorr2(E,DEM);
E_slope_xcorr = normxcorr2(E,slope);
EIC_DEM_xcorr = normxcorr2(EIC,DEM);
EIC_slope_xcorr = normxcorr2(EIC,slope);

% padY = floor(size(E, 1));
% padX = floor(size(E, 2));
% 
% E_DEM_xcorr = E_DEM_xcorr(padY:end-padY+1, padX:end-padX+1);
% E_slope_xcorr = E_slope_xcorr(padY:end-padY+1, padX:end-padX+1);
% EIC_DEM_xcorr = EIC_DEM_xcorr(padY:end-padY+1, padX:end-padX+1);
% EIC_slope_xcorr = EIC_slope_xcorr(padY:end-padY+1, padX:end-padX+1);

figure
tiledlayout(2,2)
nexttile
imagesc(E_DEM_xcorr)
colorbar

nexttile
imagesc(E_slope_xcorr)
colorbar

nexttile
imagesc(EIC_DEM_xcorr)
colorbar

nexttile
imagesc(EIC_slope_xcorr)
colorbar

%%

E = imread("../mintpy2019_v8_NW_drew_water/ReSALT/2019_E.tif");
EIC = imread("../mintpy2019_v8_NW_drew_water/ReSALT/2019_EIC_percent.tif");
DEM = single(imread("warped_elevation.dem.wgs84.tif"));
slope = imread("slope.tif");

E(E==-9999) = NaN;
EIC(EIC==-9999) = NaN;
DEM(slope==-9999) = NaN;
slope(slope==-9999) = NaN;

squeezing = squeeze(E);
rowsqueezeE = squeezing(:);

squeezing = squeeze(EIC);
rowsqueezeEIC = squeezing(:);

squeezing = squeeze(DEM);
rowsqueezeDEM = squeezing(:);

squeezing = squeeze(slope);
rowsqueezeslope = squeezing(:);

fit1 = fitlm(rowsqueezeDEM,rowsqueezeE,'linear');
fit2 = fitlm(rowsqueezeslope,rowsqueezeE,'linear');
fit3 = fitlm(rowsqueezeDEM,rowsqueezeEIC,'linear');
fit4 = fitlm(rowsqueezeslope,rowsqueezeEIC,'linear');

[r1, p1] = corr(rowsqueezeDEM,rowsqueezeE,'Rows','complete','Type','Spearman');
[r2, p2] = corr(rowsqueezeslope,rowsqueezeE,'Rows','complete','Type','Spearman');
[r3, p3] = corr(rowsqueezeDEM,rowsqueezeEIC,'Rows','complete','Type','Spearman');
[r4, p4] = corr(rowsqueezeslope,rowsqueezeEIC,'Rows','complete','Type','Spearman');

figure
tiledlayout
nexttile
%p = plot(fit1);
scatter(rowsqueezeDEM,rowsqueezeE,'filled')
xlabel('Elevation (m)')
ylabel('Seasonal Subsidence (cm)')
title('')
ax = gca;
ax.FontSize = 15;
%legend('Location','northeast')
ylim([0, 30])
text(0.5, 28, 'A.)','FontSize',15,'FontWeight','bold')
text(0.5,26,['$\rho$ = ',num2str(r1(1))],'FontSize',15,'Interpreter','latex')
text(0.5,24,['p << ',num2str(p1(1))],'fontsize',15)
% text(0.5,26,['y = ',num2str(fit1.Coefficients.Estimate(2)),'x +',num2str(fit1.Coefficients.Estimate(1))],'FontSize',15)
% text(0.5,24,['$R^{2} =$',num2str(round(fit1.Rsquared.Ordinary,2))],'Interpreter','latex','fontsize',15)

nexttile
%p = plot(fit2);
scatter(rowsqueezeslope,rowsqueezeE,'filled')
xlabel(['Slope (',char(176),')'])
ylabel('Seasonal Subsidence (cm)')
title('')
ax = gca;
ax.FontSize = 15;
legend('off')
ylim([0, 30])
text(0.5, 28, 'B.)','FontSize',15,'FontWeight','bold')
text(0.5,26,['$\rho$ = ',num2str(r2(1))],'FontSize',15,'Interpreter','latex')
text(0.5,24,['p << ',num2str(p2(1))],'fontsize',15)
% text(0.5,26,['y = ',num2str(fit2.Coefficients.Estimate(2)),'x +',num2str(fit2.Coefficients.Estimate(1))],'FontSize',15)
% text(0.5,24,['$R^{2} =$',num2str(round(fit2.Rsquared.Ordinary,2))],'Interpreter','latex','fontsize',15)

nexttile
%p = plot(fit3);
scatter(rowsqueezeDEM,rowsqueezeEIC,'filled')
xlabel('Elevation (m)')
ylabel('%EIC')
title('')
ax = gca;
ax.FontSize = 15;
legend('off')
ylim([0, 140])
text(0.5, 135, 'C.)','FontSize',15,'FontWeight','bold')
text(0.5,125,['$\rho$ = ',num2str(r3(1))],'FontSize',15,'Interpreter','latex')
text(0.5,115,['p << ',num2str(p3(1))],'fontsize',15)
% text(0.5,125,['y = ',num2str(fit3.Coefficients.Estimate(2)),'x +',num2str(fit3.Coefficients.Estimate(1))],'FontSize',15)
% text(0.5,115,['$R^{2} =$',num2str(round(fit3.Rsquared.Ordinary,2))],'Interpreter','latex','fontsize',15)

nexttile
%p = plot(fit4);
scatter(rowsqueezeslope,rowsqueezeEIC,'filled')
xlabel(['Slope (',char(176),')'])
ylabel('%EIC')
title('')
ax = gca;
ax.FontSize = 15;
legend('off')
ylim([0, 140])
text(0.5, 135, 'D.)','FontSize',15,'FontWeight','bold')
text(0.5,125,['$\rho$ = ',num2str(r4(1))],'FontSize',15,'Interpreter','latex')
text(0.5,115,['p << ',num2str(p4(1))],'fontsize',15)
% text(0.5,125,['y = ',num2str(fit4.Coefficients.Estimate(2)),'x +',num2str(fit4.Coefficients.Estimate(1))],'FontSize',15)
% text(0.5,115,['$R^{2} =$',num2str(round(fit4.Rsquared.Ordinary,2))],'Interpreter','latex','fontsize',15)