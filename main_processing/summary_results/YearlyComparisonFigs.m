summary_results = readtable("results_stats.xlsx","Sheet","Mean r2");
ALT_r2=summary_results.ALTMeanR2;
EIC_r2=summary_results.EICMeanR2;
ADDT=summary_results.ADDT;
labels=string(summary_results.Var1);
%%
f = figure;
tiledlayout(1,3)
nexttile
scatter(ALT_r2,EIC_r2,200,'blue','filled','o')
for i = 1:length(labels)
    text(ALT_r2(i)+0.1, EIC_r2(i)+0.1, labels(i),'FontSize',20);
end
text(12, 25, 'A.)','FontSize',20,'FontWeight','bold')
xscale('log')
yscale('log')
xlim([10, 1000])
ylim([0,30])
ylabel('EIC $\chi^{2}$',Interpreter='latex',FontSize=20)
xlabel('ALT $\chi^{2}$',Interpreter='latex',FontSize=20)
ax = gca;
ax.FontSize =20;
set(ax, 'Box', 'off')
legend('off')

nexttile
linfit = fitlm(ADDT,ALT_r2,'linear');
p = plot(linfit);
for i = 1:length(labels)
    text(ADDT(i)+10,ALT_r2(i)+5,labels(i),'FontSize',20)
end
p(1).MarkerSize = 12;
p(1).Marker = 'o';
p(1).MarkerFaceColor = 'b';
legend('Location',[0.3839 0.7108 0.1294, 0.1436]);
xlim([200,700])
ylim([0,510])
xlabel(['ADDT (',char(176),'C days)'],'FontSize',20)
ylabel('ALT $\chi^{2}$',Interpreter='latex',FontSize=20)
title('')
text(210,350,['y = ',num2str(linfit.Coefficients.Estimate(2)),'x +',num2str(linfit.Coefficients.Estimate(1))],'FontSize',20)
text(210,310,['$R^{2} =$',num2str(round(linfit.Rsquared.Ordinary,2))],'Interpreter','latex','fontsize',20)
text(210, 480, 'B.)','FontSize',20,'FontWeight','bold')
ax = gca;
ax.FontSize =20;
set(ax, 'Box', 'off')

nexttile
scatter(ADDT, EIC_r2, 200, 'b', 'filled', 'o');
for i = 1:length(labels)
    text(ADDT(i)+10, EIC_r2(i)+0.5, labels(i), 'FontSize', 20);
end
text(210, 19, 'C.)', 'FontSize', 20, 'FontWeight', 'bold');
xlabel(['ADDT (', char(176), 'C days)'], 'FontSize', 20);
ylabel('EIC $\chi^{2}$', 'Interpreter', 'latex', 'FontSize', 20);
xlim([200, 700])
ylim([0, 20])
ax = gca;
ax.FontSize = 20;
set(ax, 'Box', 'off')
%% ReSALT-ADDT comparisons
ADDT_sheet = readtable("results_stats.xlsx","Sheet","ADDT");
Year = string(ADDT_sheet.Year);
sub = ADDT_sheet.Median_Subsidence;
eic = ADDT_sheet.Median_EIC_percent;
ADDT = ADDT_sheet.ADDT;
precip = ADDT_sheet.Precip;

linfit1 = fitlm(ADDT,sub,'linear');
linfit2 = fitlm(ADDT,eic,'linear');
linfit3 = fitlm(precip,sub,'linear');
linfit4 = fitlm(precip,eic,'linear');

f2 = figure;
tiledlayout(2,2)
nexttile
p2 = plot(linfit1);
for i = 1:length(Year)
    text(ADDT(i)+10, sub(i)+.25, Year(i), 'FontSize',15)
end
p2(1).MarkerSize = 12;
p2(1).Marker = 'o';
p2(1).MarkerFaceColor = 'b';
xlabel(['ADDT (', char(176), 'C days)'], 'FontSize', 15);
ylabel('Median Seasonal Subsidence (cm)')
xlim([250, 700])
ylim([1.5, 6.5])
legend('Location',[0.329140054122519 0.614757901463277 0.126979170905219 0.0972735705525509])
text(260,5.8,['y = ',num2str(linfit1.Coefficients.Estimate(2)),'x +',num2str(linfit1.Coefficients.Estimate(1))],'FontSize',15)
text(260,5.5,['$R^{2} =$',num2str(round(linfit1.Rsquared.Ordinary,2))],'Interpreter','latex','fontsize',15)
text(260, 6.25, 'A.)','FontSize',15,'FontWeight','bold')
title('')
ax = gca;
ax.FontSize = 15;
set(ax,'Box', 'off')

nexttile
p3 = plot(linfit2);
for i = 1:length(Year)
    text(ADDT(i)+10, eic(i)+1, Year(i), 'FontSize',15)
end
p3(1).MarkerSize = 12;
p3(1).Marker = 'o';
p3(1).MarkerFaceColor = 'b';
xlabel(['ADDT (', char(176), 'C days)'], 'FontSize', 15);
ylabel('Median %EIC')
xlim([250, 700])
ylim([0, 14])
text(260,12,['y = ',num2str(linfit2.Coefficients.Estimate(2)),'x +',num2str(linfit2.Coefficients.Estimate(1))],'FontSize',15)
text(260,11,['$R^{2} =$',num2str(round(linfit2.Rsquared.Ordinary,2))],'Interpreter','latex','fontsize',15)
text(260, 13, 'B.)','FontSize',15,'FontWeight','bold')
legend('off')
title('')
ax = gca;
ax.FontSize = 15;
set(ax,'Box', 'off')

nexttile
p4 = plot(linfit3);
for i = 1:length(Year)
    text(precip(i)+.25, sub(i)+.25, Year(i), 'FontSize',15)
end
p4(1).MarkerSize = 12;
p4(1).Marker = 'o';
p4(1).MarkerFaceColor = 'b';
xlabel('Total Annual Precipitation (cm)', 'FontSize', 15);
ylabel('Median Seasonal Subsidence (cm)')
xlim([15, 28])
ylim([0, 8])
text(15.5,7,['y = ',num2str(linfit3.Coefficients.Estimate(2)),'x +',num2str(linfit3.Coefficients.Estimate(1))],'FontSize',15)
text(15.5,6.5,['$R^{2} =$',num2str(round(linfit3.Rsquared.Ordinary,2))],'Interpreter','latex','fontsize',15)
text(15.5, 7.5, 'C.)','FontSize',15,'FontWeight','bold')
legend('off')
title('')
ax = gca;
ax.FontSize = 15;
set(ax,'Box', 'off')

nexttile
p5 = plot(linfit4);
for i = 1:length(Year)
    text(precip(i)+.25, eic(i)+.25, Year(i), 'FontSize',15)
end
p5(1).MarkerSize = 12;
p5(1).Marker = 'o';
p5(1).MarkerFaceColor = 'b';
xlabel('Total Annual Precipitation (cm)', 'FontSize', 15);
ylabel('Median %EIC')
xlim([15, 28])
ylim([0, 20])
text(15.5,17.5,['y = ',num2str(linfit4.Coefficients.Estimate(2)),'x +',num2str(linfit4.Coefficients.Estimate(1))],'FontSize',15)
text(15.5,16,['$R^{2} =$',num2str(round(linfit4.Rsquared.Ordinary,2))],'Interpreter','latex','fontsize',15)
text(15.5, 19, 'D.)','FontSize',15,'FontWeight','bold')
legend('off')
title('')
ax = gca;
ax.FontSize = 15;
set(ax,'Box', 'off')
%% age-subsidence plots
summary_age = readtable("results_stats.xlsx","Sheet","TL_Age");

age = summary_age.Age_uncal;
sub2019 = summary_age.Med_Subsidence_2019_cm;
sub2021 = summary_age.Med_Subsidence_2021_cm;
sub2019_uncert = summary_age.Mean_E_Uncert_2019;
sub2021_uncert = summary_age.Mean_E_Uncert_2021;

coastal_age = age([1:2,4:5,7],:);
coastal_sub2019 = sub2019([1:2,4:5,7],:);
coastal_sub2021 = sub2021([1:2,4:5,7],:);

coastal_sub2019uncert = sub2019_uncert([1:2,4:5,7],:);
coastal_sub2021uncert = sub2021_uncert([1:2,4:5,7],:);
coastal_names = summary_age.Sample_basin_ID([1:2,4:5,7],:);

inland_age = age([3,6,8],:);
inland_sub2019 = sub2019([3,6,8],:);
inland_sub2021 = sub2021([3,6,8],:);
inland_sub2019uncert = sub2019_uncert([3,6,8],:);
inland_sub2021uncert = sub2021_uncert([3,6,8],:);
inland_names = summary_age.Sample_basin_ID([3,6,8],:);

fitlincoastal2019 = fitlm(coastal_age,coastal_sub2019,'linear');
fitlininland2019 = fitlm(inland_age,inland_sub2019,"linear");
fitlincoastal2021 = fitlm(coastal_age,coastal_sub2021,'linear');
fitlininland2021 = fitlm(inland_age,inland_sub2021,"linear");

figure
f3 = figure;
tiledlayout(2,2)
nexttile
yline(6.7,'-k','14C DTLB Average = 6.7','Color','k','FontSize',15)
hold on
yline(5.3,'-r','Study Area-Wide Average = 5.3','Color','r','FontSize',15)
pi9 = plot(fitlininland2019);
% for i = 1:3
%     text(inland_age(i)+10, inland_sub2019(i)+.25, inland_names(i), 'FontSize',15)
% end
pi9(1).MarkerSize = 12;
pi9(1).Marker = 'o';
pi9(1).MarkerFaceColor = 'b';
errorbar(inland_age,inland_sub2019,inland_sub2019uncert,'LineStyle','none','Color','k')
xlabel('Drainage Age', 'FontSize', 15);
ylabel('Seasonal Subsidence (cm)')
xlim([0, 2500])
ylim([0, 12])
legend(["14C DTLB Average", "Study Area-Wide Average", "Data", "Linear Fit", "95% conf.bounds"],'Location','southwest')
text(10,10.5,['y = ',num2str(fitlininland2019.Coefficients.Estimate(2)),'x +',num2str(fitlininland2019.Coefficients.Estimate(1))],'FontSize',15)
text(10,9.8,['$R^{2} =$',num2str(round(fitlininland2019.Rsquared.Ordinary,2))],'Interpreter','latex','fontsize',15)
text(10, 11.5, 'A.)','FontSize',15,'FontWeight','bold')
title('')
ax = gca;
ax.FontSize = 15;
set(ax,'Box', 'off')

nexttile
yline(5.7,'-k','14C DTLB Average = 5.7','Color','k','FontSize',15)
hold on
yline(5.3,'-r')
pc9 = plot(fitlincoastal2019);
% for i = 1:5
%     text(coastal_age(i)+10, coastal_sub2019(i)+.25, coastal_names(i), 'FontSize',15)
% end
pc9(1).MarkerSize = 12;
pc9(1).Marker = 'o';
pc9(1).MarkerFaceColor = 'b';
errorbar(coastal_age,coastal_sub2019,coastal_sub2019uncert,'LineStyle','none','Color','k')
xlabel('Drainage Age', 'FontSize', 15);
ylabel('Seasonal Subsidence (cm)')
xlim([0, 2500])
ylim([0, 12])
legend('off')
text(10,10.5,['y = ',num2str(fitlincoastal2019.Coefficients.Estimate(2)),'x +',num2str(fitlincoastal2019.Coefficients.Estimate(1))],'FontSize',15)
text(10,9.8,['$R^{2} =$',num2str(round(fitlincoastal2019.Rsquared.Ordinary,2))],'Interpreter','latex','fontsize',15)
text(10, 11.5, 'B.)','FontSize',15,'FontWeight','bold')
title('')
ax = gca;
ax.FontSize = 15;
set(ax,'Box', 'off')

%2021
nexttile
yline(2.7,'-k','14C DTLB Average = 2.7','Color','k','FontSize',15)
hold on
yline(1.8,'-r','Study Area-Wide Average = 1.8','Color','r','FontSize',15)
pi1 = plot(fitlininland2021);
% for i = 1:3
%     text(inland_age(i)+10, inland_sub2021(i)+.25, inland_names(i), 'FontSize',15)
% end
pi1(1).MarkerSize = 12;
pi1(1).Marker = 'o';
pi1(1).MarkerFaceColor = 'b';
errorbar(inland_age,inland_sub2021,inland_sub2021uncert,'LineStyle','none','Color','k')
xlabel('Drainage Age', 'FontSize', 15);
ylabel('Seasonal Subsidence (cm)')
xlim([0, 2500])
ylim([0, 5])
legend('off')
text(10,3.8,['y = ',num2str(fitlininland2021.Coefficients.Estimate(2)),'x +',num2str(fitlininland2021.Coefficients.Estimate(1))],'FontSize',15)
text(10,3.5,['$R^{2} =$',num2str(round(fitlininland2021.Rsquared.Ordinary,2))],'Interpreter','latex','fontsize',15)
text(10, 4.5, 'C.)','FontSize',15,'FontWeight','bold')
title('')
ax = gca;
ax.FontSize = 15;
set(ax,'Box', 'off')

nexttile
yline(1,'-k','14C DTLB Average = 1','Color','k','FontSize',15)
hold on
yline(1.8,'-r')
pc1 = plot(fitlincoastal2021);
% for i = 1:5
%     text(coastal_age(i)+10, coastal_sub2021(i)+.25, coastal_names(i), 'FontSize',15)
% end
pc1(1).MarkerSize = 12;
pc1(1).Marker = 'o';
pc1(1).MarkerFaceColor = 'b';
errorbar(coastal_age,coastal_sub2021,coastal_sub2021uncert,'LineStyle','none','Color','k')
xlabel('Drainage Age', 'FontSize', 15);
ylabel('Seasonal Subsidence (cm)')
xlim([0, 2500])
ylim([0, 5])
legend('off')
text(10,3.8,['y = ',num2str(fitlincoastal2021.Coefficients.Estimate(2)),'x +',num2str(fitlincoastal2021.Coefficients.Estimate(1))],'FontSize',15)
text(10,3.5,['$R^{2} =$',num2str(round(fitlincoastal2021.Rsquared.Ordinary,2))],'Interpreter','latex','fontsize',15)
text(10, 4.5, 'D.)','FontSize',15,'FontWeight','bold')
title('')
ax = gca;
ax.FontSize = 15;
set(ax,'Box', 'off')

