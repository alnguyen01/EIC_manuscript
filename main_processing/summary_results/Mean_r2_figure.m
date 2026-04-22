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
xlabel('EIC Mean $r^{2}$',Interpreter='latex',FontSize=20)
ylabel('ALT Mean $r^{2}$',Interpreter='latex',FontSize=20)
ax = gca;
ax.FontSize =20;
set(ax, 'Box', 'off')

nexttile
linfit = fitlm(ADDT,ALT_r2,'linear');
p = plot(linfit);
for i = 1:length(labels)
    text(ADDT(i)+10,ALT_r2(i)+5,labels(i),'FontSize',20)
end
p(1).MarkerSize = 12;
p(1).Marker = 'o';
p(1).MarkerFaceColor = 'b';
legend('Location',[0.373079243408212 0.722563103926399 0.0973331236577854 0.119331069544059]);
xlim([200,700])
ylim([0,510])
xlabel(['ADDT (',char(176),'C days)'],'FontSize',20)
ylabel('ALT Mean $r^{2}$',Interpreter='latex',FontSize=20)
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
ylabel('EIC Mean $r^{2}$', 'Interpreter', 'latex', 'FontSize', 20);
xlim([200, 700])
ylim([0, 20])
ax = gca;
ax.FontSize = 20;
set(ax, 'Box', 'off')