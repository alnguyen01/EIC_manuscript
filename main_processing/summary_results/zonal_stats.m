zonal_E = readtable("zonal_stats.xls","Sheet","Subsidence");
zonal_EIC = readtable("zonal_stats.xls","Sheet","EIC");
zonal_DEM = readtable("zonal_stats.xls","Sheet","DEM");
zonal_slope = readtable("zonal_stats.xls","Sheet","Slope");
%%
E_data = [zonal_E.MIN(2) zonal_E.PCT25(2) zonal_E.MEDIAN(2) zonal_E.PCT75(2) zonal_E.MAX(2);
    zonal_E.MIN(1) zonal_E.PCT25(1) zonal_E.MEDIAN(1) zonal_E.PCT75(1) zonal_E.MAX(1);
    zonal_E.MIN(3) zonal_E.PCT25(3) zonal_E.MEDIAN(3) zonal_E.PCT75(3) zonal_E.MAX(3)]; %low to high (carex aquatillis -> wet sedge -> tussock tundra

EIC_data = [zonal_EIC.MIN(2) zonal_EIC.PCT25(2) zonal_EIC.MEDIAN(2) zonal_EIC.PCT75(2) zonal_EIC.MAX(2);
    zonal_EIC.MIN(1) zonal_EIC.PCT25(1) zonal_EIC.MEDIAN(1) zonal_EIC.PCT75(1) zonal_EIC.MAX(1);
    zonal_EIC.MIN(3) zonal_EIC.PCT25(3) zonal_EIC.MEDIAN(3) zonal_EIC.PCT75(3) zonal_EIC.MAX(3)]; %low to high (carex aquatillis -> wet sedge -> tussock tundra

DEM_data = [zonal_DEM.MIN(3) zonal_DEM.PCT25(3) zonal_DEM.MEDIAN(3) zonal_DEM.PCT75(3) zonal_DEM.MAX(3);
    zonal_DEM.MIN(2) zonal_DEM.PCT25(2) zonal_DEM.MEDIAN(2) zonal_DEM.PCT75(2) zonal_DEM.MAX(2);
    zonal_DEM.MIN(6) zonal_DEM.PCT25(6) zonal_DEM.MEDIAN(6) zonal_DEM.PCT75(6) zonal_DEM.MAX(6)]; %low to high (carex aquatillis -> wet sedge -> tussock tundra

slope_data = [zonal_slope.MIN(2) zonal_slope.PCT25(2) zonal_slope.MEDIAN(2) zonal_slope.PCT75(2) zonal_slope.MAX(2);
    zonal_slope.MIN(1) zonal_slope.PCT25(1) zonal_slope.MEDIAN(1) zonal_slope.PCT75(1) zonal_slope.MAX(1);
    zonal_slope.MIN(4) zonal_slope.PCT25(4) zonal_slope.MEDIAN(4) zonal_slope.PCT75(4) zonal_slope.MAX(4)]; %low to high (carex aquatillis -> wet sedge -> tussock tundra

E_data = E_data(:, [1 2 2 3 4 4 5]);
EIC_data = EIC_data(:, [1 2 2 3 4 4 5]);
DEM_data = DEM_data(:, [1 2 2 3 4 4 5]);
slope_data = slope_data(:, [1 2 2 3 4 4 5]);

labels = {'Carex Aquatillis', 'Wet Sedge', 'Tussock Tundra'}';

%%

boxfig = figure;
tiledlayout(2,2,"TileSpacing","tight")
nexttile
boxplot(E_data.', 'Whisker', inf);
xticks(1:numel(labels))
xticklabels(labels)
ylabel('Seasonal Subsidence (cm)')
ax = gca;
ax.FontSize = 17.5;
ylim([0, 25])
text(0.6, 23, 'A.)','FontSize',17.5,'FontWeight','bold')

nexttile
boxplot(EIC_data.', 'Whisker', inf);
xticks(1:numel(labels))
xticklabels(labels)
ylabel('%EIC')
ax = gca;
ax.FontSize = 17.5;
ylim([0, 130])
text(0.6, 120, 'B.)','FontSize',17.5,'FontWeight','bold')

nexttile
boxplot(DEM_data.', 'Whisker', inf);
xticks(1:numel(labels))
xticklabels(labels)
ylabel('Elevation (m)')
ax = gca;
ax.FontSize = 17.5;
ylim([0, 15])
text(0.6, 14, 'C.)','FontSize',17.5,'FontWeight','bold')

nexttile
boxplot(slope_data.', 'Whisker', inf);
xticks(1:numel(labels))
xticklabels(labels)
ylabel(['Slope (',char(176),')'])
ax = gca;
ax.FontSize = 17.5;
text(0.6, 6.5, 'D.)','FontSize',17.5,'FontWeight','bold')

%%

E_data = [zonal_E.PCT25(2) zonal_E.MEDIAN(2) zonal_E.PCT75(2);
     zonal_E.PCT25(1) zonal_E.MEDIAN(1) zonal_E.PCT75(1);
     zonal_E.PCT25(3) zonal_E.MEDIAN(3) zonal_E.PCT75(3)]; %low to high (carex aquatillis -> wet sedge -> tussock tundra

EIC_data = [ zonal_EIC.PCT25(2) zonal_EIC.MEDIAN(2) zonal_EIC.PCT75(2) ;
    zonal_EIC.PCT25(1) zonal_EIC.MEDIAN(1) zonal_EIC.PCT75(1) ;
     zonal_EIC.PCT25(3) zonal_EIC.MEDIAN(3) zonal_EIC.PCT75(3) ]; %low to high (carex aquatillis -> wet sedge -> tussock tundra

DEM_data = [ zonal_DEM.PCT25(3) zonal_DEM.MEDIAN(3) zonal_DEM.PCT75(3) ;
     zonal_DEM.PCT25(2) zonal_DEM.MEDIAN(2) zonal_DEM.PCT75(2) ;
     zonal_DEM.PCT25(6) zonal_DEM.MEDIAN(6) zonal_DEM.PCT75(6) ]; %low to high (carex aquatillis -> wet sedge -> tussock tundra

slope_data = [ zonal_slope.PCT25(2) zonal_slope.MEDIAN(2) zonal_slope.PCT75(2) ;
     zonal_slope.PCT25(1) zonal_slope.MEDIAN(1) zonal_slope.PCT75(1) ;
     zonal_slope.PCT25(4) zonal_slope.MEDIAN(4) zonal_slope.PCT75(4) ]; %low to high (carex aquatillis -> wet sedge -> tussock tundra

E_data = E_data(:, [1 1 2 3 3 ]);
EIC_data = EIC_data(:, [1 1 2 3 3] );
DEM_data = DEM_data(:, [1 1 2 3 3 ]);
slope_data = slope_data(:, [1 1 2 3 3 ]);

boxfig2 = figure;
tiledlayout(2,2,'TileSpacing','compact')
nexttile
boxplot(E_data.', 'Whisker', 0);
xticks(1:numel(labels))
xticklabels(labels)
ylabel('Seasonal Subsidence (cm)')
ax = gca;
ax.FontSize = 17.5;
ax.LineWidth = 1.5;
text(0.6, 7, 'A.)','FontSize',17.5,'FontWeight','bold')

bx1 = findobj(gca,'Tag','boxplot');
set(bx1.Children,'LineWidth',1.5)

nexttile
boxplot(EIC_data.', 'Whisker', 0);
xticks(1:numel(labels))
xticklabels(labels)
ylabel('%EIC')
ax = gca;
ax.FontSize = 17.5;
ax.LineWidth = 1.5;
text(0.6, 20, 'B.)','FontSize',17.5,'FontWeight','bold')

bx1 = findobj(gca,'Tag','boxplot');
set(bx1.Children,'LineWidth',1.5)

nexttile
boxplot(DEM_data.', 'Whisker', 0);
xticks(1:numel(labels))
xticklabels(labels)
ylabel('Elevation (m)')
ax = gca;
ax.FontSize = 17.5;
ax.LineWidth = 1.5;
text(0.6, 4.7, 'C.)','FontSize',17.5,'FontWeight','bold')

bx1 = findobj(gca,'Tag','boxplot');
set(bx1.Children,'LineWidth',1.5)

nexttile
boxplot(slope_data.', 'Whisker', 0);
xticks(1:numel(labels))
xticklabels(labels)
ylabel(['Slope (',char(176),')'])
ax = gca;
ax.FontSize = 17.5;
ax.LineWidth = 1.5;
text(0.6, 0.45, 'D.)','FontSize',17.5,'FontWeight','bold')

bx1 = findobj(gca,'Tag','boxplot');
set(bx1.Children,'LineWidth',1.5)
