coherence = h5read('../inputs/geo_ifgramStack.h5','/coherence');
date = h5read('../inputs/geo_ifgramStack.h5','/date');
date_combined = strcat(date(1,:),'_',date(2,:))';
[incmap, ref] = readgeoraster('../../warped_files/NW_drew/20180627.tif');

% totalval = size(coherence,1)*size(coherence,2);
% numbins = sqrt(totalval);
%% remove 2018 poor coh
coherence2018 = coherence(:,:,1:22);
date2018 = date_combined(1:22);

coherence2018(:,:,1:6) =[];
coherence2018(:,:,2) =[];
coherence2018(:,:,3) =[];
coherence2018(:,:,12:14) =[];

date2018(1:6) =[];
date2018(2) =[];
date2018(3) =[];
date2018(12:14) =[];

avecc_2018 = mean(coherence2018,3);
cor_mask_2018 = avecc_2018;

[counts, edges] = histcounts(cor_mask_2018);

otsuthreshold2018 = otsuthresh(counts);
cor_mask_2018(cor_mask_2018<otsuthreshold2018)=nan;
cor_mask_2018(cor_mask_2018>otsuthreshold2018)=1;
mask_2018 = cor_mask_2018;

figure
tiledlayout(2,2)
nexttile
imagesc(avecc_2018')
colorbar
title('Average Coherence of Kept Interferograms for 2018')

nexttile
histogram(avecc_2018)
hold on
xl = xline(otsuthreshold2018,'-',{'Otsu Coherence','threshold =',num2str(otsuthreshold2018)});
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
title('Coherence Distribution')

nexttile
imagesc(mask_2018')
title('Mask')

nexttile
imagesc((avecc_2018.*mask_2018).')
colorbar
title('Masked Coherence')

save("mask_2018.mat","mask_2018",'-mat');
geotiffwrite("otsu_2018.tif",mask_2018',ref);
%% remove 2019 poor coh
coherence2019 = coherence(:,:,23:46);
date2019 = date_combined(23:46);

coherence2019(:,:,1:8) =[];
coherence2019(:,:,2) =[];
coherence2019(:,:,3) =[];
coherence2019(:,:,4) =[];
coherence2019(:,:,5) =[];
coherence2019(:,:,6) =[];
coherence2019(:,:,7) =[];
coherence2019(:,:,8:10) =[];

date2019(1:8) =[];
date2019(2) =[];
date2019(3) =[];
date2019(4) =[];
date2019(5) =[];
date2019(6) =[];
date2019(7) =[];
date2019(8:10) =[];

avecc_2019 = mean(coherence2019,3);
cor_mask_2019 = avecc_2019;

[counts, edges] = histcounts(cor_mask_2019);

otsuthreshold2019 = otsuthresh(counts);
cor_mask_2019(cor_mask_2019<otsuthreshold2019)=nan;
cor_mask_2019(cor_mask_2019>otsuthreshold2019)=1;
mask_2019 = cor_mask_2019;

figure
tiledlayout(2,2)
nexttile
imagesc(avecc_2019')
colorbar
title('Average Coherence of Kept Interferograms for 2019')

nexttile
histogram(avecc_2019)
hold on
xl = xline(otsuthreshold2019,'-',{'Otsu Coherence','threshold =',num2str(otsuthreshold2019)});
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
title('Coherence Distribution')

nexttile
imagesc(mask_2019')
title('Mask')

nexttile
imagesc((avecc_2019.*mask_2019).')
colorbar
title('Masked Coherence')

save("mask_2019.mat","mask_2019",'-mat');
geotiffwrite("otsu_2019.tif",mask_2019',ref);

%% remove 2020 poor coh
coherence2020 = coherence(:,:,47:70);
date2020 = date_combined(47:70);

coherence2020(:,:,1:6) =[];
coherence2020(:,:,2) =[];
coherence2020(:,:,3) =[];
coherence2020(:,:,4) =[];
coherence2020(:,:,5) =[];
coherence2020(:,:,12:14) =[];

date2020(1:6) =[];
date2020(2) =[];
date2020(3) =[];
date2020(4) =[];
date2020(5) =[];
date2020(12:14) =[];

avecc_2020 = mean(coherence2020,3);
cor_mask_2020 = avecc_2020;

[counts, edges] = histcounts(cor_mask_2020);

otsuthreshold2020 = otsuthresh(counts);
cor_mask_2020(cor_mask_2020<otsuthreshold2020)=nan;
cor_mask_2020(cor_mask_2020>otsuthreshold2020)=1;
mask_2020 = cor_mask_2020;


figure
tiledlayout(2,2)
nexttile
imagesc(avecc_2020')
colorbar
title('Average Coherence of Kept Interferograms for 2020')

nexttile
histogram(avecc_2020)
hold on
xl = xline(otsuthreshold2020,'-',{'Otsu Coherence','threshold =',num2str(otsuthreshold2020)});
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
title('Coherence Distribution')

nexttile
imagesc(mask_2020')
title('Mask')

nexttile
imagesc((avecc_2020.*mask_2020).')
colorbar
title('Masked Coherence')

save("mask_2020.mat","mask_2020",'-mat');
geotiffwrite("otsu_2020.tif",mask_2020',ref);

%% remove 2021 poor coh
coherence2021 = coherence(:,:,71:94);
date2021 = date_combined(71:94);

coherence2021(:,:,1:6) =[];
coherence2021(:,:,2) =[];
coherence2021(:,:,3) =[];
coherence2021(:,:,4) =[];
coherence2021(:,:,13:15) =[];

date2021(1:6)=[];
date2021(2)=[];
date2021(3)=[];
date2021(4)=[];
date2021(13:15)=[];

avecc_2021 = mean(coherence2021,3);
cor_mask_2021 = avecc_2021;

[counts, edges] = histcounts(cor_mask_2021);

otsuthreshold2021 = otsuthresh(counts);
cor_mask_2021(cor_mask_2021<otsuthreshold2021)=nan;
cor_mask_2021(cor_mask_2021>otsuthreshold2021)=1;
mask_2021 = cor_mask_2021;

figure
tiledlayout(2,2)
nexttile
imagesc(avecc_2021')
colorbar
title('Average Coherence of Kept Interferograms for 2021')

nexttile
histogram(avecc_2021)
hold on
xl = xline(otsuthreshold2021,'-',{'Otsu Coherence','threshold =',num2str(otsuthreshold2021)});
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
title('Coherence Distribution')

nexttile
imagesc(mask_2021')
title('Mask')

nexttile
imagesc((avecc_2021.*mask_2021).')
colorbar
title('Masked Coherence')

save("mask_2021.mat","mask_2021",'-mat');
geotiffwrite("otsu_2021.tif",mask_2021',ref);

%% remove 2022 poor coh
coherence2022 = coherence(:,:,95:120);
date2022 = date_combined(95:120);

coherence2022(:,:,1:11) = [];
coherence2022(:,:,3) =[];
coherence2022(:,:,12:14) =[];

date2022(1:11) =[];
date2022(3) =[];
date2022(12:14) =[];

avecc_2022 = mean(coherence2022,3);
cor_mask_2022 = avecc_2022;

[counts, edges] = histcounts(cor_mask_2022);

otsuthreshold2022 = otsuthresh(counts);
cor_mask_2022(cor_mask_2022<otsuthreshold2022)=nan;
cor_mask_2022(cor_mask_2022>otsuthreshold2022)=1;
mask_2022 = cor_mask_2022;

figure
tiledlayout(2,2)
nexttile
imagesc(avecc_2022')
colorbar
title('Average Coherence of Kept Interferograms for 2022')

nexttile
histogram(avecc_2022)
hold on
xl = xline(otsuthreshold2022,'-',{'Otsu Coherence','threshold =',num2str(otsuthreshold2022)});
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
title('Coherence Distribution')

nexttile
imagesc(mask_2022')
title('Mask')

nexttile
imagesc((avecc_2022.*mask_2022).')
colorbar
title('Masked Coherence')

save("mask_2022.mat","mask_2022",'-mat');
geotiffwrite("otsu_2022.tif",mask_2022',ref);

%% remove 2023 poor coh
coherence2023 = coherence(:,:,121:146);
date2023 = date_combined(121:146);

coherence2023(:,:,1:4) =[];
coherence2023(:,:,2) =[];
coherence2023(:,:,3) =[];
coherence2023(:,:,4) =[];
coherence2023(:,:,5) =[];
coherence2023(:,:,6) =[];
coherence2023(:,:,7) =[];
coherence2023(:,:,8) =[];
coherence2023(:,:,13:15) =[];

date2023(1:4)=[];
date2023(2)=[];
date2023(3)=[];
date2023(4)=[];
date2023(5)=[];
date2023(6)=[];
date2023(7)=[];
date2023(8)=[];
date2023(13:15)=[];

avecc_2023 = mean(coherence2023,3);
cor_mask_2023 = avecc_2023;

[counts, edges] = histcounts(cor_mask_2023);

otsuthreshold2023 = otsuthresh(counts);
cor_mask_2023(cor_mask_2023<otsuthreshold2023)=nan;
cor_mask_2023(cor_mask_2023>otsuthreshold2023)=1;
mask_2023 = cor_mask_2023;

figure
tiledlayout(2,2)
nexttile
imagesc(avecc_2023')
colorbar
title('Average Coherence of Kept Interferograms for 2023')

nexttile
histogram(avecc_2023)
hold on
xl = xline(otsuthreshold2023,'-',{'Otsu Coherence','threshold =',num2str(otsuthreshold2023)});
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
title('Coherence Distribution')

nexttile
imagesc(mask_2023')
title('Mask')

nexttile
imagesc((avecc_2023.*mask_2023).')
colorbar
title('Masked Coherence')

save("mask_2023.mat","mask_2023",'-mat');
geotiffwrite("otsu_2023.tif",mask_2023',ref);

% %% remove 2024 poor coh
% coherence2024 = coherence(:,:,147:170);
% date2024 = date_combined(147:170);
% 
% coherence2024(:,:,1:8)=[];
% coherence2024(:,:,2)=[];
% coherence2024(:,:,3)=[];
% coherence2024(:,:,4:6)=[];
% coherence2024(:,:,9:11)=[];
% 
% date2024(1:8) =[];
% date2024(2) =[];
% date2024(3) =[];
% date2024(4:6) =[];
% date2024(9:11) =[];
% 
% avecc_2024 = mean(coherence2024,3);
% cor_mask_2024 = avecc_2024;
% 
% [counts, edges] = histcounts(cor_mask_2024);
% 
% otsuthreshold2024 = otsuthresh(counts);
% cor_mask_2024(cor_mask_2024<otsuthreshold2024)=nan;
% cor_mask_2024(cor_mask_2024>otsuthreshold2024)=1;
% mask_2024 = cor_mask_2024;
% 
% figure
% tiledlayout(2,2)
% nexttile
% imagesc(avecc_2024')
% colorbar
% title('Average Coherence of Kept Interferograms for 2024')
% 
% nexttile
% histogram(avecc_2024)
% hold on
% xl = xline(otsuthreshold2024,'-',{'Otsu Coherence','threshold =',num2str(otsuthreshold2024)});
% xl.LabelVerticalAlignment = 'middle';
% xl.LabelHorizontalAlignment = 'center';
% title('Coherence Distribution')
% 
% nexttile
% imagesc(mask_2024')
% title('Mask')
% 
% nexttile
% imagesc((avecc_2024.*mask_2024).')
% colorbar
% title('Masked Coherence')
% 
% save("mask_2024.mat","mask_2024",'-mat');
% geotiffwrite("otsu_2024.tif",mask_2024',ref);
% 
% %% remove 2025 poor coh
% coherence2025 = coherence(:,:,171:191);
% date2025 = date_combined(171:191);
% 
% coherence2025(:,:,1:8) =[];
% coherence2025(:,:,2) =[];
% coherence2025(:,:,3) =[];
% 
% date2025(1:8) =[];
% date2025(2) =[];
% date2025(3) =[];
% 
% avecc_2025 = mean(coherence2025,3);
% cor_mask_2025 = avecc_2025;
% 
% [counts, edges] = histcounts(cor_mask_2025);
% 
% otsuthreshold2025 = otsuthresh(counts);
% cor_mask_2025(cor_mask_2025<otsuthreshold2025)=nan;
% cor_mask_2025(cor_mask_2025>otsuthreshold2025)=1;
% mask_2025 = cor_mask_2025;
% 
% figure
% tiledlayout(2,2)
% nexttile
% imagesc(avecc_2025')
% colorbar
% title('Average Coherence of Kept Interferograms for 2025')
% 
% nexttile
% histogram(avecc_2025)
% hold on
% xl = xline(otsuthreshold2025,'-',{'Otsu Coherence','threshold =',num2str(otsuthreshold2025)});
% xl.LabelVerticalAlignment = 'middle';
% xl.LabelHorizontalAlignment = 'center';
% title('Coherence Distribution')
% 
% nexttile
% imagesc(mask_2025')
% title('Mask')
% 
% nexttile
% imagesc((avecc_2025.*mask_2025).')
% colorbar
% title('Masked Coherence')
% 
% save("mask_2025.mat","mask_2025",'-mat');
% geotiffwrite("otsu_2025.tif",mask_2025',ref);

%% combine back
% coherence_good = cat(3,coherence2018,coherence2019,coherence2020,coherence2021,coherence2022,coherence2023,coherence2024,coherence2025);
% dates_kept = cat(1,date2018,date2019,date2020,date2021,date2022,date2023,date2024,date2025);

coherence_good = cat(3,coherence2018,coherence2019,coherence2020,coherence2021,coherence2022,coherence2023);
dates_kept = cat(1,date2018,date2019,date2020,date2021,date2022,date2023);

%checking discrepancy
keptlist = readlines('../inputs/keptlist');
diff = setdiff(keptlist,dates_kept);
%% average coh

avecc_good = mean(coherence_good,3);
[counts, edges] = histcounts(avecc_good);
otsuthresholdgood = otsuthresh(counts);

avecc_good(avecc_good<otsuthresholdgood) = nan;

cor_mask_good = avecc_good;
cor_mask_good(avecc_good<otsuthresholdgood)=nan;
cor_mask_good(avecc_good>otsuthresholdgood)=1;
mask_good = cor_mask_good;

squeeze_avecc_good = squeeze(avecc_good);
row_squeezed_avecc = (squeeze_avecc_good(:));
percent95_avecc_good = prctile(row_squeezed_avecc,99); %99

findrefpixels = avecc_good>percent95_avecc_good;
highcohpixels = findrefpixels.*avecc_good;

figure
tiledlayout(2,2)
nexttile

imagesc(avecc_good')
colorbar
% title('Average Coherence of Kept Interferograms for 2018-2025 (n=78)')
title('Average Coherence of Kept Interferograms for 2018-2023 (n=64)')

nexttile
histogram(avecc_good)
title('Coherence Distribution')

nexttile
imagesc(highcohpixels')
title('Pixels above coherence 99th percentile = ',num2str(percent95_avecc_good))
colorbar

nexttile
imagesc(mask_good')
title('Masking out pixels below Otsu coherence threshold')
%% total average coh
avecc = mean(coherence,3);
[counts, edges] = histcounts(avecc);
otsuthreshold = otsuthresh(counts);

avecc(avecc<otsuthreshold) = nan;

cor_mask = avecc;
cor_mask(avecc<otsuthreshold)=nan;
cor_mask(avecc>otsuthreshold)=1;
mask = cor_mask;

squeeze_avecc = squeeze(avecc);
row_squeezed_avecc_total = (squeeze_avecc(:));
percent95_avecc = prctile(row_squeezed_avecc_total,99); %99

findrefpixels_total = avecc>percent95_avecc;
highcohpixels_total = findrefpixels.*avecc;

figure
tiledlayout(2,2)
nexttile
imagesc(avecc')
colorbar
title('Average Coherence of ALL Interferograms for 2018-2025 (n=191)')

nexttile
histogram(avecc)
title('Coherence Distribution')

nexttile
imagesc(highcohpixels_total')
title('Pixels above coherence 99th percentile = ',num2str(percent95_avecc))
colorbar

nexttile
imagesc(mask')
title('Masking out pixels below Otsu coherence threshold')

%combined figure
figure
tiledlayout(2,4)
nexttile
imagesc(avecc_good')
colorbar
title('Average Coherence of Kept Interferograms for 2018-2023 (n=64)')

nexttile
histogram(avecc_good)
title('Coherence Distribution')

nexttile
imagesc(avecc')
colorbar
title('Average Coherence of ALL Interferograms for 2018-2025 (n=191)')

nexttile
histogram(avecc)
title('Coherence Distribution')

nexttile
imagesc(highcohpixels')
title('Pixels above coherence 99th percentile = ',num2str(percent95_avecc_good))
colorbar

nexttile
imagesc(mask_good')
title('Masking out pixels below Otsu coherence threshold')

nexttile
imagesc(highcohpixels_total')
title('Pixels above coherence 99th percentile = ',num2str(percent95_avecc))
colorbar

nexttile
imagesc(mask')
title('Masking out pixels below Otsu coherence threshold')

%% export
[inc_map,ref] = readgeoraster('../../warped_files/NW_drew/warped_S1AA_20180907T165937_20180919T165937_VVP012_INT40_G_weF_A162_inc_map.tif');
highcohpixels(highcohpixels==0)=-9999;
highcohpixels_total(highcohpixels_total==0)=-9999;

geotiffwrite('highcoh2018_2023_kept.tif',highcohpixels',ref);
geotiffwrite('highcoh2018_2025_total.tif',highcohpixels_total',ref);