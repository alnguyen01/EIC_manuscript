ts_thaw_waterMask = struct2array(load("../../mintpy2019_v7_NW_drew/ReSALT/ts_thaw_waterMask.mat"));
ts_thaw_waterMask(ts_thaw_waterMask==0) = NaN;
ts_thaw = struct2array(load("ts_thaw.mat"));
ts_thaw(ts_thaw==0) = NaN;
dates=h5read('../geo/geo_timeseries_ERA5_ramp_demErr.h5','/date');

for i=1:length(dates)
    t(i) = datetime(str2double(dates(i)),'ConvertFrom','yyyymmdd');
    doys(i,1) = day(t(i),'dayofyear');
end

%% plot present in both only
nr = size(ts_thaw,1);
naz = size(ts_thaw,2);

common_pixels = zeros(nr,naz);


for j = 1:nr
   for k = 1:naz
       if ~isnan(ts_thaw_waterMask(j,k,2))
           if ~isnan(ts_thaw(j,k,2))
               common_pixels(j,k) = 1;
           else
               common_pixels(j,k) = NaN;
           end
       else
           common_pixels(j,k) = NaN;
       end
   end
end

figure
imagesc(common_pixels')
colorbar
title('Common Pixels between waterMask and coherence mask')
%%
test_diff = (ts_thaw_waterMask - ts_thaw).*common_pixels;

figure
tiledlayout(2,4)
for i = 1:size(test_diff,3)
    nexttile
    imagesc(test_diff(:,:,i).')
    c= colorbar;
    xlabel(c, '[cm]')
    title(['Difference of waterMask - no waterMask ',string(t(i))])
end
