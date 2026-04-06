% creating a combination ALT product between Whitcomb and Yi & Kimball
% issue is that Whitcomb product is higher resolution (~30m resolution)
% but masked out a lot due to open water or no data
% Yi and Kimball have lower resolution

[W, X] = readgeoraster("warped_upscaled_alt_2015.tif","OutputType","double","CoordinateSystemType","geographic");
W(W==-9999) = NaN; %remove NoData values
W(W==1) = NaN; %remove open water pixels


[Y, Z] = readgeoraster("warped_Alaska_active_layer_thickness_1km_2001-2015_ALT.tif","OutputType","double","CoordinateSystemType","geographic");
Y(Y==-999) = NaN; %remove NoData values
Y = Y(:,:,15); %just 2015

naz = size(W,1); nr = size(W,2);

maskW = zeros(naz,nr);
maskY = zeros(naz,nr);

for i = 1:naz
    for j = 1:nr
        if isnan(W(i,j))
            maskW(i,j) = 1; %no data
        end
        if ~isnan(W(i,j))
            maskW(i,j) = 0; % real data
        end
    end
end

for i = 1:naz
    for j = 1:nr
        if isnan(Y(i,j))
            maskY(i,j) = 1; %no data
        end
        if ~isnan(Y(i,j))
            maskY(i,j) = 0; % real data
        end
    end
end

combinedmask = maskW+maskY;

%%

combinedmask(combinedmask==0) = NaN; %data in both
combinedmask(combinedmask==2) = NaN; %no data in both
filling = Y.*combinedmask;

%setup to add filling + whitcomb
fillingzero = filling;
for i = 1:naz
    for j = 1:nr
        if isnan(filling(i,j))
            fillingzero(i,j) = 0; %no data
        end
    end
end

Wzero = W;
for i = 1:naz
    for j = 1:nr
        if isnan(W(i,j))
            Wzero(i,j) = 0; %no data
        end
    end
end

combinedALTproduct = Wzero + fillingzero;
combinedALTproduct(combinedALTproduct==0) = NaN;

%% Combine Yi + Kimball w/ Whitcomb product for 2014; then do std dev between the two combined products for uncertainty

[A, B] = readgeoraster("warped_upscaled_alt_2014.tif","OutputType","double","CoordinateSystemType","geographic");
A(A==-9999) = NaN; %remove NoData values
A(A==1) = NaN; %remove open water pixels


[C, D] = readgeoraster("warped_Alaska_active_layer_thickness_1km_2001-2015_ALT.tif","OutputType","double","CoordinateSystemType","geographic");
C(C==-999) = NaN; %remove NoData values
C = C(:,:,14); %just 2014

maskA = zeros(naz,nr);
maskC = zeros(naz,nr);

for i = 1:naz
    for j = 1:nr
        if isnan(A(i,j))
            maskA(i,j) = 1; %no data
        end
        if ~isnan(A(i,j))
            maskA(i,j) = 0; % real data
        end
    end
end

for i = 1:naz
    for j = 1:nr
        if isnan(C(i,j))
            maskC(i,j) = 1; %no data
        end
        if ~isnan(C(i,j))
            maskC(i,j) = 0; % real data
        end
    end
end

combinedmask2014 = maskA+maskC;

combinedmask2014(combinedmask2014==0) = NaN; %data in both
combinedmask2014(combinedmask2014==2) = NaN; %no data in both
filling2014 = C.*combinedmask2014;

%setup to add filling + whitcomb
fillingzero2014 = filling2014;
for i = 1:naz
    for j = 1:nr
        if isnan(filling2014(i,j))
            fillingzero2014(i,j) = 0; %no data
        end
    end
end

Azero = A;
for i = 1:naz
    for j = 1:nr
        if isnan(A(i,j))
            Azero(i,j) = 0; %no data
        end
    end
end

combinedALTproduct2014 = Azero + fillingzero2014;
combinedALTproduct2014(combinedALTproduct2014==0) = NaN;

%%

combinedALTproduct2014_2015 = cat(3, combinedALTproduct2014, combinedALTproduct);
combinedALTproductUncertainty_v2 = std(combinedALTproduct2014_2015,0,3,'omitnan');

averageALTproduct = mean(combinedALTproduct2014_2015,3,'omitnan').*100;
combinedALTproductUncertainty_v2 = combinedALTproductUncertainty_v2.*100;

geotiffwrite("averageALTproduct.tif",averageALTproduct,X);
geotiffwrite("combinedALTproductUncertainty_v2.tif",combinedALTproductUncertainty_v2,X);

averageALTproduct_plot = averageALTproduct;
averageALTproduct_plot(isnan(averageALTproduct_plot))=-9999;
geotiffwrite("averageALTproduct_plot.tif",averageALTproduct_plot,X);
