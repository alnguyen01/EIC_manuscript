function [addt,Ntemp] = addt_curve(year)
% close all
%ADDT Script
%Calculate NADDT

%% Read in Temperature data, Daymet ADDT Info
%In our case, the temperature data is hourly averages, so get daily avg.

file = strcat('../../temperature_data/',string(year),'.csv');

%Will need to change depending on your temperature data format, but general
%idea is grab daily average temperature
T = readtable(file,VariableNamingRule='preserve');      %reading in data
Time = datetime(T{:,1:3},'InputFormat','yyyy,MM,dd');   %datetime vector
Time.Hour = T{:,2};
T.Time = Time;                    
T(:,1:5) = [];                      %adding datetime vector, removing cols 1-5
TT = table2timetable(T);            %converting to timetable
TTdailyMean = retime(TT,'daily','mean'); %calculating daily mean
TempDailyMean = TTdailyMean{:,4};   %only daily mean temps!

num_days = size(TempDailyMean,1);

%% Set up zero matrices and integrate for ADDT

thaw_count = zeros(365,1);          %number of thaw days
freeze_count = zeros(365,1);        %number of freeze days
ADDT = zeros(365,1);
ADDF1 = zeros(365,1);
ADDF2 = zeros(365,1);

thaw_count(1,1) = 0;
freeze_count(1,1) = 1;
ADDT(1,1) = 0;
ADDF1(1,1) = 0 - TempDailyMean(1,1);
ADDF2(1,1) = 0 - TempDailyMean(1,1);

%Go through yearly temperature record and calculate Accumulated degree days
%of thaw (ADDT) and accumulated degree days of freezing (ADDF)
for i = 2:num_days
    if TempDailyMean(i,1)>0
        thaw_count(i,1) = thaw_count(i-1,1)+1;
        ADDT(i,1) = ADDT(i-1,1)+TempDailyMean(i,1);
        ADDF1(i,1) = ADDF1(i-1,1);
        ADDF2(i,1) = ADDF1(i-1,1);
    end
    if TempDailyMean(i,1)<0
        freeze_count(i,1) = freeze_count(i-1,1)+1;
        ADDT(i,1) = ADDT(i-1,1);
        ADDF1(i,1) = ADDF1(i-1,1)+(0-TempDailyMean(i,1));
        ADDF2(i,1) = ADDF2(i-1,1)+(0-TempDailyMean(i,1));
    end
    if TempDailyMean(i,1) == 0
        ADDT(i,1) = ADDT(i-1,1);
        ADDF1(i,1) = ADDF1(i-1,1);
        ADDF2(i,1) = ADDF1(i-1,1);
        freeze_count(i,1) = freeze_count(i-1,1);
        thaw_count(i,1) = thaw_count(i-1,1);
    end
end

ADDF2(ADDF2<0) = 0;
ADDT(ADDT<0) = 0;

%generate temperature/subsidence proxy curves
temp_curve = sqrt(ADDT)-4.*sqrt(ADDF2);
temp_curve(temp_curve<0) = 0;

N_temp_curve = temp_curve./max(temp_curve);       %Temperature curve with thaw and freeze
NADDT = sqrt(ADDT./max(ADDT));                    %temperature curve just with thaw

addt = NADDT;                                     %NADDT
Ntemp = N_temp_curve;

%plot NADDT curve
figure
plot(addt)
title({'NADDT vs Time for' num2str(year)})
xlabel('Date (month)')
ylabel('NADDT')
xlim([122 275])
set(gca,'XTick',[122 153 183 214 245 275])      %~DoY corresponding to 1st of the month
set(gca,'XTickLabel',[5 6 7 8 9 10])

%output NADDT matrix to read into inversion.m
filename = strcat(string(year),'NADDT','.mat');
save(filename,"addt")


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%