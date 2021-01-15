clear;
close all;
clc;

%% Read in data

% Real Data Load
%{
    Each data set column is as follows: 1) Time(sec) 
                                        2) Pressure (psi) 
                                        3) Top of Top (C) 
                                        4) Bottom of Top (C)
                                        5) Top of Bottom (C)
                                        6) Bottom of Bottom (C)
                                        7) Current (A)
                                        8) Opt Switch
                                        9) Sterling Engine Number	
%}
data_8degDiff = load('StirlingEngineDataGroupF4_8degDiff');
%data_8degDiff = load('Engine3_Temp8');
data_9degDiff = load('StirlingEngineDataGroupF4_9degDiff');
%data_9degDiff = load('Engine3_Temp9');
data_10degDiff = load('StirlingEngineDataGroupF4_10degDiff');
%data_10degDiff = load('Engine3_Temp10');
data_8degDiff(:,2) = data_8degDiff(:,2).*4.448.*144.*9.*(1/0.83612736); % Converting psi to Pa and gage presure to absolute pressure
%data_8degDiff(:,2) = data_8degDiff(:,2) - min(data_8degDiff(:,2)); % To get positive pressure value
data_8degDiff(:,3:6) = data_8degDiff(:,3:6) + 273.15; % Converting C to K
data_9degDiff(:,2) = data_9degDiff(:,2).*4.448.*144.*9.*(1/0.83612736); % Converting psi to Pa and gage presure to absolute pressure
%data_9degDiff(:,2) = data_9degDiff(:,2) - min(data_9degDiff(:,2)); % To get positive pressure value
data_9degDiff(:,3:6) = data_9degDiff(:,3:6) + 273.15; % Converting C to K
data_10degDiff(:,2) = data_10degDiff(:,2).*4.448.*144.*9.*(1/0.83612736); % Converting psi to Pa and gage presure to absolute pressure
%data_10degDiff(:,2) = data_10degDiff(:,2) - min(data_10degDiff(:,2)); % To get positive pressure value
data_10degDiff(:,3:6) = data_10degDiff(:,3:6) + 273.15; % Converting C to K
Fs = 1650;

% Solidworks data load
%{
    Each data set column is as follows: 1) Frame
                                        2) Time (sec) 
                                        3) Linear Displacement (mm)
%}
data_8degDiff_SW = xlsread('powerpiston_displacement8','Sheet1');
data_8degDiff_SW(:,3) = (data_8degDiff_SW(:,3)-min(data_8degDiff_SW(:,3)))./1000; % Converting mm to m
data_9degDiff_SW = xlsread('powerpiston_displacement9','Sheet1');
data_9degDiff_SW(:,3) = (data_9degDiff_SW(:,3)-min(data_9degDiff_SW(:,3)))./1000; % Converting mm to m
data_10degDiff_SW = xlsread('powerpiston_displacement10','Sheet1');
data_10degDiff_SW(:,3) = (data_10degDiff_SW(:,3)-min(data_10degDiff_SW(:,3)))./1000; % Converting mm to m
Fs_SW = 165;

%% RPM calculation
% Used to size the vectors
diff8 = 0;
diff9 = 0;
diff10 = 0;
for i = 1:length(data_8degDiff)-1
    if data_8degDiff(i,8) == 0 && data_8degDiff(i+1,8) == 1
        diff8 = diff8 + 1;
    end 
end
for i = 1:length(data_9degDiff)-1
    if data_9degDiff(i,8) == 0 && data_9degDiff(i+1,8) == 1
        diff9 = diff9 + 1;
    end 
end
for i = 1:length(data_10degDiff)-1
    if data_10degDiff(i,8) == 0 && data_10degDiff(i+1,8) == 1
        diff10 = diff10 + 1;
    end 
end
t_for_rotation8 = zeros(diff8,1);
t_for_rotation9 = zeros(diff9,1);
t_for_rotation10 = zeros(diff10,1);
RPMvec8 = zeros(diff8-1,1);
RPMvec9 = zeros(diff9-1,1);
RPMvec10 = zeros(diff10-1,1);
j = 1;
for i = 1:length(data_8degDiff)-1
    if data_8degDiff(i,8) == 0 && data_8degDiff(i+1,8) == 1
        t_for_rotation8(j) = data_8degDiff(i+1,1)/60;
        if j == 1
            SW_matchup_index8 = i;
            RPMvec8(j) = 1/t_for_rotation8(j);
        else
            RPMvec8(j-1) = 1/(t_for_rotation8(j)-t_for_rotation8(j-1));
        end
        j = j + 1;
    end    
end
j = 1;
for i = 1:length(data_9degDiff)-1
    if data_9degDiff(i,8) == 0 && data_9degDiff(i+1,8) == 1
        t_for_rotation9(j) = data_9degDiff(i+1,1)/60;
        if j == 1
            SW_matchup_index9 = i;
            RPMvec9(j) = 1/t_for_rotation9(j);
        else
            RPMvec9(j-1) = 1/(t_for_rotation9(j)-t_for_rotation9(j-1));
        end
        j = j + 1;
    end    
end
j = 1;
for i = 1:length(data_10degDiff)-1
    if data_10degDiff(i,8) == 0 && data_10degDiff(i+1,8) == 1
        t_for_rotation10(j) = data_10degDiff(i+1,1)/60;
        if j == 1
            SW_matchup_index10 = i;
            RPMvec10(j) = 1/t_for_rotation10(j);
        else
            RPMvec10(j-1) = 1/(t_for_rotation10(j)-t_for_rotation10(j-1));
        end
        j = j + 1;
    end    
end

AvgRPM8 = sum(RPMvec8)/length(RPMvec8);
AvgRPM9 = sum(RPMvec9)/length(RPMvec9);
AvgRPM10 = sum(RPMvec10)/length(RPMvec10);

%% Temperature Difference Calculations
% Bottom of top - Top of bottom

tempDiff8Vec = data_8degDiff(:,5)-data_8degDiff(:,4);
tempDiff9Vec = data_9degDiff(:,5)-data_9degDiff(:,4);
tempDiff10Vec = data_10degDiff(:,5)-data_10degDiff(:,4);

avgTempDiff8 = sum(tempDiff8Vec)/length(tempDiff8Vec);
avgTempDiff9 = sum(tempDiff9Vec)/length(tempDiff9Vec);
avgTempDiff10 = sum(tempDiff10Vec)/length(tempDiff10Vec);

%% Volume Calculations

D = (150-6)/1000; % Chaning mm to m. 150 mm
H = 21/1000; % Chaning mm to m. 21 mm
H_disp = 11/1000;
D_disp = 140/1000;
V_dp = pi*((D/2)^2)*H - pi*((D_disp/2)^2)*H_disp;
D_pp = 15.5/1000; % Chaning mm to m. 15.5 mm
data8_minIndex = find(data_8degDiff(:,2)==min(data_8degDiff(:,2)));
data8_maxIndex_SW = find(data_8degDiff_SW(:,3)==max(data_8degDiff_SW(:,3)));
data9_minIndex = find(data_9degDiff(:,2)==min(data_9degDiff(:,2)));
data9_maxIndex_SW = find(data_9degDiff_SW(:,3)==max(data_9degDiff_SW(:,3)));
data10_minIndex = find(data_10degDiff(:,2)==min(data_10degDiff(:,2)));
data10_maxIndex_SW = find(data_10degDiff_SW(:,3)==max(data_10degDiff_SW(:,3)));

% most of these values were found by looking at the graphs not by
% calulations
pressureVec_RealData8 = (data_8degDiff(data8_minIndex:10:data8_minIndex+(Fs*(113/165)),2));
pressureVec_RealData9 = (data_9degDiff(data9_minIndex:10:data9_minIndex+(Fs*(34/55)),2));
pressureVec_RealData10 = (data_10degDiff(data10_minIndex:10:data10_minIndex+(Fs*(16/33)),2));
Volume8 = V_dp + (data_8degDiff_SW(:,3)); % m^3
Volume9 = V_dp + (data_9degDiff_SW(:,3));
Volume10 = V_dp + (data_10degDiff_SW(:,3));
Volume8_plot = (data_8degDiff_SW(data8_maxIndex_SW:data8_maxIndex_SW+113,3));
Volume9_plot = (data_9degDiff_SW(data9_maxIndex_SW:data9_maxIndex_SW+102,3));
Volume10_plot = (data_10degDiff_SW(data10_maxIndex_SW:data10_maxIndex_SW+80,3));


%% Work net out Calculation
% prepping for trapz analysis
trapz_8_minV_index = find(Volume8_plot==min(Volume8_plot));
trapz_9_minV_index = find(Volume9_plot==min(Volume9_plot));
trapz_10_minV_index = find(Volume10_plot==min(Volume10_plot));

trapzTop8 = [Volume8_plot(trapz_8_minV_index:end,1), (pressureVec_RealData8(trapz_8_minV_index:end,1)-min(pressureVec_RealData8(trapz_8_minV_index:end,1)))];
trapzBottom8 = flip([Volume8_plot(1:trapz_8_minV_index,1), (pressureVec_RealData8(1:trapz_8_minV_index,1)-min(pressureVec_RealData8(1:trapz_8_minV_index,1)))]);
trapzBottom9 = [Volume9_plot(trapz_9_minV_index:end,1), pressureVec_RealData9(trapz_9_minV_index:end,1)-min(pressureVec_RealData9(trapz_9_minV_index:end,1))];
trapzTop9 = flip([Volume9_plot(1:trapz_9_minV_index,1), pressureVec_RealData9(1:trapz_9_minV_index,1)-min(pressureVec_RealData9(1:trapz_9_minV_index,1))]);
trapzBottom10 = [Volume10_plot(trapz_10_minV_index:end,1), pressureVec_RealData10(trapz_10_minV_index:end,1)-min(pressureVec_RealData10(trapz_10_minV_index:end,1))];
trapzTop10 = flip([Volume10_plot(1:trapz_10_minV_index,1), pressureVec_RealData10(1:trapz_10_minV_index,1)-min(pressureVec_RealData10(1:trapz_10_minV_index,1))]);

work8 =trapz(trapzTop8(:,1),trapzTop8(:,2)) -  trapz(trapzBottom8(:,1),trapzBottom8(:,2));
work9 = trapz(trapzTop9(:,1),trapzTop9(:,2)) - trapz(trapzBottom9(:,1),trapzBottom9(:,2));
work10 = trapz(trapzTop10(:,1),trapzTop10(:,2)) - trapz(trapzBottom10(:,1),trapzBottom10(:,2));

%% Carnot Efficiency Calculations

avg_hightemp8 = mean(data_8degDiff(:,5));
avg_lowtemp8 = mean(data_8degDiff(:,4));
avg_hightemp9 = mean(data_9degDiff(:,5));
avg_lowtemp9 = mean(data_9degDiff(:,4));
avg_hightemp10 = mean(data_10degDiff(:,5));
avg_lowtemp10 = mean(data_10degDiff(:,4));

Carnot_th_effic_8 = 1 - (avg_lowtemp8/avg_hightemp8);
Carnot_th_effic_9 = 1 - (avg_lowtemp9/avg_hightemp9);
Carnot_th_effic_10 = 1 - (avg_lowtemp10/avg_hightemp10);

%% Q In
j = 1;
for i = data8_minIndex:data8_minIndex+(Fs*(113/165))
    if data_8degDiff(i,7) > 0.15
        I8(j,1) = data_8degDiff(i,7);
        j = j + 1;
    end
end
t8_1 = (1/Fs)*length(I8);
t8_2 = max(data_8degDiff(data8_minIndex:data8_minIndex+(Fs*(113/165)),1)-min(data_8degDiff(data8_minIndex:data8_minIndex+(Fs*(113/165)),1)))*(length(I8)/length(data8_minIndex:data8_minIndex+(Fs*(113/165))));
j = 1;
for i = data9_minIndex:data9_minIndex+(Fs*(34/55))
    if data_9degDiff(i,7) > 0.15
        I9(j,1) = data_9degDiff(i,7);
        j = j + 1;
    end
end
t9_1 = (1/Fs)*length(I9);
t9_2 = max(data_9degDiff(data9_minIndex:data9_minIndex+(Fs*(34/55)),1)-min(data_9degDiff(data9_minIndex:data9_minIndex+(Fs*(34/55)),1)))*(length(I9)/length(data9_minIndex:data9_minIndex+(Fs*(34/55))));
j = 1;
for i = data10_minIndex:data10_minIndex+(Fs*(16/33))
    if data_10degDiff(i,7) > 0.1
        I10(j,1) = data_10degDiff(i,7);
        j = j + 1;
    end
end
t10_2 = max(data_10degDiff(data10_minIndex:data10_minIndex+(Fs*(16/33)),1)-min(data_10degDiff(data10_minIndex:data10_minIndex+(Fs*(16/33)),1)))*(length(I10)/length(data10_minIndex:data10_minIndex+(Fs*(16/33))));
t10_1 = (1/Fs)*length(I10);
V = 5; % 5 V from the power supply

Q_in8 = (V*mean(I8))*data_10degDiff(end,1);
Q_in9 = (V*mean(I9))*data_10degDiff(end,1);
Q_in10 = (V*mean(I10))*data_10degDiff(end,1);

%% Thermal Efficiency

nth8 = work8/Q_in8;
nth9 = work9/Q_in9;
nth10 = work10/Q_in10;

%% Q Out

Q_out8 = (1-nth8)*Q_in8;
Q_out9 = (1-nth9)*Q_in9;
Q_out10 = (1-nth10)*Q_in10;

%% Work In

W_in8 = (V*mean(I8))*t8_2;
W_in9 = (V*mean(I8))*t9_2;
W_in10 = (V*mean(I8))*t10_2;

%% Work Out

W_out8 = work8 + W_in8;
W_out9 = work9 + W_in9;
W_out10 = work10 + W_in10;

%% Plots

% Pressure
figure
plot(data_8degDiff(:,1),data_8degDiff(:,2),...
    data_9degDiff(:,1),data_9degDiff(:,2),...
    data_10degDiff(:,1),data_10degDiff(:,2),...
    'linewidth',2);
xlim([0 2]);
xlabel('Time(s)');
ylabel('Pressure(Pa)');
title('Internal Pressure vs Time');
legend('8^o Temperature Differential','9^o Temperature Differential',...
    '10^o Temperature Differential','Location','northoutside');

% Current
figure
plot(data_8degDiff(:,1),data_8degDiff(:,7),...
    data_9degDiff(:,1),data_9degDiff(:,7),...
    data_10degDiff(:,1),data_10degDiff(:,7),...
    'linewidth',2);
xlim([0 0.5]);
xlabel('Time(s)');
ylabel('Current(A)');
title('Current vs Time');
legend('8^o Temperature Differential','9^o Temperature Differential',...
    '10^o Temperature Differential','Location','northoutside');

% Temperature Difference
figure
plot(data_8degDiff(:,1),(data_8degDiff(:,5)-data_8degDiff(:,4)),...
    data_9degDiff(:,1),(data_9degDiff(:,5)-data_9degDiff(:,4)),...
    data_10degDiff(:,1),(data_10degDiff(:,5)-data_10degDiff(:,4)),...
    'linewidth',2);
xlim([0 2]);
xlabel('Time(s)');
ylabel('Temperautre(^oK)');
title('Temperature Difference vs Time');
legend('8^o Temperature Differential','9^o Temperature Differential',...
    '10^o Temperature Differential','Location','northoutside');

% PV Diagram
figure
plot(Volume8_plot,pressureVec_RealData8,...
    Volume9_plot,pressureVec_RealData9,...
    Volume10_plot,pressureVec_RealData10,...
    'linewidth',2);
xlabel('Volume(m^3)');
ylabel('Gage Pressure(Pa)');
title('Internal Pressure vs Volume');
legend('8^o Temperature Differential','9^o Temperature Differential',...
    '10^o Temperature Differential','Location','northoutside');
xlim([-0.0005 0.012]);
ylim([-400 420]);

% PV Diagram (10)
figure
plot(Volume10_plot,pressureVec_RealData10,...
    'linewidth',2);
xlabel('Volume(m^3)');
ylabel('Gage Pressure(Pa)');
title('Internal Pressure vs Volume (10^o Temp Difference)');
xlim([-0.0005 0.012]);
ylim([-400 420]);
