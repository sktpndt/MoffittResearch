%% Import Data (TI Counts + LRF)
ti_all= readtable("/Users/saketpandit/Documents/Moffitt/Project/Data/TumorAndImmuneCounts_Anonymized.xlsx");
ti_outcome = readtable("/Users/saketpandit/Documents/Moffitt/Project/Data/TumorAndImmuneCounts_LungOutcomeSubset_Anonymized.xlsx");

ti_all = join(ti_outcome(1:end-9,:), ti_all);

% Extracting lung cancer (including all LC due to small N)
cond = ti_all.DiseaseSite == "Lung";
% & ti_all.HistologicalSubtypes ~= "LUSC"... 
% & ti_all.HistologicalSubtypes ~= "";

% subset of interest
ti_sub = ti_all(cond,:);

%% initCond
% % Setting initCond
% AT = ti_sub.TotalAnti_TumorCase3;
% T = ti_sub.Total_TumorCells;
% LocalFailure = ti_sub.LocalFailure;
% % initCond = [AT/max(T)*20, T/max(T)*20];
% 
% % Scaling values
% initCond = [AT/max(T)*20, T/max(T)*20];

%% Data Visualization: RSI
% Boxplot
subplot(1,2,1)
boxplot(ti_sub.RSI)
title RSI Spread

hold on

% Histogram 
subplot(1,2,2)
histogram(ti_sub.RSI)
title RSI Distribution
