%% Import Data (TI Counts + LRF)
ti_all= readtable("/Users/saketpandit/Documents/Moffitt/Project/Data/TumorAndImmuneCounts_Anonymized.xlsx");
ti_outcome = readtable("/Users/saketpandit/Documents/Moffitt/Project/Data/TumorAndImmuneCounts_LungOutcomeSubset_Anonymized.xlsx");

ti_all = join(ti_outcome(1:end-9,:), ti_all);

% % Extracting Non-small cell lung cancer
% cond = ti_all.DiseaseSite == "Lung" ...
% & ti_all.HistologicalSubtypes ~= "LUSC"...
% & ti_all.HistologicalSubtypes ~= "";

% subset of interest
ti_sub = ti_all(cond,:);

%% initCond
% Setting initCond
AT = ti_sub.TotalAnti_TumorCase3;
T = ti_sub.Total_TumorCells;
LocalFailure = ti_sub.LocalFailure;
SF_2 = ti_sub.SF_2_;
% initCond = [AT/max(T)*20, T/max(T)*20];

% Scaling values
initCond = [AT/max(T)*20, T/max(T)*20];

% Values after 1 round of radiation therapy
postRT = [AT/max(T)*20 .* (1-SF_2), T/max(T)*20 .* (1-SF_2)];

%% Scatter plot

Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);

icTbl = table(initCond(:,1)./max(x)*(Npoints-1), ...
    initCond(:,2)./max(y)*(Npoints-1),...
    LocalFailure);
icTbl.Properties.VariableNames = ["AT", "T", "LocalFailure"];
s1 = scatter(icTbl, "AT", "T", ColorVariable="LocalFailure")
s1.Marker = "o"


hold on 

rtTbl = table(postRT(:,1)./max(x)*(Npoints-1), ...
    postRT(:,2)./max(y)*(Npoints-1),...
    LocalFailure);
rtTbl.Properties.VariableNames = ["AT", "T", "LocalFailure"];
s2 = scatter(rtTbl, "AT", "T", "filled", ColorVariable="LocalFailure")
s2.Marker = "o"
colormap copper
colorbar

for i = 1:size(icTbl, 1)
    plot([icTbl{i,"AT"} rtTbl{i,"AT"}], [icTbl{i,"T"} rtTbl{i, "T"}], ...
        "k-.")
end

sgtitle("Patient-Specific Radiation Response")
