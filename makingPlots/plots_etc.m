%% load data
% Get data (already generated from LymphSF_PtSpec.m)
load("/Users/saketpandit/Documents/Moffitt/Project/MATLAB/Scripts/Lymphocytes/pt_spec_sf.mat")
%% Fixing SFends
for m = 1:size(SFends, 1)
    SFends{m, 4} = zeros(1)
    if max(SFends{m, 3}) == 0 % never resolves
%     res = max(SFends{m, 3} - min(SFends{m, 3}));
%     if res > 0 % indicating that it does resolve
        SFends{m, 4} = 1; 
    end
end
SFends


%% Plot initCond + postRT
figure(5); clf
% Together
subplot(1,1,1)
% Plot separatrix
% plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1), LineWidth=4, ...
%     Color=[0.8 0.34 0.34])
% hold on

% Plot before and after
x_I0 = I0_orig(:, 1)/max(x)*(Npoints-1);
y_I0 = I0_orig(:, 2)/max(y)*(Npoints-1);
for m = 1:size(I0_orig, 1)
    if postRT(m, 2) > 1e-10
        plot([postRT(m, 1), x_I0(m)], [postRT(m, 2), y_I0(m)],  "k:", LineWidth=0.1);
        hold on
    end
end
pre = scatter(x_I0, y_I0, 200, RSI, 'filled');
post = scatter(postRT(:, 1), postRT(:, 2), 300, RSI, 'filled');
% Set plot window
% colorbar
set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
ax = gca;
ax.FontSize = 20;
xlabel("Effector Cell Count", FontSize = 35)
ylabel("Tumor Cell Count", FontSize = 35)
xlim([1e-4, 10^8])
ylim([1e-10, 1e8])
colormap winter
a = colorbar;
h = ylabel(a, "RSI", Rotation=0, FontSize=25, Position = [5, 0.58, 0])
tit = "Modeling Radiation Therapy";
title(tit, FontSize = 40)



% initCond, separate
% subplot(2, 3, 1)
% pre = scatter(x_I0, y_I0, 200, RSI, 'filled');
% % colorbar
% % set(gca, 'xscale', 'log')
% % set(gca, 'yscale', 'log')
% ax = gca;
% ax.FontSize = 15;
% xlabel("E", FontSize = 20)
% ylabel("T", FontSize = 20)
% % xlim([1e-4, 10^8])
% % ylim([1e-18, 1e7])
% tit = "Pre Radiation Therapy";
% title(tit, FontSize = 15)
% 
% 
% 


%% pt spec SF rat

figure(8); clf  
subplot(1,1,1)
% Plot separatrix

% Npoints = 30;
% x = linspace(0,3.5,Npoints);
% y = linspace(0,450,Npoints);

% just playing around
% Npoints = 30;
% x = linspace(0,300,Npoints);
% y = linspace(0,500,Npoints);


% Kuznetzov parameters
sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;   

% Saket's adjustments
alpha = 7.6;
sigma = 1.3;

plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1), LineWidth=7, ...
    Color=[0.4 0.34 0.34])
hold on

post = gscatter(postRT(:, 1), postRT(:, 2), LocalFailure, ['b', 'r'], ".", 50)

% separatrix w/ contour plot
dx = x(2)-x(1);
dy = y(2)-y(1);
[X, Y] = meshgrid(x,y);
G = rhs([],[reshape(X,1,[]); reshape(Y,1,[])]);
U = reshape(G(1,:),Npoints,Npoints);
V = reshape(G(2,:),Npoints,Npoints)*dx/dy;
N = sqrt(U.^2+V.^2);
U = U./N; V = V./N;
[X1, Y1] = meshgrid(0:Npoints-1,0:Npoints-1);
q = streamslice(X1,Y1,U,V); 
% q.Color = [0 0 0]; 
% q.AutoScaleFactor = 0.5;

% Window settings
set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
legend(post, 'Locoregional Control', 'Locoregional Failure', 'fontsize', 30)
ax = gca
ax.FontSize = 20
xlabel("Effector Cell Count", FontSize = 30)
ylabel("Tumor Cell Count", FontSize = 30)
xlim([1e-4, 10^2.1])
ylim([1e-10, 1e3])
tit = "Patient-Specific R Values"
% subtit = "Sensitivity = " + round(17/20, 2) + ", Specificity = " + round(29/30, 2) + ", Accuracy = " + round((29+17)/(29+17+4), 2)
% tit = "SF2_E = SF2_T / " + 0.9;
title(tit, FontSize = 35)
% subtitle(subtit)


%% PT-Spec Accuracy
preds_ptsp = zeros(size(SFends, 1),1);
for m = 1:size(SFends, 1)
    ypath = ypaths{m, 1};
    if ypath(end) > 5
        preds_ptsp(m) = 1;
    end
end
figure(30);clf
ptacc = confusionchart(LocalFailure, preds_ptsp);
ptacc.Title = "LRF Classification Using Patient-Specific R"
ptacc.XLabel = "Predictions"
ptacc.YLabel = "LRF (1 = Failure)"
ax = gca
ax.FontSize = 50


%% PT-Spec AUC 
% Effector Cells
class1 = postRT(LocalFailure==1, 1); class2 = postRT(LocalFailure==0, 1);
resultsv =roc_curve(class1,class2, 0, 0);
% Tumor Cells
class1 = postRT(LocalFailure==0, 2); class2 = postRT(LocalFailure==1, 2);
resultsh =roc_curve(class1,class2, 0, 0);

% Plot
figure(2);clf
plot(1-resultsv.curve(:,2), resultsv.curve(:,1), 'm-', LineWidth=3);
hold on
plot(1-resultsh.curve(:,2), resultsh.curve(:,1), 'c-', LineWidth=3);
plot(linspace(0,1,10), linspace(0,1,10), 'k--')
% Plot settings
axis square
ax = gca;
ax.FontSize = 18;
xlabel('1 - Specificity', FontSize = 20);
ylabel('Sensitivity', FontSize = 20);
legend("Effector Cells (AUC = " + round(resultsv.param.AROC, 3) + ")", ...
    "Tumor Cells (AUC = " + round(resultsh.param.AROC, 3) + ")", ...
    "Unit Line", ...
    Location = 'southeast')
title("ROC for Patient-Specific R")
hold off





%% INITIAL CONDIITONS: T v E
figure(3); clf
subplot(1,2,1)
ice = gscatter(AT, T, LocalFailure, ['b', 'r'], ".", 60)
ax = gca
ax.FontSize = 20
xlabel("E", FontSize = 25)
ylabel("T", FontSize = 25)
legend(['Locoregional Control'; 'Locoregional Failure'], 'FontSize', 25, Location='southwest')
tit = "Initial Conditions"
% title(tit, FontSize= 25)
% sgtitle(tit, 'fontsize', 25)

subplot(1,2,2)
icp =  gscatter(AT, RSI, LocalFailure, ['b', 'r'], ".", 60)
ax = gca
ax.FontSize = 20
xlabel("E", FontSize = 25)
ylabel("RSI", FontSize = 25)
legend(['Locoregional Control'; 'Locoregional Failure'], 'FontSize', 25, Location='southwest')
% tit = "Initial Conditions"
% title(tit, FontSize= 25)

%% INITIAL CONDITIONS: RSI vs E (with divider?)

figure(4); clf
icp =  gscatter(AT, RSI, LocalFailure, ['b', 'r'], ".", 40)
xline(0.2017, Color=[0 0.4 0.5], LineWidth=5, LineStyle='--')
ex = [0.4 0.5];
why = [0.3  0.3];
annotation('arrow',ex,why, 'Color',[0 0.4 0.5], LineWidth=3 );
ex = [0.35 0.25];
why = [0.3  0.3];
annotation('arrow',ex,why, 'Color',[0 0.4 0.5], LineWidth=3 );

yline(0.37, Color=[0 0.6 0.2], LineWidth = 5, LineStyle='--')
ex = [0.75 0.75];
why = [0.6  0.75];
annotation('arrow',ex,why, 'Color',[0 0.6 0.2], LineWidth=3 );
ex = [0.75 0.75];
why = [0.5  0.35];
annotation('arrow',ex,why, 'Color',[0 0.6 0.2], LineWidth=3 );

% Window params
ax = gca
ax.FontSize = 23
xlabel("Effector Cell Count", FontSize = 30)
ylabel("RSI", FontSize = 30)
legend(['Locoregional Control'; 'Locoregional Failure'], 'FontSize', 25, Location='southwest')
tit = "Initial Conditions"
title(tit, FontSize= 35)


%% AUC from RSI v E
% len = 10;
% vlines = linspace(0, max(AT), len);
% aucs_inv = zeros(len, 4);
% for l = 1:size(vlines, 2)
%     lin = vlines(l);
%     preds = zeros(size(AT,1), 1);
%     for m = 1:size(preds)
%         if AT(m) < lin
%             preds(m) = 1;
%         end
%     end
%     resv = roc_curve(preds, LocalFailure, 0, 0);
%     aucs_inv(l, 1) = lin;
%     aucs_inv(l, 2) = 1 - resv.param.Speci;
%     aucs_inv(l, 3) = resv.param.Sensi;
% %     aucs_in{l, 4} = auc;
% end
% aucs_inv
% 
% hlines = linspace(0, max(RSI), len);
% aucs_inh = zeros(len, 4);
% for l = 1:size(hlines, 2)
%     lin = hlines(l);
%     preds = zeros(size(AT,1), 1);
%     for m = 1:size(preds)
%         if RSI(m) < lin
%             preds(m) = 1;
%         end
%     end
%     resh = roc_curve(preds, LocalFailure, 0, 0);
%     aucs_inh(l, 1) = lin;
%     aucs_inh(l, 2) = 1 - resh.param.Speci;
%     aucs_inh(l, 3) = resh.param.Sensi;
% %     aucs_in{l, 4} = auc;
% end
% aucs_inh
% 
% figure(10);clf
% plot(aucs_inv(:, 2), aucs_inv(:, 3), 'm-', LineWidth = 6)
% hold on
% plot(aucs_inh(:, 2), aucs_inh(:, 3), 'c-', LineWidth = 3)
% xlim([0, 1])
% ylim([0, 1])
% axis square
% % plot(x, y, LineWidth = 5);
% text(0.82, 0.15, " AUC" + newline + 0.508, FontSize = 20);
plot(x, x, "k--");

ax = gca;
ax.FontSize = 18;
xlabel('1 - specificity', FontSize = 20);
ylabel('Sensitivity', FontSize = 20);
tit = 'ROC Curves for Initial Conditions';
legend('RSI (AUC = 0.508)', 'Effector Cell Count (AUC = 0.514)', 'Unit Line', Location = 'southeast')
title(tit, FontSize = 25)
% RSI AUC = 0.508, EC AUC = 0.5140


%%

figure(27);clf
plot(aucs_in{:, 2}, aucs_in{:, 3})

%% Plot SF_e vs SF_t
figure(6); clf
SF2_e = SF2_t ./ rstar;
scat = gscatter(SF2_t, SF2_e, LocalFailure, ['b', 'r'], '.', 60)
hold on
for m = 1:size(I0_orig, 1)
    if SFends{m, 4} == 1
        plot(SF2_t(m), SF2_e(m), 'k.', markersize = 65, LineWidth=5);
    end
end
line(linspace(0,1,10), linspace(0,1,10), 'Color', 'black', "LineStyle", ":")
% colormap winter
% colorbar
xlabel("SF_T")
ylabel("SF_E")
axis square
legend(scat, ["Locoregional Control";"Locoregional Failure"], Location='southeast')
ax = gca;
ax.FontSize = 20;
tit = "Survival Fractions"
title(tit, FontSize = 20);

%% Plotting final separations
% Plotting bubbles
