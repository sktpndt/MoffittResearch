%% AUC plots for initial Conditions
% Effector Cells
class1 = AT(LocalFailure==1); class2 = AT(LocalFailure==0);
resultsv =roc_curve(class1,class2, 0, 0);
% RSI
class1 = RSI(LocalFailure==0); class2 = RSI(LocalFailure==1);
resultsh = roc_curve(class1,class2, 0, 0);

% Plot
figure(1);clf
plot(1-resultsv.curve(:,2), resultsv.curve(:,1), Color=[0 0.4 0.5], LineWidth=3);
hold on
plot(1-resultsh.curve(:,2), resultsh.curve(:,1), Color=[0 0.6 0.2], LineWidth=3);
plot(linspace(0,1,10), linspace(0,1,10), 'k--')
% Plot settings
axis square
ax = gca;
ax.FontSize = 18;
xlabel('1 - Specificity', FontSize = 20);
ylabel('Sensitivity', FontSize = 20);
legend("Effector Cells (AUC = " + round(resultsv.param.AROC, 3) + ")", ...
    "RSI (AUC = " + round(resultsh.param.AROC, 3) + ")", ...
    "Unit Line", ...
    Location = 'southeast')
title("ROC for Initial Conditions")
hold off

%% AUC for uniform multipliers

% R = 0.9 %

% Effector Cells
class1 = postRT_9(LocalFailure==0, 1); class2 = postRT_9(LocalFailure==1, 1);
resultsv =roc_curve(class1,class2, 0, 0);
% Tumor Cells
class1 = postRT_9(LocalFailure==0, 2); class2 = postRT_9(LocalFailure==1, 2);
resultsh =roc_curve(class1,class2, 0, 0);

% Plot
figure(2);clf
subplot(1,2,1)
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
% legend("Effector Cells (AUC = " + round(resultsv.param.AROC, 3) + ")", ...
%     "Tumor Cells (AUC = " + round(resultsh.param.AROC, 3) + ")", ...
%     "Unit Line", ...
%     Location = 'southeast')
title("ROC for R = 0.9")
hold off

% R = 2 %

% Effector Cells
class1 = postRT_2(LocalFailure==0, 1); class2 = postRT_2(LocalFailure==1, 1);
resultsv =roc_curve(class1,class2, 0, 0);
% Tumor Cells
class1 = postRT_2(LocalFailure==0, 2); class2 = postRT_2(LocalFailure==1, 2);
resultsh =roc_curve(class1,class2, 0, 0);

% Plot
% figure(3);clf
subplot(1,2,2)
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
title("ROC for R = 2")
hold off