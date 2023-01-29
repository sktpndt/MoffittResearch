%% Reading in the data + Extracting global variables
tic
% load("/Users/saketpandit/Documents/Moffitt/Project/MATLAB/Scripts/T_I_Calibration/ti_calibration.mat")
ti_all= readtable("/Users/saketpandit/Documents/Moffitt/Project/Data/TumorAndImmuneCounts_Anonymized.xlsx");
ti_outcome = readtable("/Users/saketpandit/Documents/Moffitt/Project/Data/TumorAndImmuneCounts_LungOutcomeSubset_Anonymized.xlsx");

ti_all = join(ti_outcome(1:end-9,:), ti_all);



% % Extracting Non-small cell lung cancer
cond = ti_all.DiseaseSite == "Lung";
% & ti_all.HistologicalSubtypes ~= "LUSC"...
% & ti_all.HistologicalSubtypes ~= "";

% subset of interest
ti_sub = ti_all(cond,:);

% Locoregional failure ==> global
LocalFailure = ti_sub.LocalFailure;

% RSI --> Global
RSI = ti_sub.RSI;

% Survival fraction ==> not global anymore!
SF2_t = ti_sub.SF_2_;
SF2_t = SF2_t * 1.5; 
% Multiplying SF2_t by 1.5: 
% - brings all the points up
% - allows me to sweep for SFe < SFt as well as SFe > SFt

% Patient-specific tumor v lymphocyte SF2 (West et. al.)
SF2_rat = readtable("/Users/saketpandit/Documents/Moffitt/Project/MATLAB/Scripts/Lymphocytes/Tumor_Lymphocyte_SF2.csv")
rats = SF2_rat{:,1}./SF2_rat{:,2}
rats = sortrows(rats)

%% Setting initial parameters
% Kuznetzov params
Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);
 
% Kuznetzov parameters
sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;   

% Saket's adjustments
% alpha = 6.7;
% sigma = 1.3;

mu = 0.008;

% alpha = 10;
% gamma = 1;
% beta = 0.0021;
% delta = 0.6;
% mu = 0.01;
% sigma = 1.6;
% sigma = 0.6
% rho = c3_rho_star % from ti_calibration, 0.0931

% Calculate separatrix
[curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);

% Solve ODE
rhs = @(t,x)([sigma+rho*x(1,:).*x(2,:)./(eta+x(2,:))-mu*x(1,:).*x(2,:)-delta*x(1,:);...
          alpha*x(2,:).*(1-beta*x(2,:))-gamma*x(1,:).*x(2,:)]);

options = odeset('Refine',100);
solve = @(init)(ode45(rhs,[0 100],init,options));
% Separatrix Plotting is below


% Setting up I0_orig
AT = ti_sub.TotalAnti_TumorCase3;
T = ti_sub.Total_TumorCells;
% AT = AT/max(AT)
% T = T/max(T)

% NEW SCALING
AT = AT./(T+AT);
T = T./(T+AT);
% Assuming tumor volume of 10^7.5 cells
I0_orig = [AT, T]*10^7.5;
I0 = I0_orig;

% dose vector
dose_vec = ti_sub.FxSize; % dose per fraction

% number of fractions
num_fx_vec = ti_sub.No_Fractions; % number of doses

% Time steps
fx_dt = 2; % time steps for regrowth (one time step = 12 hrs)


% CREATE A VECTOR OF PRE-POST VALUES (so I don't have to run it every time)
% SFt / SFe = 0.9
postRT_9 = zeros(size(I0_orig));
xpaths_9 = cell(size(I0_orig, 1), 1);
ypaths_9 = cell(size(I0_orig, 1), 1);

% Solution vectors
solsRT_9 = cell(size(I0_orig));

% SFt / SFe = 2
postRT_2 = zeros(size(I0_orig));
xpaths_2 = cell(size(I0_orig, 1), 1);
ypaths_2 = cell(size(I0_orig, 1), 1);

% Solution vectors
solsRT_2 = cell(size(I0_orig));

% SFt / SFe = 1
postRT_1 = zeros(size(I0_orig));
xpaths_1 = cell(size(I0_orig, 1), 1);
ypaths_1 = cell(size(I0_orig, 1), 1);

% Solution vectors
solsRT_1 = cell(size(I0_orig));


%% Getting final positions + solutions to ODE
% SFt / SFe = 0.9
% tic
"Starting Big Loop, 0.9"
SF2_e = SF2_t / 0.9;
parfor m = 1:size(I0_orig,1)
% parfor m = 1:10
% for m = 38
%     parfor m = 1
        "Point " + m
        start_time = 0;
        dose = dose_vec(m);
        num_fx = num_fx_vec(m);    
        xpost = 0;
        ypost = 0;
        xpath = 0;
        ypath = 0;
        I0 = I0_orig(m, :); % initial conditions
        x_I0 = I0(1)/max(x)*(Npoints-1);
        y_I0 = I0(2)/max(y)*(Npoints-1);
        % how many fractions to run
            for k = 1:num_fx
%                 "Num_fx = " + k
                % recalculate initCond
                initCond = [I0(1)*SF2_e(m) I0(2)*SF2_t(m)];
                sols = solve(initCond);
                y_vec = sols.y(2,:)/max(y)*(Npoints-1);
                x_vec = sols.y(1,:)/max(x)*(Npoints-1);
               if k < num_fx % is this the final fraction?
                   % updating I0 to treat w/ another dose 
                   I0 = [x_vec(fx_dt)*max(x) y_vec(fx_dt)*max(y)]./(Npoints-1);
                   start_time = start_time + fx_dt;
               elseif k == num_fx % what to do on final fraction
                   % run it out all the way
                   x_rt = initCond(1)/max(x)*(Npoints-1);
                   y_rt = initCond(2)/max(y)*(Npoints-1);
                   xpost = x_rt;
                   ypost = y_rt;
                   xpath = x_vec; 
                   ypath = y_vec;
                   start_time = start_time + fx_dt;
               end
            end
            postRT_9(m, :) = [xpost, ypost];
            xpaths_9{m} = xpath;
            ypaths_9{m} = ypath;
end
% toc
%% Getting final positions + solutions to ODE
% SFt / SFe = 2

% tic
"Starting Big Loop, 2"
SF2_e = SF2_t / 2;
parfor m = 1:size(I0_orig,1)
% parfor m = 1:10
% for m = 38
%     parfor m = 1
        "Point " + m
        start_time = 0;
        dose = dose_vec(m);
        num_fx = num_fx_vec(m);    
        xpost = 0;
        ypost = 0;
        xpath = 0;
        ypath = 0;
        I0 = I0_orig(m, :); % initial conditions
        x_I0 = I0(1)/max(x)*(Npoints-1);
        y_I0 = I0(2)/max(y)*(Npoints-1);
        % how many fractions to run
            for k = 1:num_fx
%                 "Num_fx = " + k
                % recalculate initCond
                initCond = [I0(1)*SF2_e(m) I0(2)*SF2_t(m)];
                sols = solve(initCond);
                y_vec = sols.y(2,:)/max(y)*(Npoints-1);
                x_vec = sols.y(1,:)/max(x)*(Npoints-1);
               if k < num_fx % is this the final fraction?
                   % updating I0 to treat w/ another dose 
                   I0 = [x_vec(fx_dt)*max(x) y_vec(fx_dt)*max(y)]./(Npoints-1);
                   start_time = start_time + fx_dt;
               elseif k == num_fx % what to do on final fraction
                   % run it out all the way
                   x_rt = initCond(1)/max(x)*(Npoints-1);
                   y_rt = initCond(2)/max(y)*(Npoints-1);
                   xpost = x_rt;
                   ypost = y_rt;
                   xpath = x_vec; 
                   ypath = y_vec;
                   start_time = start_time + fx_dt;
               end
            end
            postRT_2(m, :) = [xpost, ypost];
            xpaths_2{m} = xpath;
            ypaths_2{m} = ypath;
end
% toc


%%
% %% Getting final positions + solutions to ODE
% tic
% % SF2_t / SF2_e = 0.7
% "Starting Big Loop, 0.7"
% SF2_e = SF2_t / 0.7;
% for i = 1:length(SF2_e);
%     if SF2_e(i) >= 1;
%         SF2_e(i) = 1;
%     end
% end
% parfor m = 1:size(I0_orig,1)
% % parfor m = 1:10
% % for m = 38
% %     parfor m = 1
%         "Point " + m
%         start_time = 0;
%         dose = dose_vec(m);
%         num_fx = num_fx_vec(m);    
%         xpost = 0;
%         ypost = 0;
%         xpath = 0;
%         ypath = 0;
%         I0 = I0_orig(m, :); % initial conditions
%         x_I0 = I0(1)/max(x)*(Npoints-1);
%         y_I0 = I0(2)/max(y)*(Npoints-1);
%         % how many fractions to run
%             for k = 1:num_fx
% %                 "Num_fx = " + k
%                 % recalculate initCond
%                 initCond = [I0(1)*SF2_e(m) I0(2)*SF2_t(m)];
%                 sols = solve(initCond);
%                 y_vec = sols.y(2,:)/max(y)*(Npoints-1);
%                 x_vec = sols.y(1,:)/max(x)*(Npoints-1);
%                if k < num_fx % is this the final fraction?
%                    % updating I0 to treat w/ another dose 
%                    I0 = [x_vec(fx_dt)*max(x) y_vec(fx_dt)*max(y)]./(Npoints-1);
%                    start_time = start_time + fx_dt;
%                elseif k == num_fx % what to do on final fraction
%                    % run it out all the way
%                    x_rt = initCond(1)/max(x)*(Npoints-1);
%                    y_rt = initCond(2)/max(y)*(Npoints-1);
%                    xpost = x_rt;
%                    ypost = y_rt;
%                    xpath = x_vec; 
%                    ypath = y_vec;
%                    start_time = start_time + fx_dt;
%                end
%             end
%             postRT_9(m, :) = [xpost, ypost];
%             xpaths_9{m} = xpath;
%             ypaths_9{m} = ypath;
% end
% toc

%% Getting final positions + solutions to ODE

% SFt / SFe = 1
% tic
"Starting Big Loop, 1"
SF2_e = SF2_t / 1;
parfor m = 1:size(I0_orig,1)
% parfor m = 1:10
% for m = 38
%     parfor m = 1
        "Point " + m
        start_time = 0;
        dose = dose_vec(m);
        num_fx = num_fx_vec(m);    
        xpost = 0;
        ypost = 0;
        xpath = 0;
        ypath = 0;
        I0 = I0_orig(m, :); % initial conditions
        x_I0 = I0(1)/max(x)*(Npoints-1);
        y_I0 = I0(2)/max(y)*(Npoints-1);
        % how many fractions to run
            for k = 1:num_fx
%                 "Num_fx = " + k
                % recalculate initCond
                initCond = [I0(1)*SF2_e(m) I0(2)*SF2_t(m)];
                sols = solve(initCond);
                y_vec = sols.y(2,:)/max(y)*(Npoints-1);
                x_vec = sols.y(1,:)/max(x)*(Npoints-1);
               if k < num_fx % is this the final fraction?
                   % updating I0 to treat w/ another dose 
                   I0 = [x_vec(fx_dt)*max(x) y_vec(fx_dt)*max(y)]./(Npoints-1);
                   start_time = start_time + fx_dt;
               elseif k == num_fx % what to do on final fraction
                   % run it out all the way
                   x_rt = initCond(1)/max(x)*(Npoints-1);
                   y_rt = initCond(2)/max(y)*(Npoints-1);
                   xpost = x_rt;
                   ypost = y_rt;
                   xpath = x_vec; 
                   ypath = y_vec;
                   start_time = start_time + fx_dt;
               end
            end
            postRT_1(m, :) = [xpost, ypost];
            xpaths_1{m} = xpath;
            ypaths_1{m} = ypath;
end

%% Plotting for R = 0.9 and R = 2
figure(1); clf
% (SFe = SFt / 0.9)
subplot(1,1,1)
post_9 = gscatter(postRT_9(:,1), postRT_9(:,2), LocalFailure, ['b', 'r'], ".", 60)
hold on
% Plot separatrix
% plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1), LineWidth=5, ...
%     Color=[0.8 0.34 0.34])
% legend(post_9, 'Locoregional Control', 'Locoregional Failure')
legend('off')
% Set plot window
% colorbar
set(gca, 'xscale', 'log')
% set(gca, 'Xticklabel', [])
set(gca, 'yscale', 'log')
ax = gca
ax.FontSize = 20
xlabel("Effector Cell Count", FontSize = 30)
ylabel("Tumor Cell Count", FontSize = 30)
xlim([1e-4, 10^8])
ylim([1e-20, 1e4])
legend(post_9, 'Locoregional Control', 'Locoregional Failure', Location='southeast', fontsize=25)
tit = "R = " + 0.9;
title(tit, FontSize = 35)


% (SFe = SFt / 2)
figure(2);clf
subplot(1,1,1)
post_2 = gscatter(postRT_2(:,1), postRT_2(:,2), LocalFailure, ['b', 'r'], ".", 60)
hold on
% Plot separatrix
% plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1), LineWidth=5, ...
%     Color=[0.8 0.34 0.34])

% Set plot window
% colorbar
set(gca, 'xscale', 'log')
% set(gca, 'Xticklabel', [])
set(gca, 'yscale', 'log')
ax = gca
ax.FontSize = 20
xlabel("Effector Cell Count", FontSize = 30)
ylabel("Tumor Cell Count", FontSize = 30)
xlim([1e-4, 10^8])
ylim([1e-20, 1e4])
legend(post_1, 'Locoregional Control', 'Locoregional Failure', Location='southeast', fontsize=25)
tit = "R = " + 2;
title(tit, FontSize = 35)
toc

%% Plotting for R = 1
% (SFe = SFt / 1)
figure(3);clf
subplot(1,1,1)

% Scatterplot of points
post_1 = gscatter(postRT_1(:,1), postRT_1(:,2), LocalFailure, ['b', 'r'], ".", 60)
hold on

% % Plot solution path
% for m = 1:size(I0_orig, 1)
%     "Point " + m
%     % Initial points
%     x_I0 = I0_orig(m, 1)/max(x)*(Npoints-1);
%     y_I0 = I0_orig(m, 2)/max(y)*(Npoints-1);
%     % Final points
%     x_rt = postRT_1(m, 1); % scaling by Npoints already done in for loop
%     y_rt = postRT_1(m, 2);
% %     % Solutions for final points
%     x_vec = xpaths_1{m};
%     y_vec = ypaths_1{m};
%     hold on
%     plot(x_vec,y_vec,'k:','linewidth',0.1); % solsRT
% %     plot([x_rt, x_I0], [y_rt, y_I0],  "b", LineStyle=":", LineWidth=0.001) % Lines connecting the two
%     if LocalFailure(m) == 1 % bubbles for non-LRF
% %         plot(x_I0, y_I0, 'bo', markersize = 12)
%         plot(x_rt, y_rt, 'r.', markersize = 15, LineWidth=2)
%     end
% end

% Plot separatrix
plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1), LineWidth=5, ...
    Color=[0.8 0.34 0.34])

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

% Set plot window
% colorbar
set(gca, 'xscale', 'log')
% set(gca, 'Xticklabel', [])
set(gca, 'yscale', 'log')
ax = gca
ax.FontSize = 20
xlabel("Effector Cell Count", FontSize = 30)
ylabel("Tumor Cell Count", FontSize = 30)
xlim([1e-4, 10^8])
ylim([1e-20, 1e4])
legend(post_1, 'Locoregional Control', 'Locoregional Failure', Location='southeast', fontsize=25)
tit = "R = " + 1;
title(tit, FontSize = 35)


%% Calculating Sensitivity / Specificity
% figure(7); clf
% fin = zeros(size(ypaths_9, 1), 1);

% Calculating predictions for SFt / 0.9
len = 10
vlines = linspace(0, max(postRT_9(:,1)), len);
aucs_sf9v = zeros(len, 4);
for l = 1:size(vlines, 2)
    lin = vlines(l);
    preds = zeros(size(postRT_9,1), 1);
    for m = 1:size(preds)
        if postRT_9(m, 1) < lin
            preds(m) = 1;
        end
    end
    resv = roc_curve(preds, LocalFailure, 0, 0);
    aucs_sf9v(l, 1) = lin;
    aucs_sf9v(l, 2) = 1 - resv.param.Speci;
    aucs_sf9v(l, 3) = resv.param.Sensi;
%     aucs_in{l, 4} = auc;
end

hlines = logspace(0, 6, len);
aucs_sf9h = zeros(len, 4);
for l = 1:size(hlines, 2)
    lin = hlines(l);
    preds = zeros(size(postRT_9,1), 1);
    for m = 1:size(preds)
        if postRT_9(m, 2) > lin
            preds(m) = 1;
        end
    end
    resh = roc_curve(preds, LocalFailure, 0);
    aucs_sf9h(l, 1) = lin;
    aucs_sf9h(l, 2) = 1 - resh.param.Speci;
    aucs_sf9h(l, 3) = resh.param.Sensi;
%     aucs_in{l, 4} = auc;
end

aucs_sf9v
aucs_sf9h

% % Plotting ROC Curves 
% [X9, Y9, T, AUC9] = perfcurve(LocalFailure, preds_9, 1)
% sf9 = plot(X9, Y9, 'Color', [0 0.2 0], LineWidth=5)
% text(0.02, 0.155, " AUC" + newline + round(AUC9, 3), 'FontSize', 17)
% hold on
% [X2, Y2, T, AUC2] = perfcurve(LocalFailure, preds_2, 1)
% sf2 = plot(X2, Y2, 'Color', [0 0.7 0.7], LineWidth=5)
% text(0.72, 0.875, " AUC" + newline + round(AUC2, 3), 'FontSize', 17)
% % Unit Line
% plot(X9, X9, "k--")
% legend('SF_T / 0.9', 'SF_T / 2', 'AUC = 0.5', Location='northwest')
% 
% ax = gca
% ax.FontSize = 18
% xlabel('1 - specificity', FontSize = 20)
% ylabel('Sensitivity', FontSize = 20)
% tit = 'ROC Curves for Uniform Multipliers'
% title(tit, FontSize = 25)

%% Finding solutions to initCond

xpaths_I0 = cell(size(I0_orig, 1), 1);
ypaths_I0 = cell(size(I0_orig, 1), 1);
tic
parfor m = 1:size(I0_orig,1)
% parfor m = 1:10
% for m = 38
%     parfor m = 1
        "Point " + m
        start_time = 0;
%         dose = dose_vec(m);
%         num_fx = num_fx_vec(m);    
%         xpost = 0;
%         ypost = 0;
        xpath = 0;
        ypath = 0;
        I0 = I0_orig(m, :); % initial conditions
        x_I0 = I0(1)/max(x)*(Npoints-1);
        y_I0 = I0(2)/max(y)*(Npoints-1);
        % how many fractions to run
        sols = solve(I0);
        y_vec = sols.y(2,:)/max(y)*(Npoints-1);
        x_vec = sols.y(1,:)/max(x)*(Npoints-1);
        xpath = x_vec; 
        ypath = y_vec;
        xpaths_I0{m} = xpath;
        ypaths_I0{m} = ypath;
end
toc
%                 "Num_fx = " + k
                % recalculate initCond
%                 initCond = [I0(1)*SF2_e(m) I0(2)*SF2_t(m)];
                
%                if k < num_fx % is this the final fraction?
%                    % updating I0 to treat w/ another dose 
%                    I0 = [x_vec(fx_dt)*max(x) y_vec(fx_dt)*max(y)]./(Npoints-1);
%                    start_time = start_time + fx_dt;
%                elseif k == num_fx % what to do on final fraction
%                    % run it out all the way
%                    x_rt = initCond(1)/max(x)*(Npoints-1);
%                    y_rt = initCond(2)/max(y)*(Npoints-1);
%                    xpost = x_rt;
%                    ypost = y_rt;
%                    start_time = start_time + fx_dt;
%                end
%             end
%             postRT_9(m, :) = [xpost, ypost];



%% Plotting separatrix w/ initCond
subplot(1,1,1)
for m = 1:size(I0_orig, 1)
    "Point " + m
    % Initial points
    x_I0 = I0_orig(m, 1)/max(x)*(Npoints-1);
    y_I0 = I0_orig(m, 2)/max(y)*(Npoints-1);
    % Final points
%     x_rt = postRT(m, 1); % scaling by Npoints already done in for loop
%     y_rt = postRT(m, 2);
%     % Solutions for final points
    x_vec = xpaths_I0{m};
    y_vec = ypaths_I0{m};
    hold on
    plot(x_vec,y_vec,'k:','linewidth',0.1); % solsRT
end

% Separatrix
plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1), LineWidth=2)

% Plot before and after
% x_I0 = I0_orig(:, 1)/max(x)*(Npoints-1);
% y_I0 = I0_orig(:, 2)/max(y)*(Npoints-1);
% pre = scatter(x_I0, y_I0, 50, RSI, 'filled')
% hold on
% post = scatter(postRT(:, 1), postRT(:, 2), 75, 'filled')
% colormap winter
% colorbar

% Set plot window
% colorbar
set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
ax = gca
ax.FontSize = 20
xlabel("E", FontSize = 25)
ylabel("T", FontSize = 25)
% xlim([1e-4, 10^2.1])
% ylim([1e-10, 1e5])
% tit = "Effects of radiation therapy, rho = " + c3_rho_star
% tit = "Radiation therapy with Patient-Specific SF_T / SF_E Ratios"
% tit = "SF2_E = SF2_T / " + 0.9;
% title(tit, FontSize = 25)

