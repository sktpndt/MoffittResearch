%% Reading in the data + Extracting global variables
load("/Users/saketpandit/Documents/Moffitt/Project/MATLAB/Scripts/T_I_Calibration/ti_calibration.mat")
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

% Survival fraction ==> also global
SF_2 = ti_sub.SF_2_;

RSI = ti_sub.RSI;

%% tumor immune playground

% Kuznetzov params
Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);
 
    
sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;   

% alpha = 10;
% gamma = 1;
% beta = 0.0021;
% delta = 0.6;
% mu = 0.01;
% sigma = 1.6;
% sigma = 0.6
% rho = c3_rho_star % from ti_calibration, 0.0931


% Setting up I0_orig
AT = ti_sub.TotalAnti_TumorCase3;
T = ti_sub.Total_TumorCells;
AT = AT/max(AT)
T = T/max(T)

% More reasonable values in postRT
I0_orig = [AT, T]*10^6;
I0 = I0_orig;

% dose vector
dose_vec = ti_sub.FxSize; % dose per fraction

% number of fractions
num_fx_vec = ti_sub.No_Fractions; % number of doses

% Time steps
fx_dt = 2; % time steps for regrowth (one time step = 12 hrs)


% Percentage resolved vs unresolved
r = 0; nr = 0; % number resolved (r) / not resolved (nr)

"Starting Big Loop"
ronr_a = zeros(size(I0_orig, 1), 4);
bubbleSize = 11;

% CREATE A VECTOR OF PRE-POST VALUES (so I don't have to run it every time)
postRT = zeros(size(I0_orig));

% Solution vectors
solsRT = cell(size(I0_orig));

% sigma = 0.1
% alphas = linspace(1, 7, 10);
% aucs = zeros(size(alphas, 2),1);
% for z = 1:size(alphas, 2)
%     figure(z); clf
%     alpha = alphas(z);
% %     alpha_auc(z, 1) = alpha
%     [curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);
%     rhs = @(t,x)([sigma+rho*x(1,:).*x(2,:)./(eta+x(2,:))-mu*x(1,:).*x(2,:)-delta*x(1,:);...
%               alpha*x(2,:).*(1-beta*x(2,:))-gamma*x(1,:).*x(2,:)]);
%     options = odeset('Refine',100);
%     solve = @(init)(ode45(rhs,[0 100],init,options));
      for m = 1:size(I0_orig,1)
%     for m = 1:5
        "Run " + m
        start_time = 0;
        I0 = I0_orig(m, :); % initial conditions
        ronr_a(m, 1:2) = I0_orig(m, :);
        ronr_a(m, 4) = LocalFailure(m);
        dose = dose_vec(m);
        num_fx = num_fx_vec(m);
        prtk_t = 1-SF_2(m);
        prtk_e = prtk_t * 1.1 % made up value, change this one
%         subplot(1,1,1)
        % plotting the initial conditions
%         x_I0 = I0(1)/max(x)*(Npoints-1);
%         y_I0 = I0(2)/max(y)*(Npoints-1);
%         plot(x_I0, y_I0,'b.','linewidth',1.5,'markersize',10)
%         if LocalFailure(m) == 0
%             plot(x_I0, y_I0, 'bo', markersize = bubbleSize)
%         end
    
        hold on
        % how many fractions to run
       
        for k = 1:num_fx
            % recalculate initCond
            initCond = [I0(1)*prtk_e I0(2)*prtk_t]; 
            sols = solve(initCond);
            y_vec = sols.y(2,:)/max(y)*(Npoints-1);
            x_vec = sols.y(1,:)/max(x)*(Npoints-1);
            solsRT(m, :) = {x_vec, y_vec};
           if k < num_fx % is this the final fraction?
               % updating I0 to treat w/ another dose 
               I0 = [x_vec(fx_dt)*max(x) y_vec(fx_dt)*max(y)]./(Npoints-1);
               start_time = start_time + fx_dt;
%                plot(initCond(1)/max(x)*(Npoints-1),initCond(2)/max(y)*(Npoints-1),'b+','linewidth',0.5,'markersize',10)
               %plot(x_vec(1:fx_dt),y_vec(1:fx_dt),'k','linewidth',1); 
           elseif k == num_fx % what to do on final fraction
               % run it out all the way
%                plot(x_vec,y_vec,'r-','linewidth',0.1);
               x_rt = initCond(1)/max(x)*(Npoints-1);
               y_rt = initCond(2)/max(y)*(Npoints-1);
               postRT(m, :) = [x_rt, y_rt];
%                plot(x_rt, y_rt ,'r.','linewidth',20,'markersize',10)
%                if LocalFailure(m) == 0
%                    plot(x_rt, y_rt, 'ro', markersize = bubbleSize)
%                end
%                plot([x_rt, x_I0], [y_rt, y_I0],  "b-", LineWidth=0.1, LineStyle="--")
               hold on
               start_time = start_time + fx_dt;
               
               % Percent classification
               if y_vec(end) > 10 % y value for resolution
                   ronr_a(m, 3) = 1;
               end
           end
        end
    end
%     [Xa, Ya, T, AUCa] = perfcurve(ronr_a(:, 4), ronr_a(:, 3), 0);
%     aucs(z) = AUCa;
% end
% alpha_auc = [alphas', aucs]
% alpha_auc


% prtcorrect_a = 0
% for m = 1:size(ronr_a, 1)
%     if ronr_a(m, 3) == ronr_a(m, 4)
%         prtcorrect_a = prtcorrect_a + 1;
%     end
% end
% 
% "Percent correct = " + prtcorrect_a/size(ronr_a, 1)

% rt3_rho_perc = r/(r+nr)
% it3_rho_perc = c3_percs_rho(c3_percs_rho(:,1)==c3_rho_star, 2)/...
%     (c3_percs_rho(c3_percs_rho(:,1)==c3_rho_star, 2)+ ...
%     c3_percs_rho(c3_percs_rho(:,1)==c3_rho_star, 3))
% toc

%% Plotting (pre / post + RSI vals)
figure(5); clf

% Plotting
% % vector of initial conditions
% x_I0 = I0_orig(:, 1)./max(x)*(Npoints-1);
% y_I0 = I0_orig(:, 2)./max(y)*(Npoints-1);
% % vector of post RT
% x_rt = postRT(:, 1); % scaling by Npoints already done in previous for loop
% y_rt = postRT(:, 2);
% % solutions
% x_vec = solsRT{:, 1};
% y_vec = solsRT{:, 2};

% subplot(2,3,[1 2 4 5])
% plot(x_I0, y_I0, 'b.','linewidth',1.5,'markersize', 10)
% hold on
% plot(x_rt, y_rt,'r.','linewidth',20, 'MarkerSize', 10) %postRT
% plot(x_vec,y_vec,'r-','linewidth',0.1); % solsRT
% plot([x_rt, x_I0], [y_rt, y_I0],  "b-", LineWidth=0.05, LineStyle="--") % Lines connecting the two
% if LocalFailure(m) == 0 % bubbles for non-LRF
%     plot(x_I0, y_I0, 'bo', markersize = bubbleSize)
%     plot(x_rt, y_rt, 'ro', markersize = bubbleSize)
% end
% set(gca, 'xscale', 'log')
% set(gca, 'yscale', 'log')
% % tit = "Effects of radiation therapy, rho = " + c3_rho_star
% tit = "Effects of radiation therapy";
% title(tit)
% hold off

% c = 10;
% blues = zeros(c, 3);
% b = linspace(0, 1, c);
% for n = 1:c
%     blues(n, 1) = b(n);
% end
% 
% reds = zeros(c, 3);
% re = linspace(0, 1, c);
% for n = 1:c
%     reds(n, 1) = re(n);
% end


for m = 1:size(I0_orig, 1)
    "Point " + m
    % Initial points
    x_I0 = I0_orig(m, 1)/max(x)*(Npoints-1);
    y_I0 = I0_orig(m, 2)/max(y)*(Npoints-1);
    % Final points
    x_rt = postRT(m, 1); % scaling by Npoints already done in for loop
    y_rt = postRT(m, 2);
    % Solutions
    x_vec = solsRT{m, 1};
    y_vec = solsRT{m, 2};
    
    % Big plot
    subplot(2,3,[1 2 4 5])
    plot(x_I0, y_I0, 'b.','linewidth',1.5,'markersize', 10)
    hold on
    plot(x_rt, y_rt,'r.','linewidth',20, 'MarkerSize', 10) %postRT
    plot(x_vec,y_vec,'r-','linewidth',0.1); % solsRT
    plot([x_rt, x_I0], [y_rt, y_I0],  "b-", LineWidth=0.01, LineStyle="--") % Lines connecting the two
    if LocalFailure(m) == 0 % bubbles for non-LRF
        plot(x_I0, y_I0, 'bo', markersize = 12)
        plot(x_rt, y_rt, 'ro', markersize = 12)
    end
    set(gca, 'xscale', 'log')
    set(gca, 'yscale', 'log')
    xlabel("AT", FontSize = 15)
    ylabel("T", FontSize = 15)
    % tit = "Effects of radiation therapy, rho = " + c3_rho_star
    tit = "Effects of radiation therapy";
    title(tit, FontSize = 20)
%     hold off

    % Individual Plots
%     mSize = RSI(m)*50;
%     bubbleSize = mSize + 1;
    % initial condiitons
    subplot(2,3,3) 
    s1 = scatter(x_I0, y_I0, 55, 'blue', 'filled')
    s1.AlphaData = RSI(m)
    s1.MarkerFaceAlpha = 'flat'
%     plot(x_I0, y_I0, 'b.', 'linewidth', 1.5, 'Color'=)
    hold on
    if LocalFailure(m) == 0 % bubbles for non-LRF
        plot(x_I0, y_I0, 'bo', markersize = 12)
    end
    set(gca, 'xscale', 'log')
    set(gca, 'yscale', 'log')
    xlabel("AT", FontSize = 15)
    ylabel("T", FontSize = 15)
    tit = "Pre-radiation therapy";
    title(tit, FontSize = 20)
%     hold off
    % post RT
    subplot(2,3,6)
    s2 = scatter(x_rt, y_rt, 55, 'red', 'filled')
    s2.AlphaData = RSI(m)
    s2.MarkerFaceAlpha = 'flat'
%     plot(x_rt, y_rt, 'r.', 'linewidth', 1.5, 'markersize', mSize)
    hold on
    if LocalFailure(m) == 0 % bubbles for non-LRF
        plot(x_rt, y_rt, 'ro', markersize = 12)
    end
    set(gca, 'xscale', 'log')
    set(gca, 'yscale', 'log')
    xlabel("AT", FontSize = 15)
    ylabel("T", FontSize = 15)
    tit = "Post-radiation therapy";
    title(tit, FontSize = 20)
%     hold off
end
hold off


%toc

%% Regulatory Cells - Data Viz

% Data Visualization
figure(2);clf
RC = ti_sub.TotalPro_TumorCase3;
RC = RC/max(RC)
subplot(1,2,1)
histogram(RC)
subplot(1,2,2)
boxplot(RC)

%% 3D Plot
% 3D plot: Tumor, AntiTumor, Regulatory Cells
figure(3);clf
for m = 1:size(T)
    plot3(T(m), AT(m), RC(m), 'b.', 'markersize', 15)
    hold on
    if LocalFailure(m) == 0 % bubbles for non-LRF
        plot3(T(m), AT(m), RC(m), 'bo', markersize = 18)
    end
end
clear xlabel ylabel zlabel
xlabel('T', FontSize=20)
ylabel('AT', FontSize=20)
zlabel('RC', FontSize=20)
title('3D Data Visualization')
grid on
hold off


figure(4);clf
% T v AT
subplot(2,2,1)
for m = 1:size(T)
    plot(AT(m), T(m), 'b.', 'markersize', 15)
    hold on
    if LocalFailure(m) == 0 % bubbles for non-LRF
        plot(AT(m), T(m), 'bo', markersize = 18)
    end
end
xlabel("AT", FontSize = 15)
ylabel("T", FontSize = 15)

% T v RC
subplot(2,2,2)
for m = 1:size(T)
    plot(RC(m), T(m), 'b.', 'markersize', 15)
    hold on
    if LocalFailure(m) == 0 % bubbles for non-LRF
        plot(RC(m), T(m), 'bo', markersize = 18)
    end
end
xlabel("RC", FontSize = 15)
ylabel("T", FontSize = 15)


% RC v AT
subplot(2,2,3)
for m = 1:size(T)
    plot(AT(m), RC(m), 'b.', 'markersize', 15)
    hold on
    if LocalFailure(m) == 0 % bubbles for non-LRF
        plot(AT(m), RC(m), 'bo', markersize = 18)
    end
end
xlabel("AT", FontSize = 15)
ylabel("RC", FontSize = 15)
hold off

% AT/RC vs T/AT
subplot(2,2,4)
for m = 1:size(T)
    loglog(AT(m)/RC(m), T(m)/AT(m), ...
        'b.', markersize=15)
    hold on
    if LocalFailure(m) == 0 % bubbles for non-LRF
        loglog(AT(m)/RC(m), T(m)/AT(m), ...
            'bo', markersize = 18)
    end
end
xlabel("AT/RC", FontSize = 15)
ylabel("T/AT", FontSize = 15)
hold off




