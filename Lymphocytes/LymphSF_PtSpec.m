%% 
% WHOLE SCRIPT TAKES ~ 1 HR TO RUN (len = 20)

%% Reading in the data + Extracting global variables
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

% Survival fraction ==> not global anymore!
SF2_t = ti_sub.SF_2_;
SF2_t = SF2_t * 1.5;
% Patient-specific tumor v lymphocyte SF2 (West et. al. 1998)
SF2_rat = readtable("/Users/saketpandit/Documents/Moffitt/Project/MATLAB/Scripts/Lymphocytes/Tumor_Lymphocyte_SF2.csv");
rats = SF2_rat{:,1}./SF2_rat{:,2};
rats = sortrows(rats);
%[min(SF2_t) max(SF2_t)] % / 0.89 is the lowest you can go without SF > 1

% Selecting compatible survival fractions
lowR = 0.8; % 0.6 leaves in 23 rows, 0.4 leaves in 3 rows, 0.5 leaves in 13 rows, 0.8 leaves 43 rows
cond = SF2_t / lowR  <= 1; 
ti_sub = ti_sub(cond,:);


% Locoregional failure ==> global
LocalFailure = ti_sub.LocalFailure;

% RSI --> Global
RSI = ti_sub.RSI;

% Setting up I0_orig
AT = ti_sub.TotalAnti_TumorCase3;
T = ti_sub.Total_TumorCells;

%% Setting initial parameters
% dose vector
dose_vec = ti_sub.FxSize; % dose per fraction

% number of fractions
num_fx_vec = ti_sub.No_Fractions; % number of doses

% Time steps
fx_dt = 2; % time steps for regrowth (one time step = 12 hrs)

% Kuznetzov params
Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);

% dx = x(2)-x(1);
%  
% dy = y(2)-y(1);
% [X, Y] = meshgrid(x,y);
% G = rhs([],[reshape(X,1,[]); reshape(Y,1,[])]);
% U = reshape(G(1,:),Npoints,Npoints);
% V = reshape(G(2,:),Npoints,Npoints)*dx/dy;
% N = sqrt(U.^2+V.^2);
% U = U./N; V = V./N;
% [X1, Y1] = meshgrid(0:Npoints-1,0:Npoints-1);
% q = quiver(X1,Y1,U,V); 
% q.Color = [0 0 0]; 
% q.AutoScaleFactor = 0.5;

 
% Kuznetzov Params
sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;   

% Saket's adjustments
alpha = 7.6;
sigma = 1.3;

% gamma = 1;
% beta = 0.0021;
% delta = 0.6;
% mu = 0.008;
% sigma = 0.3;
% rho = c3_rho_star % from ti_calibration, 0.0931


% Calculate separatrix
[curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);

% Solve ODE
rhs = @(t,x)([sigma+rho*x(1,:).*x(2,:)./(eta+x(2,:))-mu*x(1,:).*x(2,:)-delta*x(1,:);...
          alpha*x(2,:).*(1-beta*x(2,:))-gamma*x(1,:).*x(2,:)]);

options = odeset('Refine',100);
solve = @(init)(ode45(rhs,[0 100],init,options));
    
% NEW SCALING
AT = AT./(T+AT);
T = T./(T+AT);
% Assuming tumor volume of 10^7.5 cells
I0_orig = [AT, T]*10^7.5;
I0 = I0_orig;


% CREATE A VECTOR OF PRE-POST VALUES (so I don't have to run it every time)
postRT = zeros(size(I0_orig));

% Solution vectors
solsRT = cell(size(I0_orig));
xpaths = cell(size(I0_orig, 1), 1);
ypaths = cell(size(I0_orig, 1), 1);

% Ratios
len = 20;
% SFrats = linspace(0.9, 2, len);
SFrats = linspace(lowR, 2, len); % corresponding with low end in previous block
SFends = cell(size(I0_orig, 1), 4);
for i = 1:size(SFends, 1)
    SFends{i, 1} = I0_orig(i, :);
    SFends{i, 2} = SFrats;
    SFends{i, 3} = zeros(len, 1);
    SFends{i, 4} = zeros(1);
end

% Cutoff for escape / resolve
threshold = 5;

%% Sweeping across SF ratios
tic
"Starting Big Loop"
parfor m = 1:size(I0_orig,1)
% for m = 1:size(I0_orig, 1) % for debugging
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
        rnr = 0;
        for rt = 1:size(SFrats, 2)
            q = SFrats(rt);
            SF2_e = SF2_t(m)/q;
            if SF2_e < 1
                I0 = I0_orig(m, :); % initial conditions
                x_I0 = I0(1)/max(x)*(Npoints-1);
                y_I0 = I0(2)/max(y)*(Npoints-1);
            % how many fractions to run
                for k = 1:num_fx
    %                 "Num_fx = " + k
                    % recalculate initCond
                    initCond = [I0(1)*SF2_e I0(2)*SF2_t(m)];
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
                postRT(m, :) = [xpost, ypost];
                xpaths{m} = xpath;
                ypaths{m} = ypath;
                if ypath
                    if ypath(end) > threshold
                         rnr = 1; % represents LRF
                    end
                end
                SFends{m, 3}(rt) = rnr;
            end
        end
    end
toc % ~30min for 13 points

%% Visualizing SFends
figure(3);clf
for m = 1:size(SFends, 1)
    subplot(10, 5, m)
    plot(SFends{m, 2}, SFends{m, 3})
end

%% Does the point make a switch?
for m = 1:size(SFends, 1)
    res = max(SFends{m, 3} - min(SFends{m, 3}));
    if res > 0
        SFends{m, 4} = 1;
    end
end
%% Finding R*
rstar = zeros(size(SFends, 1), 1);
for m = 1:size(SFends, 1)
% for m = 1:50
% for m = 4:6
    m
    ratios = SFends{m, 2};
    lrf = SFends{m, 3};
    preds = [ratios.' lrf];
    sortrows(preds, 1);
    if LocalFailure(m) == 1
        for k = 1:size(preds, 1)
            if preds(k, 2) == 1
                rstar(m) = preds(k, 1); % first value that shifts point to escape
                break
            else
                rstar(m) = 1;
            end
        end
    else
        for k = 1:size(preds, 1)
            if preds(k, 2) == 1
                if k == 1
                    rstar(m) = preds(1,1); % just using the first ratio value (0.9)
                    break
                else
                    rstar(m) = preds(k-1, 1); % largest ratio that tumor still resolves
                    break
                end
            else
                rstar(m) = 1; % stubborn point, never resolves
            end
        end
    end
end
rstar

%% Distribution of rstar values
figure(1); clf
histogram(rstar)
ax = gca
ax.FontSize = 18
title("Patient-Specific Effector Cell Survival Fractions", FontSize = 20)

%% Solving ODE with patient specific SF ratios (rstar)
%tic
"Starting Big Loop"
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
        rnr = 0;
        SF2_e = SF2_t/rstar(m)
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
            postRT(m, :) = [xpost, ypost];
            xpaths{m} = xpath;
            ypaths{m} = ypath;
    end
%toc


%% Plotting with pt-spec SF ratios %% <-- DONT USE THIS, USE THE SCRIPT IN PLOTS_ETC

figure(2); clf  
subplot(1,1,1)
for m = 1:size(I0_orig, 1)
    "Point " + m
    % Initial points
    x_I0 = I0_orig(m, 1)/max(x)*(Npoints-1);
    y_I0 = I0_orig(m, 2)/max(y)*(Npoints-1);
    % Final points
    x_rt = postRT(m, 1); % scaling by Npoints already done in for loop
    y_rt = postRT(m, 2);
%     % Solutions for final points
    x_vec = xpaths{m};
    y_vec = ypaths{m};
    hold on
    plot(x_vec,y_vec,'k:','linewidth',0.1); % solsRT
%     plot([x_rt, x_I0], [y_rt, y_I0],  "b", LineStyle=":", LineWidth=0.001) % Lines connecting the two
    if LocalFailure(m) == 1 % bubbles for non-LRF
%         plot(x_I0, y_I0, 'bo', markersize = 12)
        plot(x_rt, y_rt, 'r.', markersize = 15, LineWidth=2)
    end
end
% Plot separatrix
plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1), LineWidth=2)

% Plot before and after
% x_I0 = I0_orig(:, 1)/max(x)*(Npoints-1);
% y_I0 = I0_orig(:, 2)/max(y)*(Npoints-1);
% pre = scatter(x_I0, y_I0, 50, RSI, 'filled')
% hold on
post = scatter(postRT(:, 1), postRT(:, 2), 75, 'filled')
colormap winter
colorbar

% Set plot window
% colorbar
set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
ax = gca
ax.FontSize = 20
xlabel("E", FontSize = 25)
ylabel("T", FontSize = 25)
xlim([1e-4, 10^2.1])
ylim([1e-10, 1e5])
% tit = "Effects of radiation therapy, rho = " + c3_rho_star
tit = "Radiation therapy with Patient-Specific SF_T / SF_E Ratios"
% tit = "SF2_E = SF2_T / " + 0.9;
title(tit, FontSize = 25)
% sgtitle("Radiation therapy with Patient-Specific SF_T Ratios", fontsize=20)
