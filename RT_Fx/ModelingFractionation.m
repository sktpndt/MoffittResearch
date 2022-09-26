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

%figure(1); clf

%% modelling radiation + regrowth: Case 1, mu

% Kuznetzov params
Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);
 
    
sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;   


sigma = 0.6
mu = c1_mu_star % from ti_calibration, 0.1969


"Solving the curve"
[curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);

rhs = @(t,x)([sigma+rho*x(1,:).*x(2,:)./(eta+x(2,:))-mu*x(1,:).*x(2,:)-delta*x(1,:);...
              alpha*x(2,:).*(1-beta*x(2,:))-gamma*x(1,:).*x(2,:)]);
  options = odeset('Refine',100);

solve = @(init)(ode45(rhs,[0 100],init,options));


% Setting up I0_orig
AT = ti_sub.TotalAnti_TumorCase1;
T = ti_sub.Total_TumorCells;

% Scaling values
I0_orig = [AT/max(T)*20, T/max(T)*20];
I0 = I0_orig;

% dose vector
dose_vec = ti_sub.FxSize; % dose per fraction

% number of fractions
num_fx_vec = ti_sub.No_Fractions; % number of doses

% Time steps
fx_dt = 2; % time steps for regrowth (one time step = 12 hrs)


% Percentage resolved vs unresolved
r = 0; nr = 0; % number resolved (r) / not resolved (nr)

% big for loop
"Starting Big Loop"
ronr_m1 = zeros(size(I0_orig, 1), 4);
subplot(3,2,1)
for m = 1:size(I0_orig,1)
    "Run " + m
    start_time = 0;
    I0 = I0_orig(m, :); % initial conditions
    ronr_m1(m, 1:2) = I0_orig(m, :);
    ronr_m1(m, 4) = LocalFailure(m);
    dose = dose_vec(m);
    num_fx = num_fx_vec(m);
    prtk = 1-SF_2(m);
%     subplot(1,1,1)
    % plotting the initial conditions
    plot(I0(1)/max(x)*(Npoints-1),I0(2)/max(y)*(Npoints-1),'b.','linewidth',1.5,'markersize',10)
    hold on
    % how many fractions to run
   
    for k = 1:num_fx
        % recalculate initCond
        initCond = I0 .* prtk;
        sols = solve(initCond);
        y_vec = sols.y(2,:)/max(y)*(Npoints-1);
        x_vec = sols.y(1,:)/max(x)*(Npoints-1);

       if k < num_fx % is this the final fraction?
           % updating I0 to treat w/ another dose 
           I0 = [x_vec(fx_dt)*max(x) y_vec(fx_dt)*max(y)]./(Npoints-1);
           start_time = start_time + fx_dt;
           plot(initCond(1)/max(x)*(Npoints-1),initCond(2)/max(y)*(Npoints-1),'b+','linewidth',0.5,'markersize',0.01)
%            plot(x_vec(1:fx_dt),y_vec(1:fx_dt),'k-.','linewidth',0.01); 
       elseif k == num_fx % what to do on final fraction
           % run it out all the way
           plot(x_vec,y_vec,'k','linewidth',0.4); 
           plot(initCond(1)/max(x)*(Npoints-1),initCond(2)/max(y)*(Npoints-1),'r.','linewidth',20,'markersize',12)
           hold on
           start_time = start_time + fx_dt;
           
           % Percent classification
           if y_vec(end) > 10 % y value for resolution
               ronr_m1(m, 3) = 1;
           end
       end
    end
end

prtcorrect_m1 = 0
for m = 1:size(ronr_m1, 1)
    if ronr_m1(m, 3) == ronr_m1(m, 4)
        prtcorrect_m1 = prtcorrect_m1 + 1;
    end
end
"Percent correct = " + prtcorrect_m1/size(ronr_m1, 1)


% rt1_mu_perc = r/(r+nr)

set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
tit = "Effects of radiation therapy, mu = " + c1_mu_star
title(tit)
% "Percentage resolved, c1_mu_star = " + r/(r+nr)


%% Case 2, c2_mu_star

% Kuznetzov params
Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);
 
    
sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;   

sigma = 0.6;
mu = c2_mu_star % from ti_calibration, 0.0603



[curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);

rhs = @(t,x)([sigma+rho*x(1,:).*x(2,:)./(eta+x(2,:))-mu*x(1,:).*x(2,:)-delta*x(1,:);...
              alpha*x(2,:).*(1-beta*x(2,:))-gamma*x(1,:).*x(2,:)]);
  options = odeset('Refine',100);

solve = @(init)(ode45(rhs,[0 100],init,options));


% Setting up I0_orig
AT = ti_sub.TotalAnti_TumorCase2;
T = ti_sub.Total_TumorCells;

% Scaling values
I0_orig = [AT/max(T)*20, T/max(T)*20];
I0 = I0_orig;

% dose vector
dose_vec = ti_sub.FxSize; % dose per fraction

% number of fractions
num_fx_vec = ti_sub.No_Fractions; % number of doses

% Time steps
fx_dt = 2; % time steps for regrowth (one time step = 12 hrs)


% Percentage resolved vs unresolved
r = 0; nr = 0; % number resolved (r) / not resolved (nr)

% big for loop
% figure(2); clf
subplot(3,2,3)
"Starting Big Loop"
ronr_m2 = zeros(size(I0_orig, 1), 4);
for m = 1:size(I0_orig,1)
    "Run " + m
    start_time = 0;
    I0 = I0_orig(m, :); % initial conditions
    ronr_m2(m, 1:2) = I0_orig(m, :);
    ronr_m2(m, 4) = LocalFailure(m);
    dose = dose_vec(m);
    num_fx = num_fx_vec(m);
    prtk = 1-SF_2(m);
%     subplot(1,1,1)
    % plotting the initial conditions
    plot(I0(1)/max(x)*(Npoints-1),I0(2)/max(y)*(Npoints-1),'b.','linewidth',1.5,'markersize',10)
    hold on
    % how many fractions to run
   
    for k = 1:num_fx
        % recalculate initCond
        initCond = I0 .* prtk;
        sols = solve(initCond);
        y_vec = sols.y(2,:)/max(y)*(Npoints-1);
        x_vec = sols.y(1,:)/max(x)*(Npoints-1);

       if k < num_fx % is this the final fraction?
           % updating I0 to treat w/ another dose 
           I0 = [x_vec(fx_dt)*max(x) y_vec(fx_dt)*max(y)]./(Npoints-1);
           start_time = start_time + fx_dt;
           plot(initCond(1)/max(x)*(Npoints-1),initCond(2)/max(y)*(Npoints-1),'b+','linewidth',0.5,'markersize',0.01)
%            plot(x_vec(1:fx_dt),y_vec(1:fx_dt),'k-.','linewidth',0.5); 
       elseif k == num_fx % what to do on final fraction
           % run it out all the way
           plot(x_vec,y_vec,'k','linewidth',0.4); 
           plot(initCond(1)/max(x)*(Npoints-1),initCond(2)/max(y)*(Npoints-1),'r.','linewidth',20,'markersize',12)
           hold on
           start_time = start_time + fx_dt;
           
           % Percent classification
           if y_vec(end) > 10 % y value for resolution
               ronr_m2(m, 3) = 1;
           end
       end
    end
end

prtcorrect_m2 = 0
for m = 1:size(ronr_m2, 1)
    if ronr_m2(m, 3) == ronr_m2(m, 4)
        prtcorrect_m2 = prtcorrect_m2 + 1;
    end
end
"Percent correct = " + prtcorrect_m2/size(ronr_m2, 1)


% rt2_mu_perc = r/(r+nr)
% it2_mu_perc = c2_percs_mu(c2_percs_mu(:,1)==c2_mu_star, 2)/...
%     (c2_percs_mu(c2_percs_mu(:,1)==c2_mu_star, 2)+ ...
%     c2_percs_mu(c2_percs_mu(:,1)==c2_mu_star, 3))


set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
tit = "Effects of radiation therapy, mu = " + c2_mu_star
title(tit)


% "Percentage resolved, c2_mu_star = " + r/(r+nr)

%% Case 3, c3_mu_star
% Kuznetzov params
Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);
 
    
sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;   

sigma = 0.6;
mu = c3_mu_star % from ti_calibration, 0.0831

[curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);

rhs = @(t,x)([sigma+rho*x(1,:).*x(2,:)./(eta+x(2,:))-mu*x(1,:).*x(2,:)-delta*x(1,:);...
              alpha*x(2,:).*(1-beta*x(2,:))-gamma*x(1,:).*x(2,:)]);
  options = odeset('Refine',100);

solve = @(init)(ode45(rhs,[0 100],init,options));


% Setting up I0_orig
AT = ti_sub.TotalAnti_TumorCase3;
T = ti_sub.Total_TumorCells;

% Scaling values
I0_orig = [AT/max(T)*20, T/max(T)*20];
I0 = I0_orig

% dose vector
dose_vec = ti_sub.FxSize; % dose per fraction

% number of fractions
num_fx_vec = ti_sub.No_Fractions; % number of doses

% Time steps
fx_dt = 2; % time steps for regrowth (one time step = 12 hrs)


% Percentage resolved vs unresolved
r = 0; nr = 0; % number resolved (r) / not resolved (nr)

% big for loop
% figure(3); clf
subplot(3,2,5)
"Starting Big Loop"
ronr_m3 = zeros(size(I0_orig, 1), 4);
for m = 1:size(I0_orig,1)
    "Run " + m
    start_time = 0;
    I0 = I0_orig(m, :); % initial conditions
    ronr_m3(m, 1:2) = I0_orig(m, :);
    ronr_m3(m, 4) = LocalFailure(m);
    dose = dose_vec(m);
    num_fx = num_fx_vec(m);
    prtk = 1-SF_2(m);
%     subplot(1,1,1)
    % plotting the initial conditions
    plot(I0(1)/max(x)*(Npoints-1),I0(2)/max(y)*(Npoints-1),'b.','linewidth',1.5,'markersize',10)
    hold on
    % how many fractions to run
   
    for k = 1:num_fx
        % recalculate initCond
        initCond = I0 .* prtk;
        sols = solve(initCond);
        y_vec = sols.y(2,:)/max(y)*(Npoints-1);
        x_vec = sols.y(1,:)/max(x)*(Npoints-1);

       if k < num_fx % is this the final fraction?
           % updating I0 to treat w/ another dose 
           I0 = [x_vec(fx_dt)*max(x) y_vec(fx_dt)*max(y)]./(Npoints-1);
           start_time = start_time + fx_dt;
           plot(initCond(1)/max(x)*(Npoints-1),initCond(2)/max(y)*(Npoints-1),'b+','linewidth',0.5,'markersize',0.01)
%            plot(x_vec(1:fx_dt),y_vec(1:fx_dt),'k-.','linewidth',0.5); 
       elseif k == num_fx % what to do on final fraction
           % run it out all the way
           plot(x_vec,y_vec,'k','linewidth',0.4); 
           plot(initCond(1)/max(x)*(Npoints-1),initCond(2)/max(y)*(Npoints-1),'r.','linewidth',20,'markersize',12)
           hold on
           start_time = start_time + fx_dt;
           
           % Percent classification
           if y_vec(end) > 10 % y value for resolution
               ronr_m3(m, 3) = 1;
           end
       end
    end
end

prtcorrect_m3 = 0
for m = 1:size(ronr_m3, 1)
    if ronr_m3(m, 3) == ronr_m3(m, 4)
        prtcorrect_m3 = prtcorrect_m3 + 1;
    end
end
"Percent correct = " + prtcorrect_m3/size(ronr_m3, 1)

% rt3_mu_perc = r/(r+nr)
% it3_mu_perc = c3_percs_mu(c3_percs_mu(:,1)==c3_mu_star, 2)/...
%     (c3_percs_mu(c3_percs_mu(:,1)==c3_mu_star, 2)+ ...
%     c3_percs_mu(c3_percs_mu(:,1)==c3_mu_star, 3))

set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
tit = "Effects of radiation therapy, mu = " + c3_mu_star
title(tit)

% "Percentage resolved, c3_mu_star = " + r/(r+nr)

%% Case 1, c1_rho_star
% Kuznetzov params
Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);
 
    
sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;   

sigma = 0.579; % changed value from 0.6 bc there were too many ICs not being resolved
rho = c1_rho_star % from ti_calibration, 0.0517


[curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);

rhs = @(t,x)([sigma+rho*x(1,:).*x(2,:)./(eta+x(2,:))-mu*x(1,:).*x(2,:)-delta*x(1,:);...
              alpha*x(2,:).*(1-beta*x(2,:))-gamma*x(1,:).*x(2,:)]);
  options = odeset('Refine',100);

solve = @(init)(ode45(rhs,[0 100],init,options));


% Setting up I0_orig
AT = ti_sub.TotalAnti_TumorCase1;
T = ti_sub.Total_TumorCells;

% Scaling values
I0_orig = [AT/max(T)*20, T/max(T)*20];
I0 = I0_orig

% dose vector
dose_vec = ti_sub.FxSize; % dose per fraction

% number of fractions
num_fx_vec = ti_sub.No_Fractions; % number of doses

% Time steps
fx_dt = 2; % time steps for regrowth (one time step = 12 hrs)


% Percentage resolved vs unresolved
r = 0; nr = 0; % number resolved (r) / not resolved (nr)

% big for loop
% figure(4); clf
subplot(3,2,2)
"Starting Big Loop"
ronr_r1 = zeros(size(I0_orig, 1), 4);
for m = 1:size(I0_orig,1)
    "Run " + m
    start_time = 0;
    I0 = I0_orig(m, :); % initial conditions
    ronr_r1(m, 1:2) = I0_orig(m, :);
    ronr_r1(m, 4) = LocalFailure(m);
    dose = dose_vec(m);
    num_fx = num_fx_vec(m);
    prtk = 1-SF_2(m);
%     subplot(1,1,1)
    % plotting the initial conditions
    plot(I0(1)/max(x)*(Npoints-1),I0(2)/max(y)*(Npoints-1),'b.','linewidth',0.5,'markersize',10)
    hold on
    % how many fractions to run
   
    for k = 1:num_fx
        % recalculate initCond
        initCond = I0 .* prtk;
        sols = solve(initCond);
        y_vec = sols.y(2,:)/max(y)*(Npoints-1);
        x_vec = sols.y(1,:)/max(x)*(Npoints-1);

       if k < num_fx % is this the final fraction?
           % updating I0 to treat w/ another dose 
           I0 = [x_vec(fx_dt)*max(x) y_vec(fx_dt)*max(y)]./(Npoints-1);
           start_time = start_time + fx_dt;
           a = plot(initCond(1)/max(x)*(Npoints-1),initCond(2)/max(y)*(Npoints-1),'b+','linewidth',0.5,'markersize',0.01);
%            b = plot(x_vec(1:fx_dt),y_vec(1:fx_dt),'k-.','linewidth',0.01);
       elseif k == num_fx % what to do on final fraction
           % run it out all the way
           c = plot(x_vec,y_vec,'k','linewidth',0.4); 
           d = plot(initCond(1)/max(x)*(Npoints-1),initCond(2)/max(y)*(Npoints-1),'r.','linewidth',20,'markersize',15);
%            uistack(d, "top")
           hold on
           start_time = start_time + fx_dt;
           
           % Percent classification
           if y_vec(end) > 10 % y value for resolution
               ronr_r1(m, 3) = 1;
           end
       end
    end
end

prtcorrect_r1 = 0
for m = 1:size(ronr_r1, 1)
    if ronr_r1(m, 3) == ronr_r1(m, 4)
        prtcorrect_r1 = prtcorrect_r1 + 1;
    end
end
"Percent correct = " + prtcorrect_r1/size(ronr_r1, 1)

% rt1_rho_perc = r/(r+nr)
% it1_rho_perc = c1_percs_rho(c1_percs_rho(:,1)==c1_rho_star, 2)/...
%     (c1_percs_rho(c1_percs_rho(:,1)==c1_rho_star, 2)+ ...
%     c1_percs_rho(c1_percs_rho(:,1)==c1_rho_star, 3))

set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
tit = "Effects of radiation therapy, rho = " + c1_rho_star
title(tit)


% "Percentage resolved, c1_rho_star = " + r/(r+nr)

%% Case 2, c2_rho_star
% Kuznetzov params
Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);
 
    
sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;   

sigma = 0.6;
rho = c2_rho_star % from ti_calibration, 0.1034


[curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);

rhs = @(t,x)([sigma+rho*x(1,:).*x(2,:)./(eta+x(2,:))-mu*x(1,:).*x(2,:)-delta*x(1,:);...
              alpha*x(2,:).*(1-beta*x(2,:))-gamma*x(1,:).*x(2,:)]);
  options = odeset('Refine',100);

solve = @(init)(ode45(rhs,[0 100],init,options));


% Setting up I0_orig
AT = ti_sub.TotalAnti_TumorCase2;
T = ti_sub.Total_TumorCells;

% Scaling values
I0_orig = [AT/max(T)*20, T/max(T)*20];
I0 = I0_orig

% dose vector
dose_vec = ti_sub.FxSize; % dose per fraction

% number of fractions
num_fx_vec = ti_sub.No_Fractions; % number of doses

% Time steps
fx_dt = 2; % time steps for regrowth (one time step = 12 hrs)


% Percentage resolved vs unresolved
r = 0; nr = 0; % number resolved (r) / not resolved (nr)

% big for loop
% figure(5); clf
subplot(3,2,4)

"Starting Big Loop"
ronr_r2 = zeros(size(I0_orig, 1), 4);
for m = 1:size(I0_orig,1)
    "Run " + m
    start_time = 0;
    I0 = I0_orig(m, :); % initial conditions
    ronr_r2(m, 1:2) = I0_orig(m, :);
    ronr_r2(m, 4) = LocalFailure(m);
    dose = dose_vec(m);
    num_fx = num_fx_vec(m);
    prtk = 1-SF_2(m);
%     subplot(1,1,1)
    % plotting the initial conditions
    plot(I0(1)/max(x)*(Npoints-1),I0(2)/max(y)*(Npoints-1),'b.','linewidth',1.5,'markersize',10)
    hold on
    % how many fractions to run
   
    for k = 1:num_fx
        % recalculate initCond
        initCond = I0 .* prtk;
        sols = solve(initCond);
        y_vec = sols.y(2,:)/max(y)*(Npoints-1);
        x_vec = sols.y(1,:)/max(x)*(Npoints-1);

       if k < num_fx % is this the final fraction?
           % updating I0 to treat w/ another dose 
           I0 = [x_vec(fx_dt)*max(x) y_vec(fx_dt)*max(y)]./(Npoints-1);
           start_time = start_time + fx_dt;
           plot(initCond(1)/max(x)*(Npoints-1),initCond(2)/max(y)*(Npoints-1),'b+','linewidth',0.5,'markersize',0.01)
%            plot(x_vec(1:fx_dt),y_vec(1:fx_dt),'k-.','linewidth',0.5); 
       elseif k == num_fx % what to do on final fraction
           % run it out all the way
           plot(x_vec,y_vec,'k','linewidth',0.4); 
           plot(initCond(1)/max(x)*(Npoints-1),initCond(2)/max(y)*(Npoints-1),'r.','linewidth',20,'markersize',12)
           hold on
           start_time = start_time + fx_dt;
           
           % Percent classification
           if y_vec(end) > 10 % y value for resolution
               ronr_r2(m, 3) = 1;
           end
       end
    end
end

prtcorrect_r2 = 0
for m = 1:size(ronr_r2, 1)
    if ronr_r2(m, 3) == ronr_r2(m, 4)
        prtcorrect_r2 = prtcorrect_r2 + 1;
    end
end

"Percent correct = " + prtcorrect_r2/size(ronr_r2, 1)



% rt2_rho_perc = r/(r+nr)
% it2_rho_perc = c2_percs_rho(c2_percs_rho(:,1)==c2_rho_star, 2)/...
%     (c2_percs_rho(c2_percs_rho(:,1)==c2_rho_star, 2)+ ...
%     c2_percs_rho(c2_percs_rho(:,1)==c2_rho_star, 3))

set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
tit = "Effects of radiation therapy, rho = " + c2_rho_star
title(tit)


% "Percentage resolved, c2_rho_star = " + r/(r+nr)

%% Case 3, c3_rho_star
figure(1);clf

% Kuznetzov params
Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);
 
    
sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;   

% alpha = 10;
% % gamma = 1;
% % beta = 0.0021;
% % delta = 0.6;
% % mu = 0.01;
% sigma = 1.6;
sigma = 0.6
rho = c3_rho_star % from ti_calibration, 0.0931


[curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);

rhs = @(t,x)([sigma+rho*x(1,:).*x(2,:)./(eta+x(2,:))-mu*x(1,:).*x(2,:)-delta*x(1,:);...
              alpha*x(2,:).*(1-beta*x(2,:))-gamma*x(1,:).*x(2,:)]);
  options = odeset('Refine',100);

solve = @(init)(ode45(rhs,[0 100],init,options));


% Setting up I0_orig
AT = ti_sub.TotalAnti_TumorCase3;
T = ti_sub.Total_TumorCells;

% Scaling values
I0_orig = [AT/max(T)*20, T/max(T)*20];

% AT = AT*max(x)/max(AT)
% T = T*max(y)/max(T)
% 
% I0_orig = [AT, T]*10^6;
% 
% 
% I0 = I0_orig

% dose vector
dose_vec = ti_sub.FxSize; % dose per fraction

% number of fractions
num_fx_vec = ti_sub.No_Fractions; % number of doses

% Time steps
fx_dt = 2; % time steps for regrowth (one time step = 12 hrs)


% Percentage resolved vs unresolved
r = 0; nr = 0; % number resolved (r) / not resolved (nr)

% big for loop
% figure(6); clf
% subplot(3,2,6)
subplot(1,1,1)
% dx = x(2)-x(1);
%  
% dy = y(2)-y(1);
% [X, Y] = meshgrid(x,y);
% G = rhs([],[reshape(X,1,[]); reshape(Y,1,[])]);
% U = reshape(G(1,:),Npoints,Npoints);
% V = reshape(G(2,:),Npoints,Npoints)*dx/dy;
% N = sqrt(U.^2+V.^2);
% U = U./N; V = V./N;

% % Plot separatrix
% plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1), LineWidth=2)

hold on

% Plot vector field
% [X1, Y1] = meshgrid(0:Npoints-1,0:Npoints-1);
% q = quiver(X1,Y1,U,V); 
% q.Color = [0 0 0]; 
% q.AutoScaleFactor = 0.5;
% tic
"Starting Big Loop"
ronr_r3 = zeros(size(I0_orig, 1), 4);
for m = 1:size(I0_orig,1)
%for m = 1:10
    "Run " + m
    start_time = 0;
    I0 = I0_orig(m, :); % initial conditions
    ronr_r3(m, 1:2) = I0_orig(m, :);
    ronr_r3(m, 4) = LocalFailure(m);
    dose = dose_vec(m);
    num_fx = num_fx_vec(m);
    prtk = 1-SF_2(m);
%     subplot(1,1,1)
    % plotting the initial conditions
   plot(I0(1)/max(x)*(Npoints-1),I0(2)/max(y)*(Npoints-1),'b.','linewidth',1.5,'markersize',10)

    hold on
    % how many fractions to run
   
    for k = 1:num_fx
        % recalculate initCond
        initCond = I0 .* prtk;
        sols = solve(initCond);
        y_vec = sols.y(2,:)/max(y)*(Npoints-1);
        x_vec = sols.y(1,:)/max(x)*(Npoints-1);

       if k < num_fx % is this the final fraction?
           % updating I0 to treat w/ another dose 
           I0 = [x_vec(fx_dt)*max(x) y_vec(fx_dt)*max(y)]./(Npoints-1);
           start_time = start_time + fx_dt;
           plot(initCond(1)/max(x)*(Npoints-1),initCond(2)/max(y)*(Npoints-1),'b+','linewidth',0.5,'markersize',0.01)
           %plot(x_vec(1:fx_dt),y_vec(1:fx_dt),'k','linewidth',1); 
       elseif k == num_fx % what to do on final fraction
           % run it out all the way
           plot(x_vec,y_vec,'k','linewidth',0.4); 
           plot(initCond(1)/max(x)*(Npoints-1),initCond(2)/max(y)*(Npoints-1),'r.','linewidth',20,'markersize',10)
           hold on
           start_time = start_time + fx_dt;
           
           % Percent classification
           if y_vec(end) > 10 % y value for resolution
               ronr_r3(m, 3) = 1;
           end
       end
    end
end

prtcorrect_r3 = 0
for m = 1:size(ronr_r3, 1)
    if ronr_r3(m, 3) == ronr_r3(m, 4)
        prtcorrect_r3 = prtcorrect_r3 + 1;
    end
end

"Percent correct = " + prtcorrect_r3/size(ronr_r3, 1)

% rt3_rho_perc = r/(r+nr)
% it3_rho_perc = c3_percs_rho(c3_percs_rho(:,1)==c3_rho_star, 2)/...
%     (c3_percs_rho(c3_percs_rho(:,1)==c3_rho_star, 2)+ ...
%     c3_percs_rho(c3_percs_rho(:,1)==c3_rho_star, 3))


set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
tit = "Effects of radiation therapy, rho = " + c3_rho_star
title(tit)


% "Percentage resolved, c3_rho_star = " + r/(r+nr)

% toc



%% Barplot, percentage that resolve
% Case 1, mu
rt1_mu_perc = 1-sum(ronr_m1(:,3))/size(ronr_m1,1);

it1_mu_perc = c1_percs_mu(c1_percs_mu(:,1)==c1_mu_star, 2)/...
    (c1_percs_mu(c1_percs_mu(:,1)==c1_mu_star, 2)+ ...
    c1_percs_mu(c1_percs_mu(:,1)==c1_mu_star, 3));

% Case 2, mu
rt2_mu_perc = 1-sum(ronr_m2(:,3))/size(ronr_m2,1);

it2_mu_perc = c2_percs_mu(c2_percs_mu(:,1)==c2_mu_star, 2)/...
    (c2_percs_mu(c2_percs_mu(:,1)==c2_mu_star, 2)+ ...
    c2_percs_mu(c2_percs_mu(:,1)==c2_mu_star, 3));

% Case 3, mu
rt3_mu_perc = 1-sum(ronr_m3(:,3))/size(ronr_m3,1);

it3_mu_perc = c3_percs_mu(c3_percs_mu(:,1)==c3_mu_star, 2)/...
    (c3_percs_mu(c3_percs_mu(:,1)==c3_mu_star, 2)+ ...
    c3_percs_mu(c3_percs_mu(:,1)==c3_mu_star, 3));


% Case 1, rho
rt1_rho_perc = 1-sum(ronr_r1(:,3))/size(ronr_r1,1);

it1_rho_perc = c1_percs_rho(c1_percs_rho(:,1)==c1_rho_star, 2)/...
    (c1_percs_rho(c1_percs_rho(:,1)==c1_rho_star, 2)+ ...
    c1_percs_rho(c1_percs_rho(:,1)==c1_rho_star, 3));


% Case 2, rho
rt2_rho_perc = 1-sum(ronr_r2(:,3))/size(ronr_r2,1);

it2_rho_perc = c2_percs_rho(c2_percs_rho(:,1)==c2_rho_star, 2)/...
    (c2_percs_rho(c2_percs_rho(:,1)==c2_rho_star, 2)+ ...
    c2_percs_rho(c2_percs_rho(:,1)==c2_rho_star, 3));


% Case 3, rho
rt3_rho_perc = 1-sum(ronr_r3(:,3))/size(ronr_r3,1);

it3_rho_perc = c3_percs_rho(c3_percs_rho(:,1)==c3_rho_star, 2)/...
    (c3_percs_rho(c3_percs_rho(:,1)==c3_rho_star, 2)+ ...
    c3_percs_rho(c3_percs_rho(:,1)==c3_rho_star, 3));

names = categorical({'Case 1','Case 2','Case 3'});
figure(7);clf

subplot(1,2,1)
vals = [it1_mu_perc it2_mu_perc it3_mu_perc; 
    rt1_mu_perc rt2_mu_perc rt3_mu_perc];
b = bar(names, vals)
title("Mu")

subplot(1,2,2)
vals = [it1_rho_perc rt1_rho_perc; 
    it2_rho_perc rt2_rho_perc; 
    it3_rho_perc rt3_rho_perc;]
bar(names, vals)
title("Rho")
sgtitle("Percentage resolution before and after radiation therapy")
legend(["Before" "After"])

%% ROC Plots
figure(8);clf

lw = 4
titsize = 15
axsize = 15

subplot(1,2,1)
[Xr1, Yr1, T, AUCr1] = perfcurve(ronr_r1(:, 4), ronr_r1(:, 3), 0);
plot(Xr1, Yr1, LineWidth=lw)
hold on
% tit = "AUC for Rho 3 = " + AUCr3


subplot(1,2,1)
[Xr2, Yr2, T, AUCr2] = perfcurve(ronr_r2(:, 4), ronr_r2(:, 3), 0);
plot(Xr2, Yr2, LineWidth=lw)
% tit = "AUC for Rho 2 = " + AUCr2


subplot(1,2,1)
[Xr3, Yr3, T, AUCr3] = perfcurve(ronr_r3(:, 4), ronr_r3(:, 3), 0);
plot(Xr3, Yr3, LineWidth=lw)
hold off
tit = "AUC for Rho [1, 2, 3] = [" + AUCr1 + ", " + AUCr2 + ", " + AUCr3 + "]"
title(tit, FontSize=titsize)
xlabel("False positive rate", FontSize=axsize)
ylabel("True positive rate", FontSize=axsize)
legend('r1', 'r2', 'r3')

subplot(1,2,2)
[Xm1, Ym1, T, AUCm1] = perfcurve(ronr_m1(:, 4), ronr_m1(:, 3), 0);
plot(Xm1, Ym1, LineWidth=lw)
hold on

subplot(1,2,2)
[Xm2, Ym2, T, AUCm2] = perfcurve(ronr_m2(:, 4), ronr_m2(:, 3), 0);
plot(Xm2, Ym2, LineWidth=lw)

subplot(1,2,2)
[Xm3, Ym3, T, AUCm3] = perfcurve(ronr_m3(:, 4), ronr_m3(:, 3), 0);
plot(Xm3, Ym3, LineWidth=lw)
hold off
tit = "AUC for Mu [1, 2, 3] = [" + AUCm1 + ", " + AUCm2 + ", " + AUCm3 + "]"
title(tit, FontSize=titsize)
xlabel("False positive rate", FontSize=axsize)
ylabel("True positive rate", FontSize=axsize)
legend('m1', 'm2', 'm3')




