%% Import
% ti_all= readtable("/Users/saketpandit/Documents/Moffitt/Project/Data/TumorAndImmuneCounts_Anonymized.xlsx");
% 
% % Extracting Non-small cell lung cancer
% cond = ti_all.DiseaseSite == "Lung" ...
% & ti_all.HistologicalSubtypes ~= "LUSC"...
% & ti_all.HistologicalSubtypes ~= "";
% % this ends up selecting for just LUAD, lung adenocarcinoma
% 
% % subset of interest
% ti_sub = ti_all(cond,:);

% data for those with failure
ti_all= readtable("/Users/saketpandit/Documents/Moffitt/Project/Data/TumorAndImmuneCounts_Anonymized.xlsx");
ti_outcome = readtable("/Users/saketpandit/Documents/Moffitt/Project/Data/TumorAndImmuneCounts_LungOutcomeSubset_Anonymized.xlsx");

ti_all = join(ti_outcome(1:end-9,:), ti_all);

% % Extracting Non-small cell lung cancer
cond = ti_all.DiseaseSite == "Lung";
% & ti_all.HistologicalSubtypes ~= "LUSC"...
% & ti_all.HistologicalSubtypes ~= "";

% subset of interest
ti_sub = ti_all(cond,:);

%% Percentage sweeps
figure(1);clf
values = 30; % having it at one place so I can tweak everything at once
thresh_mu = 1
thresh_rho = 2 % rho doesn't get down to 0

%% Case 1: Mu sweep

% Setting initCond
AT = ti_sub.TotalAnti_TumorCase1;
T = ti_sub.Total_TumorCells;
initCond = [AT/max(T)*20, T/max(T)*20]; % correction with scalar value


%initCond = [c1_AT/max(c1_AT), c1_T/max(c1_T)*1000]; % correction to get it on a good scale

% sweep params
first = 0.03;
last = 0.25;

% Mu sweep
c1_percs_mu = param_sweep(first, last, values, "mu", initCond);

% finding mu_star (first value that brings percentage to 0)
% percs = c1_percs_mu;
for i = 1:size(c1_percs_mu, 1)
    c1_percs_mu = sortrows(c1_percs_mu,2, 'descend');
    if c1_percs_mu(i,2) < thresh_mu
        c1_mu_star = c1_percs_mu(i,1);
        break
    end
end


% Plotting the percentages
subplot(3,2,1)
% plot(percs(:, 1), percs(:,2)./(percs(:,2)+percs(:,3)), ...
%     LineWidth=3)
plot(c1_percs_mu(:, 1), c1_percs_mu(:,2)./(c1_percs_mu(:,2)+c1_percs_mu(:,3)), ...
    LineWidth=3)
xl = xline(c1_mu_star, "--r", "\mu* = " + c1_mu_star);
xl.LabelHorizontalAlignment = "center";
xl.LabelVerticalAlignment = "middle";
%xlabel("\mu", FontSize = 18)
ylabel("% Tumor resolution", FontSize=18)
title("Case 1 \mu sweep", FontSize=20)
% c1_mu_star = 

%% Case 2: Mu sweep

% Setting initCond
AT = ti_sub.TotalAnti_TumorCase2;
T = ti_sub.Total_TumorCells;
initCond = [AT/max(T)*20, T/max(T)*20]; % correction with scalar value
%initCond = [c1_AT/max(c1_AT), c1_T/max(c1_T)*1000]; % correction to get it on a good scale

% sweep params
% first = 0.006;
first = 0.03;
last = 0.25;

% Mu sweep
c2_percs_mu = param_sweep(first, last, values, "mu", initCond);

% finding mu_star (first value that brings percentage to 0)
% percs = c2_percs_mu;
for i = 1:size(c2_percs_mu, 1)
    c2_percs_mu = sortrows(c2_percs_mu,2, 'descend');
    if c2_percs_mu(i,2) < thresh_mu
        c2_mu_star = c2_percs_mu(i,1);
        break
    end
end


% Plotting the percentages
subplot(3,2,3)
% plot(percs(:, 1), percs(:,2)./(percs(:,2)+percs(:,3)), ...
%     LineWidth=3)
plot(c2_percs_mu(:, 1), c2_percs_mu(:,2)./(c2_percs_mu(:,2)+c2_percs_mu(:,3)), ...
    LineWidth=3)
xl = xline(c2_mu_star, "--r", "\mu* = " + c2_mu_star);
xl.LabelHorizontalAlignment = "center";
xl.LabelVerticalAlignment = "middle";
%xlabel("\mu", FontSize = 18)
ylabel("% Tumor resolution", FontSize=18)
title("Case 2 \mu sweep", FontSize=20)

%% Case 3: Mu sweep
% Setting initCond
AT = ti_sub.TotalAnti_TumorCase3;
T = ti_sub.Total_TumorCells;
initCond = [AT/max(T)*20, T/max(T)*20]; % correction with scalar value
%initCond = [c1_AT/max(c1_AT), c1_T/max(c1_T)*1000]; % correction to get it on a good scale

% sweep params
% first = 0.006;
first = 0.03;
last = 0.25;

% Mu sweep
c3_percs_mu = param_sweep(first, last, values, "mu", initCond);

% finding mu_star (first value that brings percentage to 0)
% percs = c3_percs_mu;
for i = 1:size(c3_percs_mu, 1)
    c3_percs_mu = sortrows(c3_percs_mu,2, 'descend');
    if c3_percs_mu(i,2) < thresh_mu
        c3_mu_star = c3_percs_mu(i,1);
        break
    end
end

% Plotting the percentages
subplot(3,2,5)
% plot(percs(:, 1), percs(:,2)./(percs(:,2)+percs(:,3)), ...
%     LineWidth=3)
plot(c3_percs_mu(:, 1), c3_percs_mu(:,2)./(c3_percs_mu(:,2)+c3_percs_mu(:,3)), ...
    LineWidth=3)
xl = xline(c3_mu_star, "--r", "\mu* = " + c3_mu_star);
xl.LabelHorizontalAlignment = "center";
xl.LabelVerticalAlignment = "middle";
xlabel("\mu", FontSize = 18)
ylabel("% Tumor resolution", FontSize=18)
title("Case 3 \mu sweep", FontSize=20)

%% Case 1: Rho sweep

% Setting initCond
AT = ti_sub.TotalAnti_TumorCase1;
T = ti_sub.Total_TumorCells;
initCond = [AT/max(T)*20, T/max(T)*20]; % correction with scalar value
%initCond = [c1_AT/max(c1_AT), c1_T/max(c1_T)*1000]; % correction to get it on a good scale

% sweep params
first = 0;
last = 0.3;

% rho sweep
c1_percs_rho = param_sweep(first, last, values, "rho", initCond);

% finding rho_star (first value that brings percentage to 0)
% percs = c1_percs_rho;
c1_rho_star = 0;
for i = 1:size(c1_percs_rho, 1)
    c1_percs_rho = sortrows(c1_percs_rho,[2 1], 'descend');
    if c1_percs_rho(i,2) < thresh_rho % varying rho does not allow us to reach 0%
        c1_rho_star = c1_percs_rho(i,1);
        break
    end
end
c1_rho_star

% Plotting the percentages
subplot(3,2,2)
% plot(percs(:, 1), percs(:,2)./(percs(:,2)+percs(:,3)), ...
%     LineWidth=3)
plot(c1_percs_rho(:, 1), c1_percs_rho(:,2)./(c1_percs_rho(:,2)+c1_percs_rho(:,3)), ...
    LineWidth=3)
xl = xline(c1_rho_star, "--r", "\rho* = " + c1_rho_star);
xl.LabelHorizontalAlignment = "center";
xl.LabelVerticalAlignment = "middle";
%xlabel("\rho", FontSize = 18)
%ylabel("% Tumor resolution", FontSize=18)
title("Case 1 \rho sweep", FontSize=20)

%% Case 2: Rho sweep
% Setting initCond
AT = ti_sub.TotalAnti_TumorCase2;
T = ti_sub.Total_TumorCells;
initCond = [AT/max(T)*20, T/max(T)*20]; % correction with scalar value
%initCond = [c1_AT/max(c1_AT), c1_T/max(c1_T)*1000]; % correction to get it on a good scale

% sweep params
first = 0;
% last = 0.3;
last = 0.3;

% rho sweep
c2_percs_rho = param_sweep(first, last, values, "rho", initCond); %sigma = 0.579

% finding rho_star (first value that brings percentage to 0)
% percs = c2_percs_rho;
for i = 1:size(c2_percs_rho, 1)
    c2_percs_rho = sortrows(c2_percs_rho, [2 1], 'descend');
    if c2_percs_rho(i,2) < thresh_rho
        c2_rho_star = c2_percs_rho(i,1);
        break
    end
end

% Plotting the percentages
subplot(3,2,4)
% plot(percs(:, 1), percs(:,2)./(percs(:,2)+percs(:,3)), ...
%     LineWidth=3)
% subplot(1,1,1)
plot(c2_percs_rho(:, 1), c2_percs_rho(:,2)./(c2_percs_rho(:,2)+c2_percs_rho(:,3)), ...
    LineWidth=3)
xl = xline(c2_rho_star, "--r", "\rho* = " + c2_rho_star);
xl.LabelHorizontalAlignment = "center";
xl.LabelVerticalAlignment = "middle";
%xlabel("\rho", FontSize = 18)
%ylabel("% Tumor resolution", FontSize=18)
title("Case 2 \rho sweep", FontSize=20)

%% Case 3: Rho sweep
% Setting initCond
AT = ti_sub.TotalAnti_TumorCase3;
T = ti_sub.Total_TumorCells;
initCond = [AT/max(T)*20, T/max(T)*20]; % correction with scalar value
%initCond = [c1_AT/max(c1_AT), c1_T/max(c1_T)*1000]; % correction to get it on a good scale

% sweep params
first = 0;
% last = 0.3;
last = 0.3;

% rho sweep
c3_percs_rho = param_sweep(first, last, values, "rho", initCond);

% finding rho_star (first value that brings percentage to 0)
% percs = c3_percs_rho;
for i = 1:size(c3_percs_rho, 1)
    c3_percs_rho = sortrows(c3_percs_rho, [2 1], 'descend');
    if c3_percs_rho(i,2) < thresh_rho
        c3_rho_star = c3_percs_rho(i,1);
        break
    end
end

% Plotting the percentages
subplot(3,2,6)
% plot(percs(:, 1), percs(:,2)./(percs(:,2)+percs(:,3)), ...
%     LineWidth=3)
plot(c3_percs_rho(:, 1), c3_percs_rho(:,2)./(c3_percs_rho(:,2)+c3_percs_rho(:,3)), ...
    LineWidth=3)
xl = xline(c3_rho_star, "--r", "\rho* = " + c3_rho_star);
xl.LabelHorizontalAlignment = "center";
xl.LabelVerticalAlignment = "middle";
xlabel("\rho", FontSize = 18)
%ylabel("% Tumor resolution", FontSize=18)
title("Case 3 \rho sweep", FontSize=20)


%% Scatterplots %%

param_stars = [
    c1_mu_star c1_rho_star;
    c2_mu_star c2_rho_star;
    c3_mu_star c3_rho_star;
    ];

%% Case 1, mu
% Setting initCond
AT = ti_sub.TotalAnti_TumorCase1;
T = ti_sub.Total_TumorCells;
initCond = [AT/max(T)*20, T/max(T)*20]; % correction with scalar value

% initializing parameter values
rho = 0.95;  eta = 20.19; 
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;

% Changed from Kuznetzov
sigma = 0.6; % sigma = 0.118 --> 0.6
mu = c1_mu_star; % mu = 0.00311 --> 0.0485


figure(2);clf
subplot(3,2,1)
% Calculate separatrix
[curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);

% Solve ODE
rhs = @(t,x)([sigma+rho*x(1,:).*x(2,:)./(eta+x(2,:))-mu*x(1,:).*x(2,:)-delta*x(1,:);...
          alpha*x(2,:).*(1-beta*x(2,:))-gamma*x(1,:).*x(2,:)]);

options = odeset('Refine',100);
solve = @(init)(ode45(rhs,[0 100],init,options));
Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);
 
dx = x(2)-x(1);
 
dy = y(2)-y(1);
[X, Y] = meshgrid(x,y);
G = rhs([],[reshape(X,1,[]); reshape(Y,1,[])]);
U = reshape(G(1,:),Npoints,Npoints);
V = reshape(G(2,:),Npoints,Npoints)*dx/dy;
N = sqrt(U.^2+V.^2);
U = U./N; V = V./N;

plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1), LineWidth=2)

hold on

% % Plot vector field
% [X1, Y1] = meshgrid(0:Npoints-1,0:Npoints-1);
% q = quiver(X1,Y1,U,V); 
% q.Color = [0 0 0]; 
% q.AutoScaleFactor = 0.5;

% initCond w/ "scaling factor"
scatter(initCond(:,1)./max(x)*(Npoints-1), initCond(:,2)./max(y)*(Npoints-1), ...
        30, 'filled', 'red')


% Make the plot a little prettier
axis([0 22 0 30]);
% set(gca,'xscale', 'log')
% set(gca, 'yscale', 'log')
ttl = "Case 1 Separatrix, \mu = " + c1_mu_star
title(ttl, FontSize=15);

%% Case 2, mu
% Setting initCond
AT = ti_sub.TotalAnti_TumorCase2;
T = ti_sub.Total_TumorCells;
initCond = [AT/max(T)*20, T/max(T)*20]; % correction with scalar value

% initializing parameter values
rho = 0.95;  eta = 20.19; 
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;

% Changed from Kuznetzov
sigma = 0.6; % sigma = 0.118 --> 0.6
mu = c2_mu_star; % mu = 0.00311 --> 0.0485


subplot(3,2,3)
% Calculate separatrix
[curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);

% Solve ODE
rhs = @(t,x)([sigma+rho*x(1,:).*x(2,:)./(eta+x(2,:))-mu*x(1,:).*x(2,:)-delta*x(1,:);...
          alpha*x(2,:).*(1-beta*x(2,:))-gamma*x(1,:).*x(2,:)]);

options = odeset('Refine',100);
solve = @(init)(ode45(rhs,[0 100],init,options));
Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);
 
dx = x(2)-x(1);
 
dy = y(2)-y(1);
[X, Y] = meshgrid(x,y);
G = rhs([],[reshape(X,1,[]); reshape(Y,1,[])]);
U = reshape(G(1,:),Npoints,Npoints);
V = reshape(G(2,:),Npoints,Npoints)*dx/dy;
N = sqrt(U.^2+V.^2);
U = U./N; V = V./N;

plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1), LineWidth=2)

hold on

% % Plot vector field
% [X1, Y1] = meshgrid(0:Npoints-1,0:Npoints-1);
% q = quiver(X1,Y1,U,V); 
% q.Color = [0 0 0]; 
% q.AutoScaleFactor = 0.5;

% initCond w/ "scaling factor"
scatter(initCond(:,1)./max(x)*(Npoints-1), initCond(:,2)./max(y)*(Npoints-1), ...
        30, 'filled', 'red')


% Make the plot a little prettier
axis([0 22 0 30]);
% set(gca,'xscale', 'log')
% set(gca, 'yscale', 'log')
ttl = "Case 2 Separatrix, \mu = " + c2_mu_star
title(ttl, FontSize=15);

%% Case 3, mu
% Setting initCond
AT = ti_sub.TotalAnti_TumorCase3;
T = ti_sub.Total_TumorCells;
initCond = [AT/max(T)*20, T/max(T)*20]; % correction with scalar value

% initializing parameter values
rho = 0.95;  eta = 20.19; 
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;

% Changed from Kuznetzov
sigma = 0.6; % sigma = 0.118 --> 0.6
mu = c3_mu_star; % mu = 0.00311 --> 0.0485


subplot(3,2,5)
% Calculate separatrix
[curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);

% Solve ODE
rhs = @(t,x)([sigma+rho*x(1,:).*x(2,:)./(eta+x(2,:))-mu*x(1,:).*x(2,:)-delta*x(1,:);...
          alpha*x(2,:).*(1-beta*x(2,:))-gamma*x(1,:).*x(2,:)]);

options = odeset('Refine',100);
solve = @(init)(ode45(rhs,[0 100],init,options));
Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);
 
dx = x(2)-x(1);
 
dy = y(2)-y(1);
[X, Y] = meshgrid(x,y);
G = rhs([],[reshape(X,1,[]); reshape(Y,1,[])]);
U = reshape(G(1,:),Npoints,Npoints);
V = reshape(G(2,:),Npoints,Npoints)*dx/dy;
N = sqrt(U.^2+V.^2);
U = U./N; V = V./N;

plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1), LineWidth=2)

hold on

% % Plot vector field
% [X1, Y1] = meshgrid(0:Npoints-1,0:Npoints-1);
% q = quiver(X1,Y1,U,V); 
% q.Color = [0 0 0]; 
% q.AutoScaleFactor = 0.5;

% initCond w/ "scaling factor"
scatter(initCond(:,1)./max(x)*(Npoints-1), initCond(:,2)./max(y)*(Npoints-1), ...
        30, 'filled', 'red')


% Make the plot a little prettier
axis([0 22 0 30]);
% set(gca,'xscale', 'log')
% set(gca, 'yscale', 'log')
ttl = "Case 3 Separatrix, \mu = " + c3_mu_star
title(ttl, FontSize=15);

%% Case 1, rho
% Setting initCond
AT = ti_sub.TotalAnti_TumorCase1;
T = ti_sub.Total_TumorCells;
initCond = [AT/max(T)*20, T/max(T)*20]; % correction with scalar value

% initializing parameter values
eta = 20.19; mu = 0.00311
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;

% Changed from Kuznetzov
sigma = 0.6; % sigma = 0.118 --> 0.6
rho = c1_rho_star; 

subplot(3,2,2)
% Calculate separatrix
[curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);

% Solve ODE
rhs = @(t,x)([sigma+rho*x(1,:).*x(2,:)./(eta+x(2,:))-mu*x(1,:).*x(2,:)-delta*x(1,:);...
          alpha*x(2,:).*(1-beta*x(2,:))-gamma*x(1,:).*x(2,:)]);

options = odeset('Refine',100);
solve = @(init)(ode45(rhs,[0 100],init,options));
Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);
 
dx = x(2)-x(1);
 
dy = y(2)-y(1);
[X, Y] = meshgrid(x,y);
G = rhs([],[reshape(X,1,[]); reshape(Y,1,[])]);
U = reshape(G(1,:),Npoints,Npoints);
V = reshape(G(2,:),Npoints,Npoints)*dx/dy;
N = sqrt(U.^2+V.^2);
U = U./N; V = V./N;

plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1), LineWidth=2)

hold on

% % Plot vector field
% [X1, Y1] = meshgrid(0:Npoints-1,0:Npoints-1);
% q = quiver(X1,Y1,U,V); 
% q.Color = [0 0 0]; 
% q.AutoScaleFactor = 0.5;

% initCond w/ "scaling factor"
scatter(initCond(:,1)./max(x)*(Npoints-1), initCond(:,2)./max(y)*(Npoints-1), ...
        30, 'filled', 'red')


% Make the plot a little prettier
axis([0 22 0 30]);
% set(gca,'xscale', 'log')
% set(gca, 'yscale', 'log')
ttl = "Case 1 Separatrix, \rho = " + c1_rho_star
title(ttl, FontSize=15);

%% Case 2, rho
% Setting initCond
AT = ti_sub.TotalAnti_TumorCase2;
T = ti_sub.Total_TumorCells;
initCond = [AT/max(T)*20, T/max(T)*20]; % correction with scalar value

% initializing parameter values
eta = 20.19; mu = 0.00311
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;

% Changed from Kuznetzov
sigma = 0.6; % sigma = 0.118 --> 0.6
rho = c2_rho_star; 

subplot(3,2,4)
% Calculate separatrix
[curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);

% Solve ODE
rhs = @(t,x)([sigma+rho*x(1,:).*x(2,:)./(eta+x(2,:))-mu*x(1,:).*x(2,:)-delta*x(1,:);...
          alpha*x(2,:).*(1-beta*x(2,:))-gamma*x(1,:).*x(2,:)]);

options = odeset('Refine',100);
solve = @(init)(ode45(rhs,[0 100],init,options));
Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);
 
dx = x(2)-x(1);
 
dy = y(2)-y(1);
[X, Y] = meshgrid(x,y);
G = rhs([],[reshape(X,1,[]); reshape(Y,1,[])]);
U = reshape(G(1,:),Npoints,Npoints);
V = reshape(G(2,:),Npoints,Npoints)*dx/dy;
N = sqrt(U.^2+V.^2);
U = U./N; V = V./N;

plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1), LineWidth=2)

hold on

% % Plot vector field
% [X1, Y1] = meshgrid(0:Npoints-1,0:Npoints-1);
% q = quiver(X1,Y1,U,V); 
% q.Color = [0 0 0]; 
% q.AutoScaleFactor = 0.5;

% initCond w/ "scaling factor"
scatter(initCond(:,1)./max(x)*(Npoints-1), initCond(:,2)./max(y)*(Npoints-1), ...
        30, 'filled', 'red')


% Make the plot a little prettier
axis([0 22 0 30]);
% set(gca,'xscale', 'log')
% set(gca, 'yscale', 'log')
ttl = "Case 2 Separatrix, \rho = " + c2_rho_star
title(ttl, FontSize=15);

%% Case 3, rho
% Setting initCond
AT = ti_sub.TotalAnti_TumorCase3;
T = ti_sub.Total_TumorCells;
initCond = [AT/max(T)*20, T/max(T)*20]; % correction with scalar value

% initializing parameter values
eta = 20.19; mu = 0.00311; sigma = 0.118; rho = 0.95
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;

% Changed from Kuznetzov
sigma = 0.6; % sigma = 0.118 --> 0.6
rho = c3_rho_star;

subplot(3,2,6)
% Calculate separatrix
[curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);

% Solve ODE
rhs = @(t,x)([sigma+rho*x(1,:).*x(2,:)./(eta+x(2,:))-mu*x(1,:).*x(2,:)-delta*x(1,:);...
          alpha*x(2,:).*(1-beta*x(2,:))-gamma*x(1,:).*x(2,:)]);

options = odeset('Refine',100);
solve = @(init)(ode45(rhs,[0 100],init,options));
Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);
 
dx = x(2)-x(1);
 
dy = y(2)-y(1);
[X, Y] = meshgrid(x,y);
G = rhs([],[reshape(X,1,[]); reshape(Y,1,[])]);
U = reshape(G(1,:),Npoints,Npoints);
V = reshape(G(2,:),Npoints,Npoints)*dx/dy;
N = sqrt(U.^2+V.^2);
U = U./N; V = V./N;

plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1), LineWidth=2)

hold on

% % Plot vector field
% [X1, Y1] = meshgrid(0:Npoints-1,0:Npoints-1);
% q = quiver(X1,Y1,U,V); 
% q.Color = [0 0 0]; 
% q.AutoScaleFactor = 0.5;

% initCond w/ "scaling factor"
scatter(initCond(:,1)./max(x)*(Npoints-1), initCond(:,2)./max(y)*(Npoints-1), ...
        30, 'filled', 'red')


% Make the plot a little prettier
axis([0 22 0 30]);
% set(gca,'xscale', 'log')
% set(gca, 'yscale', 'log')
ttl = "Case 3 Separatrix, \rho = " + c3_rho_star
title(ttl, FontSize=15);