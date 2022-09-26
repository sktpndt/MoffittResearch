%% Import
ti_all= readtable("/Users/saketpandit/Documents/Moffitt/Project/Data/TumorAndImmuneCounts_Anonymized.xlsx");

% Extracting Non-small cell lung cancer
cond = ti_all.DiseaseSite == "Lung" ...
& ti_all.HistologicalSubtypes ~= "LUSC"...
& ti_all.HistologicalSubtypes ~= "";

% subset of interest
ti_sub = ti_all(cond,:);
%%%
% Setting initCond
c3_AT = ti_sub.TotalAnti_TumorCase3;
c3_T = ti_sub.Total_TumorCells;
initCond = [c3_AT/max(c3_T)*20, c3_T/max(c3_T)*20]; % correction with scalar value
%initCond = [c1_AT/max(c1_AT), c1_T/max(c1_T)*1000]; % correction to get it on a good scale


%% Mu sweep w/ % classification
% Initializing parameters --> Kuznetzov paper
%sigma = 0.118; ---> changed from paper

sigma = 0.6; rho = 0.95;   eta = 20.19; 
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;    

% Values found from param sweep
sigma = 0.6; %%; allows me to sweep over many values of mu
% mu = 0.006; %% vary this one


Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);

% Mu sweep
n_vals = 30;
param_start = 0.006;
param_end = 0.07;
params = linspace(param_start, param_end, n_vals);

% Classifications
percs_mu = zeros(10,3);

for k = 1:size(params, 2)
    
    % sweeping mu
    mu = params(k);
    
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
    
    % Solving for the initial conditions
    sols = cell(1,size(initCond,1));
    for i = 1:size(initCond,1)
    sols{i} = solve(initCond(i,:));
    end
    
    % Percent classification
    r = 0; nr = 0;
    for j = 1:size(initCond,1)
        if sols{j}.y(2,end) < 65 % y value for resolution
            r = r+1; % number resolved
        else
            nr = nr+1; % number not resolved
        end
    end

    percs_mu(k, :) = [mu, r, nr];
    "Run " + k
end

% plot
figure(1);clf
percs = percs_mu
plot(percs(:, 1), percs(:,2)./(percs(:,2)+percs(:,3)), ...
    LineWidth=3)
xlabel("\mu", FontSize = 18)
ylabel("% Tumor resolution", FontSize=18)
title("\mu sweep", FontSize=20)


%% rho sweep
% Initializing parameters --> Kuznetzov paper
%sigma = 0.118 --> changed from paper

sigma = 0.6; eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;    

% Values found from param sweep
sigma = 0.6; %%; allows me to sweep over many values of rho
% rho = 0; %% vary this one --> 0.4


Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);

% Rho sweep
n_vals = 30;
param_start = 0;
param_end = 0.3;
params = linspace(param_start, param_end, n_vals);

% Classifications
percs_rho = zeros(10,3);

for k = 1:size(params, 2)
    
    % sweeping rho
    rho = params(k);
    
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
    
    % Solving for the initial conditions
    sols = cell(1,size(initCond,1));
    for i = 1:size(initCond,1)
    sols{i} = solve(initCond(i,:));
    end
    
    % Percent classification
    r = 0; nr = 0;
    for j = 1:size(initCond,1)
        if sols{j}.y(2,end) < 65 % y value for resolution
            r = r+1; % number resolved
        else
            nr = nr+1; % number not resolved
        end
    end

    % storing the results
    percs_rho(k, :) = [rho, r, nr];
    "Run " + k
end

% plot
figure(2);clf
percs = percs_rho;
plot(percs(:, 1), percs(:,2)./(percs(:,2)+percs(:,3)), ...
    LineWidth=3)
xlabel("\rho", FontSize = 18)
ylabel("% Tumor resolution", FontSize=18)
title("\rho sweep", FontSize=20)

