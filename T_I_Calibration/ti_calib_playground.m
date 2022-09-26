%% Import
ti_all= readtable("/Users/saketpandit/Documents/Moffitt/Project/Data/TumorAndImmuneCounts_Anonymized.xlsx");

%% Extracting Non-small cell lung cancer
cond = ti_all.DiseaseSite == "Lung" ...
& ti_all.HistologicalSubtypes ~= "LUSC"...
& ti_all.HistologicalSubtypes ~= "";

% subset of interest
ti_sub = ti_all(cond,:);
%% Plotting the initial conditions

% Initializing parameters --> Kuznetzov paper
sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;    

Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);
 
% tweaking parameters -- starting point
sigma = 0.6; %%; 
%rho = 3.87;    
%eta = 20.19;  
mu = 0.0413;
%delta = 2.8; 
%alpha = 1.636; 
%beta = 0.002; 
%gamma = 1;    



figure(1);clf
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

% Plot separatrix
plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1), LineWidth=2)
%plot(curve.x/max(c1_T)*40,curve.y/max(c1_T)*40, LineWidth=2)

hold on

% Plot vector field
[X1, Y1] = meshgrid(0:Npoints-1,0:Npoints-1);
q = quiver(X1,Y1,U,V); 
q.Color = [0 0 0]; 
q.AutoScaleFactor = 0.5;

% Make the plot a little prettier
%axis([0 30000 0 500]); %--> zoomed out plot
axis([0 22 0 30]);
% set(gca,'xscale', 'log')
% set(gca, 'yscale', 'log')
title('Real Data', FontSize=15);

% ACTUAL Initial Conditions
c1_AT = ti_sub.TotalAnti_TumorCase1;
c1_T = ti_sub.Total_TumorCells;
initCond = [c1_AT, c1_T];
%initCond = initCond(1, :); % 2366, 5094.8 ==> X = 19604.3, Y = 328.334
%initCond = initCond(1100, :);
%xind = find(initCond == min(initCond(:,1))) % xind = 97
%initCond = initCond(97,:) % --> 138.3, 3855.1 ==> X = 1145.73, Y = 248.44

initCond = [c1_AT/max(c1_T)*40, c1_T/max(c1_T)*40];
%initCond = initCond(1:5,:)

% Solving for the initial conditions
sols = cell(1,size(initCond,1));
for i = 1:size(initCond,1)
sols{i} = solve(initCond(i,:));
end

% % Plotting the initial conditions
% for i = 1:size(initCond,1)
%  h=plot(sols{i}.y(1,:)/max(x)*(Npoints-1),sols{i}.y(2,:)/max(y)*(Npoints-1),'k');
%  set(h,'linewidth',1.5);
% end

for i = 1:size(initCond,1)
 h=plot(sols{i}.y(1,:),sols{i}.y(2,:),'k');
 set(h,'linewidth',1.5);
end


% Percent classification
% r = 0; nr = 0
% for i = 1:size(initCond,1)
%     if sols{i}.y(2,end) < 65 % tumor resolution
%         r = r+1;
%     else
%         nr = nr+1;
%     end
% end

%% Sweeps w/ % classification
% Initializing parameters --> Kuznetzov paper
sigma = 0.118; rho = 0.95;   eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;    

% Values found from param sweep
sigma = 0.6 %%;
% mu = 0.006; %% vary this one


Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);

% Mu sweep
mus = linspace(0.006, 0.01, 10);

% Classifications
percs = zeros(10,3);

for k = 1:size(mus, 2)
    % sweeping mu
    mu = mus(k);
    
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

    % Setting Initial Conditions
    c1_AT = ti_sub.TotalAnti_TumorCase1;
    c1_T = ti_sub.Total_TumorCells;
    initCond = [c1_AT/max(c1_AT), c1_T/max(c1_T)*1000]; % correction to get it on a good scale
    
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

    percs(k, :) = [mu, r, nr];
    "Run " + k
end
percs



%% 
percs = c1_percs_mu
for i = 1:size(percs, 1)
    if percs(i,2) == 0
        mu_star = percs(i,1);
        break
    end
end
mu_star

