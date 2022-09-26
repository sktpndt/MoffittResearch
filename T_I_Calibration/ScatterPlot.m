% Plot separatrix
%% Import
ti_all= readtable("/Users/saketpandit/Documents/Moffitt/Project/Data/TumorAndImmuneCounts_Anonymized.xlsx");
ti_outcome = readtable("/Users/saketpandit/Documents/Moffitt/Project/Data/TumorAndImmuneCounts_LungOutcomeSubset_Anonymized.xlsx");

ti_all = join(ti_outcome(1:end-9,:), ti_all);

% Extracting Non-small cell lung cancer
cond = ti_all.DiseaseSite == "Lung" ...
& ti_all.HistologicalSubtypes ~= "LUSC"...
& ti_all.HistologicalSubtypes ~= "";

% subset of interest
ti_sub = ti_all(cond,:);
%% initCond
% Setting initCond
AT = ti_sub.TotalAnti_TumorCase3;
T = ti_sub.Total_TumorCells;
LocalFailure = ti_sub.LocalFailure;
% initCond = [AT/max(T)*20, T/max(T)*20];

% Scaling values
initCond = [AT/max(T)*20, T/max(T)*20];

% 
% initCond.Properties.VariableNames = ["AT", "T", "LocalFailure"];
% initCond = initCond(100:200,:);

%% Scatter plot, mu*


% initializing parameter values
rho = 0.95;  eta = 20.19; 
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;

% Changed from Kuznetzov
sigma = 0.6; % sigma = 0.118 --> 0.6
% mu = c1_mu_star;  
mu = 0.0424


figure(2);clf
% subplot(3,2,1)
subplot(1,1,1)
% Calculate separatrix
[curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);

% Solve ODE
rhs = @(t,x)([sigma+rho*x(1,:).*x(2,:)./(eta+x(2,:))-mu*x(1,:).*x(2,:)-delta*x(1,:);...
          alpha*x(2,:).*(1-beta*x(2,:))-gamma*x(1,:).*x(2,:)]);

options = odeset('Refine',100);
solve = @(init)(ode45(rhs,[0 100],init,options));
Npoints = 30;
x = linspace(0,3.5,Npoints);
% x = linspace(0,500,Npoints);
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

% Plot vector field
[X1, Y1] = meshgrid(0:Npoints-1, 0:Npoints-1);
q = quiver(X1,Y1,U,V); 
q.Color = [0 0 0]; 
q.AutoScaleFactor = 0.3;

% % initCond w/o "scaling factor"
% scatter(initCond(:,1), initCond(:,2), ...
%         30, 'filled', 'red')

% initCond w/ "scaling factor"
scTbl = table(initCond(:,1)./max(x)*(Npoints-1), ...
    initCond(:,2)./max(y)*(Npoints-1),...
    LocalFailure);
scTbl.Properties.VariableNames = ["AT", "T", "LocalFailure"];
s = scatter(scTbl, "AT", "T", "filled", ColorVariable="LocalFailure")
colormap winter
colorbar
% s.MarkerEdgeColor = "flat"
% s.MarkerFaceColor = ["blue"]
% s = scatter(initCond.Variables, 'AT'./max(x)*(Npoints-1), "T"./max(y)*(Npoints-1), ...
%         30, 'filled', 'ColorVariable','LocalFailure')


% getting initCond solutions
sols = cell(1,size(initCond,1));
for i = 1:size(initCond,1)
    sols{i} = solve(initCond(i,1:2));
end


% % Plotting the soluntions w/o "scaling factor"
% for i = 1:size(initCond,1)
%  h=plot(sols{i}.y(1,:),sols{i}.y(2,:),'k');
%  set(h,'linewidth',1.5);
% end

% plotting the solutions w/ "scaling factor"
for i = 1:size(initCond,1)
 h=plot(sols{i}.y(1,:)/max(x)*(Npoints-1),sols{i}.y(2,:)/max(y)*(Npoints-1),'k');
 set(h,'linewidth',1.5);
end



% Make the plot a little prettier
%axis([0 30000 0 500]); %--> zoomed out plot
axis([0 22 0 30]);
% set(gca,'xscale', 'log')
set(gca, 'yscale', 'log')
% ttl = "Separatrix Calibration, \mu = " + c1_mu_star
ttl = "Separatrix Calibration, \mu = " + mu
title(ttl, FontSize=15);


% PROBLEM %

% * without scaling factor *
% Additionally, the "scaling factor" I pointed out seems to allow the 
% initial conditions to follow the proper trajectory. If I remove this,
% then the initial conditions don't follow the vector field. 

% * with scaling factor *
% While the lines are all pointing upwards, the scatterplot shows that
% these values are clearly below the separatrix. I think the separatrix 
% calculation breaks down when the values are this far out.
% The scaling factor seems to be necessary, but including it pushes the
% values too far out.


%% Scatter plot, rho*

% initializing parameter values
eta = 20.19; mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;

% Changed from Kuznetzov
sigma = 0.6; % sigma = 0.118 --> 0.6
% rho = c1_rho_star; 
rho = 0.066667
% rho = 0.3;


% subplot(3,2,2)
subplot(1,1,1)
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

% % initCond w/o "scaling factor"
% scatter(initCond(:,1), initCond(:,2), ...
%         30, 'filled', 'red')

% initCond w/ "scaling factor"
scatter(initCond(:,1)./max(x)*(Npoints-1), initCond(:,2)./max(y)*(Npoints-1), ...
        30, 'filled', 'red')



% getting initCond solutions
sols = cell(1,size(initCond,1));
for i = 1:size(initCond,1)
    sols{i} = solve(initCond(i,:));
end


% % Plotting the soluntions w/o "scaling factor"
% for i = 1:size(initCond,1)
%  h=plot(sols{i}.y(1,:),sols{i}.y(2,:),'k');
%  set(h,'linewidth',1.5);
% end

% plotting the solutions w/ "scaling factor"
for i = 1:size(initCond,1)
    if sols{i}.y(2,end) > 10
%         i
        h=plot(sols{i}.y(1,:)/max(x)*(Npoints-1),sols{i}.y(2,:)/max(y)*(Npoints-1),'k');
        set(h,'linewidth',1.5);
    end
end



% Make the plot a little prettier
%axis([0 30000 0 500]); %--> zoomed out plot
axis([0 22 0 30]);
% set(gca,'xscale', 'log')
set(gca, 'yscale', 'log')
% ttl = "Separatrix Calibration, \rho = " + c1_rho_star
ttl = "Separatrix Calibration, \rho = " + rho
title(ttl, FontSize=15);

