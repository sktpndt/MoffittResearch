sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;    

Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);
 

figure(1);clf
% use m_vec = [0.8000, 0.9333, 1.0667, 1.2000]

%% m_vec = 0.8
subplot(2,2,1)
%subplot(1,1,1)
m_vec = 0.8

mu = m_vec*mu;

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

hold on

% Plot vector field
[X1, Y1] = meshgrid(0:Npoints-1,0:Npoints-1);
q = quiver(X1,Y1,U,V); 
q.Color = [0 0 0]; 
q.AutoScaleFactor = 0.5;

% Make the plot a little prettier
axis([0 22 0 30]);
title('\mu = \mu x 0.8', FontSize=15);

% simulates initial conditions 
initCond = [0.2, 100; % Tumor escape
            0.5, 300; % Tumor escape
            1.0, 80;  % Tumor resolution
            1.2, 450; % Tumor resolution
            1.4, 380; % Tumor resolution
            1.8, 300; % Tumor resolution
            ];

        
sols = cell(1,size(initCond,1));
for i = 1:size(initCond,1)
sols{i} = solve(initCond(i,:));
end

for i = 1:size(initCond,1)
 h=plot(sols{i}.y(1,:)/max(x)*(Npoints-1),sols{i}.y(2,:)/max(y)*(Npoints-1),'k')
 set(h,'linewidth',1.5)
end

%% m_vec = 0.9333
subplot(2,2,2)
%subplot(1,1,1)
m_vec = 0.9333

mu = m_vec*mu;

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

hold on

% Plot vector field
[X1, Y1] = meshgrid(0:Npoints-1,0:Npoints-1);
q = quiver(X1,Y1,U,V); 
q.Color = [0 0 0]; 
q.AutoScaleFactor = 0.5;

% Make the plot a little prettier
axis([0 22 0 30]);
title('\mu = \mu x 0.9333', FontSize=15);

% simulates initial conditions 
initCond = [0.2, 100; % Tumor escape
            0.5, 300; % Tumor escape
            1.0, 80;  % Tumor resolution
            1.2, 450; % Tumor resolution
            1.4, 380; % Tumor resolution
            1.8, 300; % Tumor resolution
            ];

        
sols = cell(1,size(initCond,1));
for i = 1:size(initCond,1)
sols{i} = solve(initCond(i,:));
end

for i = 1:size(initCond,1)
 h=plot(sols{i}.y(1,:)/max(x)*(Npoints-1),sols{i}.y(2,:)/max(y)*(Npoints-1),'k')
 set(h,'linewidth',1.5)
end

%% m_vec = 1.0667
subplot(2,2,3)
%subplot(1,1,1)
m_vec = 1.0667

mu = m_vec*mu;

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

hold on

% Plot vector field
[X1, Y1] = meshgrid(0:Npoints-1,0:Npoints-1);
q = quiver(X1,Y1,U,V); 
q.Color = [0 0 0]; 
q.AutoScaleFactor = 0.5;

% Make the plot a little prettier
axis([0 22 0 30]);
title('\mu = \mu x 1.0667', FontSize=15);

% simulates initial conditions 
initCond = [0.2, 100; % Tumor escape
            0.5, 300; % Tumor escape
            1.0, 80;  % Tumor resolution
            1.2, 450; % Tumor resolution
            1.4, 380; % Tumor resolution
            1.8, 300; % Tumor resolution
            ];

        
sols = cell(1,size(initCond,1));
for i = 1:size(initCond,1)
sols{i} = solve(initCond(i,:));
end

for i = 1:size(initCond,1)
 h=plot(sols{i}.y(1,:)/max(x)*(Npoints-1),sols{i}.y(2,:)/max(y)*(Npoints-1),'k')
 set(h,'linewidth',1.5)
end

%% m_vec = 1.2000
subplot(2,2,4)
%subplot(1,1,1)
m_vec = 1.200

mu = m_vec*mu;

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

hold on

% Plot vector field
[X1, Y1] = meshgrid(0:Npoints-1,0:Npoints-1);
q = quiver(X1,Y1,U,V); 
q.Color = [0 0 0]; 
q.AutoScaleFactor = 0.5;

% Make the plot a little prettier
axis([0 22 0 30]);
title('\mu = \mu x 1.200', FontSize=15);

% simulates initial conditions 
initCond = [0.2, 100; % Tumor escape
            0.5, 300; % Tumor escape
            1.0, 80;  % Tumor resolution
            1.2, 450; % Tumor escape
            1.4, 380; % Tumor resolution
            1.8, 300; % Tumor resolution
            ];

        
sols = cell(1,size(initCond,1));
for i = 1:size(initCond,1)
sols{i} = solve(initCond(i,:));
end

for i = 1:size(initCond,1)
 h=plot(sols{i}.y(1,:)/max(x)*(Npoints-1),sols{i}.y(2,:)/max(y)*(Npoints-1),'k')
 set(h,'linewidth',1.5)
end