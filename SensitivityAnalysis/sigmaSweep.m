% Overall, the separatrix looks to be stable across the values that I've
% selected (1, 1.1, 1.2, 1.3)

sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;    

Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);
 

figure(1);clf

%m_vec = logspace(-3, 1, 5); [1, 1.1667, 1.333, 1.5000]
% pretty stable behavior with several orders of magnitude change

%m_vec = linspace(1, 1.5, 4); 
% all separatrices are continuous

%% m_vec = 1
subplot(2,2,1)
m_vec = 1

sigma = m_vec*sigma ;

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
%axis([0 40 0 100]) % had to change axes for high Npoints
title('\sigma = \sigma x 1', FontSize=15);
% legendString = "\sigma = " + m_vec;
% legend(legendString);
% Weird behavior between x = 0.1 --> x = 30. Spikes up and either stays flat or comes back down 
% Behavior is more pronounced with higher N-values (spikes up higher)

% simulates initial conditions 
initCond = [0.2, 100; % Tumor escape
            0.5, 300; % Tumor escape
            1.0, 80;  % Tumor resolution
            1.2, 450; % Tumor escape
            1.4, 380; % Tumor escape
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

%% m_vec = 1.1667
subplot(2,2,2)
m_vec = 1.1667

sigma = m_vec*sigma ;

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
title('\sigma = \sigma x 1.1667', FontSize=15);

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


%% m_vec = 1.3333
subplot(2,2,3)
m_vec = 1.3333

sigma = m_vec*sigma ;

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
title('\sigma = \sigma x 1.3333', FontSize=15);

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


%% m_vec = 1.5
subplot(2,2,4)
m_vec = 1.5

sigma = m_vec*sigma ;

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
%axis([0 40 0 100]) % had to change axes for high Npoints
title('\sigma = \sigma x 1.5', FontSize=15);


% simulates initial conditions 
initCond = [0.2, 100; % Tumor escape
            0.5, 300; % Tumor escape --> right on the line!!
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

%% Weird separatrix
% Create the figure
% B = panel();
%     
%     B.pack([0.5 0.25 0.25],1)

% Set m_vec + sigma
m_vec = 0.4
sigma = m_vec*sigma ;

% Calculate separatrix
% B(1,1).select() % first plot
subplot(1,1,1)
[curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);

% Solve ODE
rhs = @(t,x)([sigma+rho*x(1,:).*x(2,:)./(eta+x(2,:))-mu*x(1,:).*x(2,:)-delta*x(1,:);...
          alpha*x(2,:).*(1-beta*x(2,:))-gamma*x(1,:).*x(2,:)]);

solve = @(init)(ode45(rhs,[0 100],init,options));
 options = odeset('Refine',100);
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
%axis([0 40 0 100]) % had to change axes for high Npoints
title('\sigma = \sigma x 0.4', FontSize=15);

% simulates initial conditions 
initCond = [0.15, 300;
            0.15, 100;
            0.8, 400;
            2, 250;
            1, 125]       ; 

        
sols = cell(1,size(initCond,1));
for i = 1:size(initCond,1)
sols{i} = solve(initCond(i,:));
end

for i = 1:size(initCond,1)
 h=plot(sols{i}.y(1,:)/max(x)*(Npoints-1),sols{i}.y(2,:)/max(y)*(Npoints-1),'k')
 set(h,'linewidth',1.5)
end

% B(2,1).select()
% for i = 1:size(initCond,1)
%  h=plot(sols{i}.y(2,:)/max(y)*(Npoints-1))
%  set(h,'linewidth',1.5)
%  hold on
% end
% set(gca,'linewidth',1.5,'tickdir','out','fontsize',14); xlabel('time');ylabel('T')
% axis([1 45 0 40])
% 
% B(3,1).select()
% for i = 1:size(initCond,1)
%  h=plot(sols{i}.y(1,:)/max(x)*(Npoints-1))
%  set(h,'linewidth',1.5)
%  hold on
% end
% set(gca,'linewidth',1.5,'tickdir','out','fontsize',14); xlabel('time');ylabel('E')
% axis([1 45 0 20]) 