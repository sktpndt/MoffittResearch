%% increasing the tumor resolution space (green area)
% Initializing parameter values

sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;    

Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);
 

figure(1);clf


%subplot(2,2,1)
subplot(1,1,1)

% Tweaking Parameters

% rho and alpha

%rho = 1.3 % more recruitment 
%alpha = 1.6 % (slightly) less tumor growth
% only works within certain bounds
% - increasing by too much pushes the separatrix into the negative
% - increasing outside of the bounds creats a discontinuous separatrix

% calculation failure
% rho = 1.4 % --> rho is too high
% alpha = 1.6

% rho = 1.3
% alpha = 1.5 % --> alpha is too low


% discontinuous separatrix
% rho = 1.3
% alpha = 2.8 --> alpha is too high

% rho = 0.9 % --> rho is too low
% alpha = 1.6

% sigma and delta --> not strong enough to push the separatrix
% sigma = 0.118, delta = 0.374
%sigma = 0.45
%delta = 0.65

% rho and mu
% rho = 0.95, mu = 0.00311
%rho = 0.802 % less tumor recruitment (if tumor recruitment is too high, calculation fails)
%mu = 0.00179 % less immune exhaustion
% separatrix is EXTREMELY sensitive to changes in mu. 



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
title('Extreme Separatrices', FontSize=15);

% simulates initial conditions --> just to see if they're behaving properly
initCond = [0.2, 100; 
            0.5, 300; 
            1.0, 80;  
            1.2, 450; 
            1.4, 380; 
            1.8, 300; 
            ];

        
sols = cell(1,size(initCond,1));
for i = 1:size(initCond,1)
sols{i} = solve(initCond(i,:));
end

for i = 1:size(initCond,1)
 h=plot(sols{i}.y(1,:)/max(x)*(Npoints-1),sols{i}.y(2,:)/max(y)*(Npoints-1),'k')
 set(h,'linewidth',1.5)
end

%% increasing the tumor escape space (red area)

% Initializing parameters
sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;    

Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);

subplot(1,1,1);
 

figure(1);clf

% Tweaking parameters
% beta not strong enough
% alpha is weird 
% - (need to increase alpha as mu increases to keep a
%   continuous separatrix)
% gamma makes a weird separatrix 
%   - (initial conditions do not go where they
%     should)

% rho and mu
% rho = 0.95, mu = 0.00311
% need to increase rho as mu increases
% - counterintuitive results
%rho = 3.0 % --> increased tumor recruitment
%mu = 0.03 % --> increased immune exhaustion

% sigma and delta
% sigma = 0.118, delta = 0.374
% waaay too finicky
%   - Each sigma has a certain range of delta that produces a continuous
%   separatrix
%   - unable to get the separatrix down much further
%sigma = 0.7
%delta = 1.0

% delta and rho
% delta = 0.374, rho = 0.95
%delta = 0.001 % --> decreased efflux
%rho = 0.1 % --> decreased tumor recruitment
% pushes the separatrix out to the side

% Calculate separatrix
[curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);

% Solve ODE
rhs = @(t,x)([sigma+rho*x(1,:).*x(2,:)./(eta+x(2,:))-mu*x(1,:).*x(2,:)-delta*x(1,:);...
          alpha*x(2,:).*(1-beta*x(2,:))-x(1,:).*x(2,:)]);

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
title('Extreme Separatrices', FontSize=15);
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

%% Changing all the values

% Initializing parameters
sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;    

Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);

subplot(1,1,1);
 

figure(1);clf


% Tweaking everything
% for i = 1:100
%     % initialize parameters
%     sigma = sigma+0.1
%     if Kuznetsov_SeparatrixCalc
%         list_of_params = params
%     else
%         rho = rho + 0.1
%         if Kuznetsov_SeparatrixCalc
%             list_of_params = params
%         else
%             mu = mu + 0.001
%             if Kuznetsov_SeparatrixCalc
%                 list_of_params = params
%             else
%                 delta = delta + 0.1
%                 if Kuznetsov_SeparatrixCalc
%                     list_of_params = params
%                 else
%                     stop
%                 end


% Calculate separatrix
[curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);

% Solve ODE
rhs = @(t,x)([sigma+rho*x(1,:).*x(2,:)./(eta+x(2,:))-mu*x(1,:).*x(2,:)-delta*x(1,:);...
          alpha*x(2,:).*(1-beta*x(2,:))-x(1,:).*x(2,:)]);

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
title('Extreme Separatrices', FontSize=15);

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