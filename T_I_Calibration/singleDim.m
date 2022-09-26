% B = panel();
%     
%     B.pack([0.5 0.25 0.25],3)
% 
% B(1,1).select()

sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;    

Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);



% Calculate separatrix
[curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);

% Solve ODE
rhs = @(t,x)([sigma+rho*x(1,:).*x(2,:)./(eta+x(2,:))-mu*x(1,:).*x(2,:)-delta*x(1,:);...
          alpha*x(2,:).*(1-beta*x(2,:))-gamma*x(1,:).*x(2,:)]);

% rhs = @(t,x)([t; ... % varies over time t
%     (alpha*x(2,:).*(1-beta*x(2,:))-gamma*x(1,:).*x(2,:)) / ... % numerator (T-dot)
%     (sigma+rho*x(1,:).*x(2,:)./(eta+x(2,:))-mu*x(1,:).*x(2,:)-delta*x(1,:)) ... % denominator (E-dot)
%     ]);

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

% dx = (x(2)-x(1))/y(2)-y(1);
% dy = 1;

% [X, Y] = meshgrid(x,y);
% G = rhs([],[linspace(0,30,900);reshape(X,1,[]), reshape(Y,1,[])]); % don't know what to change here
% U = reshape(G(1,:),Npoints, Npoints); % or here
% V = reshape(G(2,:),Npoints, Npoints)*dx/dy; 
% N = sqrt(U.^2+V.^2);
% U = U./N; V = V./N;


% figure(1);clf
% subplot(1,1,1)
% Plot separatrix
figure(1);clf
subplot(3,3,1)
plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1), LineWidth=2)
plot((curve.x/max(x)*(Npoints-1))/(curve.x/max(x)*(Npoints-1)), ...
   (curve.y/max(y)*(Npoints-1))/(curve.x/max(x)*(Npoints-1)), ...
   LineWidth=2)

hold on

% Plot vector field

[X1, Y1] = meshgrid(0:Npoints-1,0:Npoints-1);
q = quiver(X1,Y1,U,V); 
q.Color = [0 0 0]; 
q.AutoScaleFactor = 0.5;

% Make the plot a little prettier
axis([0 22 0 30]);
title('Testing', FontSize=15);

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

% for i = 1:size(initCond,1)
%     h=plot(sols{i}.y(1,:)/max(x)*(Npoints-1),sols{i}.y(2,:)/max(y)*(Npoints-1),'k')
%  h = plot((sols{i}.y(1,:)/max(x)*(Npoints-1))/(sols{i}.y(1,:)/max(x)*(Npoints-1)), ...
%      (sols{i}.y(2,:)/max(y)*(Npoints-1))/(sols{i}.y(1,:)/max(x)*(Npoints-1)), ...
%      'k');
%     set(h,'linewidth',1.5)
% end

% plots solution on vector field
for i = 1:size(initCond,1)
 h=plot(sols{i}.y(1,:)/max(x)*(Npoints-1),sols{i}.y(2,:)/max(y)*(Npoints-1),'k')
 set(h,'linewidth',1.5)
end
%%
% Plotting T v time
% B(2,1).select()
subplot(3,3,5)
for i = 1:size(initCond,1)
 h=plot(sols{i}.y(2,:)/max(y)*(Npoints-1))
 set(h,'linewidth',1.5)
 hold on
end
set(gca,'linewidth',1.5,'tickdir','out','fontsize',14); xlabel('time');ylabel('T')
axis([1 55 0 30])

% Plotting E v time
% B(3,1).select()
subplot(3,3,9)
for i = 1:size(initCond,1)
 h=plot(sols{i}.y(1,:)/max(x)*(Npoints-1))
 set(h,'linewidth',1.5)
 hold on
end
set(gca,'linewidth',1.5,'tickdir','out','fontsize',14); xlabel('time');ylabel('E')
axis([1 55 0 50])


% Plotting T/E v time
figure(2);clf
subplot(1,2,1)
s = cell(1,size(initCond, 1));
for i = 1:size(initCond,1)
 s{i} = (sols{i}.y(2,:)/max(y)*(Npoints-1))./(sols{i}.y(1,:)/max(x)*(Npoints-1));
 h=plot(s{i})
 set(h,'linewidth',1.5)
 hold on
end
set(gca,'linewidth',1.5,'tickdir','out','fontsize',14); xlabel('time');ylabel('T/E')
axis([1 25 0 30])

subplot(1,2,2)
% Getting approximate derivative (diff)
d = cell(1, size(initCond, 1))
for i = 1:size(initCond,1)
    d{i} = diff(s{i})
    h=plot(d{i})
    set(h,'linewidth',1.5)
    hold on
end
set(gca,'linewidth',1.5,'tickdir','out','fontsize',14); xlabel('time');ylabel('Approx. d/dt')
axis([1 25 0 30])

