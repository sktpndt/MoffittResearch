function [percs] = param_sweep(first, last, values, param, initCond)

% Initializing sweep values
params = linspace(first, last, values);

% Percentages
percs = zeros(values,3);

% Initializing parameter values
%sigma = 0.118; ---> changed from Kuznetzov
rho = 0.95;  eta = 20.19; mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;

for k = 1:size(params, 2)
    
    % sweep parameter
    switch param
        case "mu"
            mu = params(k);
            sigma = 0.6;
        case "rho"
            rho = params(k);
            sigma = 0.579;
    end
    
    % Solve ODE
    rhs = @(t,x)([sigma+rho*x(1,:).*x(2,:)./(eta+x(2,:))-mu*x(1,:).*x(2,:)-delta*x(1,:);...
          alpha*x(2,:).*(1-beta*x(2,:))-gamma*x(1,:).*x(2,:)]);

    options = odeset('Refine',100);
    solve = @(init)(ode45(rhs,[0 100],init,options));

%     rhs = @(t,x)([sigma+rho*x(1,:).*x(2,:)./(eta+x(2,:))-mu*x(1,:).*x(2,:)-delta*x(1,:);...
%               alpha*x(2,:).*(1-beta*x(2,:))-gamma*x(1,:).*x(2,:)]);
%     options = odeset('Refine',100);
%     solve = @(init)(ode45(rhs,[0 100],init,options));


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

    switch param
        case "mu"
            percs(k, :) = [mu, r, nr];
        case "rho"
            percs(k, :) = [rho, r, nr];
    end
    "Run " + k
end

end
