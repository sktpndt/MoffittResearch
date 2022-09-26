function [curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma)



    
%% CALCULATING BIGGEST VALUE OF STEADY STATE FOR Y
Smax = max(roots([mu*beta; -mu+beta*(eta*mu+delta-rho); ...
sigma/alpha+rho-eta*mu-delta+beta*delta*eta; ...
eta*(sigma/alpha-delta)]));

 
%% DEFINING THE MODEL (INLINE FUNCTION)
rhs = @(t,x)([sigma+rho*x(1,:).*x(2,:)./(eta+x(2,:))-mu*x(1,:).*x(2,:)-delta*x(1,:);...
alpha*x(2,:).*(1-beta*x(2,:))-gamma.*x(1,:).*x(2,:)]);
 
%% FUNCTION RETURNING MODEL SOLUTION ON [0,100] FOR GIVEN INITIAL CONDITION
solve = @(init)(ode45(rhs,[0 100],init));
 
%% CALCULATING LOWER PART OF THE SEPARATING CURVE
yInit = 50; %initial value
yTol = 1e-5; %tolerance
dy = 5; %initial step size
 
while dy > yTol
sol = solve([0 yInit]);
if abs(sol.y(2,end)-Smax) < 30 %reached maximal steady state 
    yInitPrev = yInit; yInit = yInit-dy; 
else
    dy = dy/2; yInit = yInitPrev; 
end
end

%% CALCULATING UPPER PART OF THE SEPARATING CURVE 
xInit = 0.5; %initial value 
xTol = 1e-5; %tolerance 
dx = 0.5; %initial step size 
while dx > xTol
sol = solve([xInit Smax+100]);
if abs(sol.y(2,end)-Smax) < 30 
    xInitPrev = xInit; xInit = xInit+dx; 
else
    dx = dx/2; xInit = xInitPrev; 
end
end
%% CALCULATING THE FINAL CURVE 
sol2 = solve([0 yInit]); 
sol1 = solve([xInit Smax+100]); 
sol1.y = sol1.y(:,1:find(diff(sol1.y(2,:))>0,1,'first'));
sol2.y = sol2.y(:,1:find(diff(sol2.y(1,:))<0,1,'first'));
 
curve.x = [sol2.y(1,:) sol1.y(1,end:-1:1) ];
curve.y = [sol2.y(2,:) sol1.y(2,end:-1:1) ];
 


end

