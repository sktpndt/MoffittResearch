%% Plotting for sigma
sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;    

Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);
 

%figure(1);clf

subplot(1,1,1)

%m_vec = logspace(-3, 1, 5);
% pretty stable behavior with several orders of magnitude change

m_vec = linspace(1, 1.5, 4); % all separatrices are continuous
% m_vec < 1 --> discontinuous separatrix
% m_vec > 1.5 --> calculation failure


for j = 1:length(m_vec)
    
    sigma = m_vec(j)*sigma ;

    [curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);
    
    plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1), LineWidth=2)
    hold on

end
axis([0 22 0 30]);
%axis([0 40 0 100]) % had to change axes for high Npoints
title('\sigma sweep', FontSize=15);
legendString = "m\_vec = " + m_vec;
legend(legendString);
% Weird behavior between x = 0.1 --> x = 30. Spikes up and either stays flat or comes back down 
% Behavior is more pronounced with higher N-values (spikes up higher)


%% Plotting for rho
sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;    

Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);

subplot(1, 1, 1)

%m_vec = logspace(-3, 1, 5); % logspace is breaking it
% m_vec = [0.7, 1, 1.4]; % works fine
% m_vec = 1.4 % broken

%m_vec = linspace(0.8, 1.3, 5);
% rho is highly sensitive to changes
% In fact, some values "break" the Kuznetzov calculation (<0.8, >1.3)*
% * this specific range only applies if I put it into linspace
% - For example, if m_vec = [0.7, 1, 1.4], runs fine. But, if m_vec = 1.4,
%   then it breaks. 

%m_vec = linspace(1, 1.2, 3); % use these values, they give all continuous separatrices + no "beep"
% m_vec < 1 --> discontinuous separatrix
% m_vec > 1.3 --> calculation failure

for j = 1:length(m_vec)
    
    rho = m_vec(j)*rho ;

    [curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);
    
    plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1))
    hold on

end
axis([0 22 0 30])
title('\rho sweep', FontSize=15)
legendString = "m\_vec = " + m_vec
legend(legendString)

%% Plotting for Eta
sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;    

Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);

% Spanning 5 orders of magnitude
%m_vec = logspace(-3, 1, 5); 
m_vec = logspace(-2, 1, 4); % plot only has 4 spaces ;)

subplot(1,1,1)
for j = 1:length(m_vec)
    
    eta = m_vec(j)*eta ;

    [curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);
    
    plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1))
    hold on

end
axis([0 22 0 30])
title('\eta sweep')
legendString = "m\_vec = " + m_vec
legend(legendString)

% Eta looks to be remarkably stable, over several orders of magnitude

%% Plotting for Mu

sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;    

Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);


%m_vec = logspace(-3, 1, 5); % broken!
%m_vec = [0.7, 1, 1.4];
%m_vec = linspace(0.8, 1.4, 7)
% Mu is sensitive to small changes. 
% - From 0.8*mu --> ~1.2*mu, separatrix is a smooth curve
% - after ~1.2*mu, separatrix gets a spike/plateau behavior
%       * It's difficult to find the point that this occurs, I don't know
%         how many decimal places the kuznetzov calculation can handle
%         (doesn't seem to be more than 4)
% - Kuznetzov calculation breaks for m_vec < 0.8 (get the beep)

m_vec = linspace(0.8, 1.2, 4); % all separatrices are continuous

subplot(1,1,1)
for j = 1:length(m_vec)
    
    mu = m_vec(j)*mu ;

    [curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);
    
    plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1))
    hold on

end
axis([0 22 0 30])
title('\mu sweep')
legendString = "m\_vec = " + m_vec
legend(legendString)


%% Plotting for delta
sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;    


Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);


%m_vec = logspace(-3, 1, 5); %broken!
%m_vec = linspace(0.3, 1.4, 8); 
% delta is sensitive to changes
% - increases in m_vec produce a downward shift
% - for m_vec > 2, separatrix becomes discontinuous
% - for m_vec < 0.1, calculation fails

%m_vec = linspace(0.1, 2, 4); % all separatrices are continuous
m_vec = 0.01

subplot(1,1,1)
for j = 1:length(m_vec)
    
    delta = m_vec(j)*delta ;

    [curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);
    
    plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1))
    hold on

end
axis([0 8 0 30])
title('\delta sweep')
legendString = "m\_vec = " + m_vec
legend(legendString)

%% Plotting for alpha

sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;    

Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);


%m_vec = logspace(-3, 1, 5); %broken!
%m_vec = linspace(0.9, 1.4, 6)
% alpha is not very sensitive to changes. 
%   - Separatrix becomes discontinuous for alpha>1.3*alpha
%   - separatrix calculation fails for alpha < 0.7*alpha
%   - *something* fails when length(m_vec)>7

m_vec = linspace(0.7, 1.3, 4);

subplot(1,1,1)
for j = 1:length(m_vec)
    
    alpha = m_vec(j)*alpha ;

    [curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);
    
    plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1))
    hold on

end
axis([0 22 0 30])
title('\alpha sweep')
legendString = "m\_vec = " + m_vec
legend(legendString)

%% Plotting for beta

sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;    

Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);

%m_vec = logspace(-3, 1, 5); %broken!
% m_vec = linspace(0.1, 2.3, 5) --> discontinuous separatrix
% - m_vec > 1.3 causes a calculation failure (annoying "beep")
% - m_vec < 0.9 produces a discontinuous separatrix

m_vec = linspace(0.9, 1.3, 4);

subplot(1,1,1)
for j = 1:length(m_vec)
    
    beta = m_vec(j)*beta ;

    [curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);
    
    plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1))
    hold on

end
axis([0 22 0 30])
title('\beta sweep')
legendString = "m\_vec = " + m_vec
legend(legendString)

%% Plotting for gamma

sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;    

Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);

%m_vec = logspace(-3, 1, 5); %broken!
%m_vec = [0.7, 1, 1.4]
%m_vec = linspace(0.5, 1.2, 5); --> discontinuous separatrices
m_vec = linspace(1, 1.3, 3)
% m_vec < 1 produces a discontinuous separatrix
% m_vec > 1.3 leads to calculation failure
%   - more than 4 values in the linspace for a max of 1.3 leads to failure
%   - 3 values works fine for a max of 1.3
%   - 4 values works fine for max of 1.2

subplot(1,1,1)
for j = 1:length(m_vec)
    
    gamma = m_vec(j)*gamma ;

    [curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);
    
    plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1))
    hold on

end
axis([0 30 0 30])
title('\gamma sweep')
legendString = "m_vec = " + m_vec
legend(legendString)


%%
% figure
% 
% sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
% delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;    
% 
% %sigma = 0; 
% %mu = 0.02
% %delta = .2;
% 
%     [curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);
%     
%     plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1))
%     hold on
% 
% axis([0 22 0 30])
% %title([])
% 
% 
% hold on
% 
% sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
% delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;    
% 
% %sigma = 0; 
% %mu = 0.02
% delta = 0.03;
% 
%     [curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);
%     
%     plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1))
%     hold on
% 
% axis([0 22 0 30])
% %title([])
% 
% 
% hold on
% 
% 
