sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;    

Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);
 

figure(1);clf

subplot(2,4,1)

m_vec = [0.7 1 1.4];
for j = 1:length(m_vec)
    
    sigma = m_vec(j)*sigma ;

    [curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);
    
    plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1))
    hold on

end
axis([0 22 0 30])
title('\sigma sweep')


sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;    

subplot(2,4,2)
for j = 1:length(m_vec)
    
    rho = m_vec(j)*rho ;

    [curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);
    
    plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1))
    hold on

end
axis([0 22 0 30])
title('\rho sweep')

sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;    


subplot(2,4,3)
for j = 1:length(m_vec)
    
    eta = m_vec(j)*eta ;

    [curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);
    
    plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1))
    hold on

end
axis([0 22 0 30])
title('\eta sweep')

sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;    


subplot(2,4,4)
for j = 1:length(m_vec)
    
    mu = m_vec(j)*mu ;

    [curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);
    
    plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1))
    hold on

end
axis([0 22 0 30])
title('\mu sweep')

sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;    


subplot(2,4,5)
for j = 1:length(m_vec)
    
    delta = m_vec(j)*delta ;

    [curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);
    
    plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1))
    hold on

end
axis([0 22 0 30])
title('\delta sweep')

sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;    


subplot(2,4,6)
for j = 1:length(m_vec)
    
    alpha = m_vec(j)*alpha ;

    [curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);
    
    plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1))
    hold on

end
axis([0 22 0 30])
title('\alpha sweep')

sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;    


subplot(2,4,7)
for j = 1:length(m_vec)
    
    beta = m_vec(j)*beta ;

    [curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);
    
    plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1))
    hold on

end
axis([0 22 0 30])
title('\beta sweep')

sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;    


subplot(2,4,8)
for j = 1:length(m_vec)
    
    gamma = m_vec(j)*gamma ;

    [curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);
    
    plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1))
    hold on

end
axis([0 22 0 30])
title('\gamma sweep')

%%
figure

sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;    

%sigma = 0; 
%mu = 0.02
%delta = .2;

    [curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);
    
    plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1))
    hold on

axis([0 22 0 30])
%title([])


hold on

sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;    

%sigma = 0; 
%mu = 0.02
delta = 0.03;

    [curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);
    
    plot(curve.x/max(x)*(Npoints-1),curve.y/max(y)*(Npoints-1))
    hold on

axis([0 22 0 30])
%title([])


hold on


