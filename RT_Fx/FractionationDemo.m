%% loading data, setting parameters
load('DoseResponse_PDTChemo_params.mat')
Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);
 
    
sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;   

[curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);

rhs = @(t,x)([sigma+rho*x(1,:).*x(2,:)./(eta+x(2,:))-mu*x(1,:).*x(2,:)-delta*x(1,:);...
              alpha*x(2,:).*(1-beta*x(2,:))-gamma*x(1,:).*x(2,:)]);
  options = odeset('Refine',100);

solve = @(init)(ode45(rhs,[0 100],init,options));

% Specific to PDT therapy. Linear quadratic for PDT
xrange1 = 0:60; xrange2 = 0:100;
fun = @(x,xdata)x(1)./(1+exp(x(2).*(xdata-x(3))));
fun2 = @(x,xdata)x(1).*exp((-xdata.*x(2)));


%% setting up fractionation parameters & ICs

% initCond
I0 = [0.5 450];
I0_orig = I0;

fx_dt = 2; % time steps for regrowth
% one time step = 12 hrs
    
%        dose_vec = [15 3 70 7]; % dose per fraction
dose_vec = 15
num_fx_vec = 1

%        num_fx_vec = [1 5 1 10]; % number of doses
       

%num_sims = length(dose_vec);


figure(1); clf

B = panel();
B.pack(3,2)

% panel parameters
pos_array = [1,1;
             1,2;
             2,1;
             2,2];
%% testing different fractionation schedules
         
         
for j = 1:length(dose_vec)
 
   I0 = I0_orig; % initial conditions 
start_time = 0;

   dose = dose_vec(j);
   num_fx = num_fx_vec(j);

   % panel stuff
   subB = B(pos_array(j,1), pos_array(j,2));
   subB.pack(1,2)
   subB(1,1).select();
        % plotting the initial conditions
        plot(I0_orig(1)/max(x)*(Npoints-1),I0_orig(2)/max(y)*(Npoints-1),'b.','linewidth',1.5,'markersize',12)
        hold on
   
   subC=subB(1,2);
   subC.pack(2,1)
   subC(2,1).select();
        % plotting subplots (Evt, Tvt)
        plot(start_time,I0_orig(1)/max(x)*(Npoints-1),'b.','linewidth',1.5,'markersize',12)
        hold on
   subC(1,1).select();
        plot(start_time,I0_orig(2)/max(y)*(Npoints-1),'b.','linewidth',1.5,'markersize',12)
        hold on

   % how many fractions to run
   for k = 1:num_fx
       %start_time = (k-1)*fx_dt+1;
       % recalculate initCond
       initCond = I0.*[fun(p_Tcell_CetBPD,dose) fun(p_Tumor_CetBPD,dose)];

       % running through ODE w/ new initCond
       sols = solve(initCond);
       y_vec = sols.y(2,:)/max(y)*(Npoints-1);
       x_vec = sols.y(1,:)/max(x)*(Npoints-1);
       
       if k < num_fx % is this the final fraction
   
    % plotting for however many timesteps        
   subB(1,1).select();
            plot(initCond(1)/max(x)*(Npoints-1),initCond(2)/max(y)*(Npoints-1),'r.','linewidth',1.5,'markersize',12)
            plot(x_vec(1:fx_dt),y_vec(1:fx_dt),'k','linewidth',1.5); 

   subC(2,1).select();
            plot(start_time,initCond(1)/max(x)*(Npoints-1),'r.','linewidth',1.5,'markersize',12)
            plot(start_time:fx_dt:start_time+fx_dt,x_vec(1:fx_dt),'k','linewidth',1.5); 

   subC(1,1).select();
            plot(start_time,initCond(2)/max(y)*(Npoints-1),'r.','linewidth',1.5,'markersize',12)
            plot(start_time:fx_dt:start_time+fx_dt,y_vec(1:fx_dt),'k','linewidth',1.5); 
   start_time = start_time + fx_dt; % updating start_time

           
       elseif k == num_fx % what to do on final fraction
   % run it out all the way        
   subB(1,1).select();
              plot(x_vec,y_vec,'k','linewidth',1.5); 
              plot(initCond(1)/max(x)*(Npoints-1),initCond(2)/max(y)*(Npoints-1),'r.','linewidth',1.5,'markersize',12)

   subC(2,1).select();
              plot(start_time,initCond(1)/max(x)*(Npoints-1),'r.','linewidth',1.5,'markersize',12)
              plot(start_time+(0:length(x_vec)-1),x_vec,'k','linewidth',1.5); 

   subC(1,1).select();
              plot(start_time,initCond(2)/max(y)*(Npoints-1),'r.','linewidth',1.5,'markersize',12)
              plot(start_time+(0:length(y_vec)-1),y_vec,'k','linewidth',1.5); 
   start_time = start_time + fx_dt;

       end

        % updating I0, which will then be treated with another radiation
        % dose
        I0 = [x_vec(fx_dt)*max(x) y_vec(fx_dt)*max(y)]./(Npoints-1);
   
       
   end
% aesthetics
   
   subB(1,1).select(); % separatrix
          set(gca,'tickdir','out','linewidth',1.5,'xtick',[0 15 30],'ytick',[0 15 30])
           g = plot((curve.x)./max(x)*(Npoints-1),(curve.y)./max(y)*(Npoints-1),'g--','linewidth',1);
           set(g,'color',[0 0.5 0])
           axis([0 33 0 31]); box off; axis square
           title([num2str(dose) ' J/cm^2 x ' num2str(num_fx) ])

   subC(2,1).select();
          set(gca,'tickdir','out','linewidth',1.5)
          ylabel('E'); %axis square
          
          if j == 1
             set(gca,'xlim',[0 20],'xtick',[0 10 20],'ylim',[0 5],'ytick',[0 2.5 5])
          elseif j == 2
              set(gca,'xlim',[0 60],'xtick',[0 30 60],'ylim',[0 30],'ytick',[0 15 30])
          elseif j == 3
              set(gca,'xlim',[0 20],'xtick',[0 10 20],'ylim',[0 10],'ytick',[0 5 10])              
          elseif j == 4
              set(gca,'xlim',[0 60],'xtick',[0 30 60],'ylim',[0 33],'ytick',[0 15 30])
          end
% 
   subC(1,1).select();
          set(gca,'tickdir','out','linewidth',1.5)
          ylabel('T'); %axis square
          
          if j == 1
             set(gca,'xlim',[0 20],'xtick',[0 10 20],'ylim',[0 30],'ytick',[0 15 30])
          elseif j == 2
             set(gca,'xlim',[0 60],'xtick',[0 30 60],'ylim',[0 30],'ytick',[0 15 30])              
          elseif j == 3
              set(gca,'xlim',[0 20],'xtick',[0 10 20],'ylim',[0 30],'ytick',[0 15 30])
          elseif j == 4
              set(gca,'xlim',[0 60],'xtick',[0 30 60],'ylim',[0 30],'ytick',[0 15 30])
          end
%    
end



B.fontsize = 14;

   B.de.margin = 15 ;
    B.export('Fig4_FxPIT.tiff','-w200','-h150', '-rp')
    