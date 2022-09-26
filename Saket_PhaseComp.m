figure(100);clf
B = panel();
    
    B.pack([0.5 0.25 0.25],1)

B(1,1).select()
    
sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;   
    [curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);

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


r = area(curve.x./max(x)*(Npoints-1),curve.y./max(y)*(Npoints-1)+100)

set(r,'FaceColor',[1 0 0],'EdgeColor',[1 0 0],'FaceAlpha',0.5,'EdgeAlpha',0.1,'linewidth',0.01)
hold on
    
   g1 = area(curve.x./max(x)*(Npoints-1),curve.y./max(y)*(Npoints-1))

set(g1,'FaceColor',[1 1 1],'EdgeColor',[1 1 1])

   g2 = area(curve.x./max(x)*(Npoints-1),curve.y./max(y)*(Npoints-1))

set(g2,'FaceColor',[0 1 0],'EdgeColor',[0 1 0],'FaceAlpha',0.5,'EdgeAlpha',0.5,'linewidth',0.01)
a = [max(curve.x./max(x)*(Npoints-1))  30 30 max(curve.x./max(x)*(Npoints-1))]
b = [0                                  0 30  30]
%patch([xu(:)' fliplr(xu(:)')], [ymin fliplr(ymax)], [0.1 0.9 0.1], 'FaceAlpha',0.3) % Plot Filled Background
q = patch(a,b,'b')
set(q,'FaceColor',[0 1 0],'EdgeColor',[0 1 0],'FaceAlpha',0.5,'linewidth',0.01)
hold on

g = plot(curve.x./max(x)*(Npoints-1), curve.y./max(y)*(Npoints-1),'--')

   set(g,'linewidth',1.25,'color',[0 0 0.75])

xlabel('Immune Effector Cells (10^6 cells)'); ylabel('Tumor Cells (10^6 cells)')
axis([0 22 0 30]); box off; axis square
set(gca,'tickdir','out','linewidth',1,'fontsize',8,'xtick',[],'ytick',[],...
    'xticklabel',[],'yticklabel',[])


hold on
[X1, Y1] = meshgrid(0:Npoints-1,0:Npoints-1);

 q = quiver(X1,Y1,U,V); %plotting vector field
 q.Color = [0 0 0]; 
 q.AutoScaleFactor = 0.5;
 
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

B(2,1).select()
for i = 1:size(initCond,1)
 h=plot(sols{i}.y(2,:)/max(y)*(Npoints-1))
 set(h,'linewidth',1.5)
 hold on
end
set(gca,'linewidth',1.5,'tickdir','out','fontsize',14); xlabel('time');ylabel('T')
axis([1 45 0 30])
B(3,1).select()
for i = 1:size(initCond,1)
 h=plot(sols{i}.y(1,:)/max(x)*(Npoints-1))
 set(h,'linewidth',1.5)
 hold on
end
set(gca,'linewidth',1.5,'tickdir','out','fontsize',14); xlabel('time');ylabel('E')
axis([1 45 0 25])