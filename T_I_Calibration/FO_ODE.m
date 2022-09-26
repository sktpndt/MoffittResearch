%% stable vs unstable fixed points
f = @(t,x) 9-x^2
[ts,xs] = ode45(f,[-1,5],1);
plot(ts,xs,'o-')
% stable point at x = 3 as t --> inf