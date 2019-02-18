alpha = 1;  %constants input
beta = 0.2;
mu1 = 0.1;
mu2 = 0.3;
sigma = 4;
xlimit = 150;
ylimit = 300; %end of constants


%search of stationary point
ys = alpha;
D = alpha*alpha*beta*beta + 4*alpha*mu2*(sigma - alpha*mu1);
xs = (-alpha*beta - sqrt(D))/(-2*alpha*mu2); 


xlimit = xs+1;
ylimit = ys+1;

%making Stapanova system enable to used by ode45 solver
Stepa = @(t, y)Stepanova_sys(t, y, alpha, beta, mu1, mu2, sigma); 
%end of function

%making vector representation of phase portret
[x1,y1] = meshgrid(linspace(0,15,40), linspace(0,15,80));

u = zeros(size(x1));
v = zeros(size(x1));

for i = 1:numel(x1)
    uv = Stepa(0, [x1(i); y1(i)]);
    u(i) = uv(1);
    v(i) = uv(2);
end

quiver(x1,y1,u,v);
%here vector representation ends


t1 = [0 1]; %%time is specified by system


%drawing of phase portrait
hold on
plot(xs, ys, 'k*');%stationary point
Phase_portrait(Stepa, t1, xlimit, ylimit, 14, 8);
hold off
%end of drawing

axis([xs-1 xs+1 ys-1 ys+1]);

xlabel('x')
ylabel('y')
