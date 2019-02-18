k = 0.8;
m = 0.078;
g = 9.8;
r = 0.78;
x_initial = 0.5;
z_initial = 0.5;
alpha_initial = 0;
W_initial = 0;
gamma = (k)/(2*m*r);
z0 = 3;
x0 = 3;
sample_time = 0.1;
time_end = 100;
t = 0:0.1:time_end;

scale = 5;
hold all

axis equal
axis([-scale scale -scale scale]);


%sim('pos_control.slx');
sim('almost_round.slx');
plot(4*sin((1/10)*t), 4*cos(1/10*t));
plot(x_out.data, y_out.data, 'r--')
dron = plot(x_out.data(1), y_out.data(1), 'k');
centre = plot(x_out.data(1), y_out.data(1), 'ko', 'MarkerSize', 4, 'MarkerFaceColor', [0 0 0]);
for i=1:numel(x_out.data)
   x = zeros(1, 5);
   x(3) = x_out.data(i); %very bad code
   x(2) = x(3) - r*cos(alpha_out.data(i));
   x(4) = x(3) + r*cos(alpha_out.data(i));
   x(1) = x(2) + cos(alpha_out.data(i)+pi/2)/5;
   x(5) = x(4) + cos(alpha_out.data(i)+pi/2)/5;
   y(3) = y_out.data(i);
   y(2) = y(3) - r*sin(alpha_out.data(i));
   y(4) = y(3) + r*sin(alpha_out.data(i));
   y(1) = y(2) + sin(alpha_out.data(i)+pi/2)/5;
   y(5) = y(4) + sin(alpha_out.data(i)+pi/2)/5;
   set(dron, 'XData', x, 'YData', y);
   set(centre, 'XData', x_out.data(i), 'YData', y_out.data(i));
   pause(0.04);
end    
