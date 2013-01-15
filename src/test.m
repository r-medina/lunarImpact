function test
% Function to plot orbit of a satellite
% launched from position (R,0) with velocity (0,V)

GM=1;
R=1;
V=1;
time=1000;
w0 = [R,.1,0,0,V,.1]; % Initial conditions

options = odeset('RelTol',0.00001);
[t_values,w_values] = ode45(@odefunc,[0,time],w0,options);

plot3(w_values(:,1),w_values(:,2),w_values(:,3))
function dwdt = odefunc(t,w)
x=w(1); y=w(2); z=w(3);
vx=w(4); vy=w(5); vz=w(6);
r = sqrt(x^2+y^2+z^2);
dwdt = [vx;vy;vz;-GM*x/r^3;-GM*y/r^3;-GM*z/r^3];
end
end