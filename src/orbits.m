%{
-------------------------------------------------------------------
Author: ramedina
Date: 3.14.11
function: orbits
-------------------------------------------------------------------
%}

function orbits

% Gravitational parameter
mu = 3.986012e5;
% A week in seconds. Gets multiplied by various scalors where
% appropriate
time = 604800;
% Radius of earth
re = 6378.145;
% Apogee of satellite
apo_sat = 35950;
% Difference in velocity when satellite is going to second impact
% point from perigee
%dv1 = 0.701;
dv1 = 0.703;
% Differnce in velocity when satellite is going to first impact
% point from apogee
dv2 = 2.544;
% Radius of moon
moon_radius = 1737.5;
% Ration of moon to earth mass
kappa = 0.01229;

% Sets moon's initial conditions and stores its perigee and apogee
[rp_moon,vp_moon,peri_moon,apo_moon] = ...
    moondata('../data/apoperi.csv',mu,moon_radius);
w_moon = [rp_moon,vp_moon];
% Sets the satellite's initial conditions and stores the normal
% vector to its orbital plane
[rp_sat,vp_sat,n] = satellite(7,250,apo_sat,178,180,re,mu);
w_sat = [rp_sat,vp_sat];
% Sets the initial conditions for when the satellite approaches the
% second impact point from perigee
[rp_sattomoon,vp_sattomoon] = sattomoon(dv1,rp_sat,vp_sat);
w_sattomoonperi = [rp_sattomoon,vp_sattomoon];

close all;

% Graphs the Earth
fig1 = figure;
hold on;
[e1,e2,e3] = sphere;
surf(e1*re,e2*re,e3*re);
axis equal;
grid on;

%{
-------------------------------------------------------------------
WITHOUT CORRECTION ON THE SATELLITE VELOCITY
-------------------------------------------------------------------
%}


% Graphs moon orbit and figures out the potential impact points
options = odeset('RelTol',0.00000001,'Event',@detect_impact_point1);
[t_moon,w_moonvals,t_impmoon,w_impmoon] = ...
    ode45(@eom,[0,time*10],w_moon,options);
plot3(w_moonvals(:,1),w_moonvals(:,2),w_moonvals(:,3));

% Graphs the two potential impact points
plot3(w_impmoon(1,1),w_impmoon(1,2),w_impmoon(1,3),'ro', ...
      'MarkerSize',7,'MarkerFaceColor','r')
plot3(w_impmoon(2,1),w_impmoon(2,2),w_impmoon(2,3),'ko', ...
      'MarkerSize',7,'MarkerFaceColor','k')

% Graphs the satellite and figures out its coordinates at apogee
options = odeset('RelTol',0.000000001,'Event',@satapogee1);
[t_sat,w_satvals,t_satopogee,w_satapogee] = ...
    ode45(@eom,[0,time*10],w_sat,options);
plot3(w_satvals(:,1),w_satvals(:,2),w_satvals(:,3));

% Graphs the satellite's trajectory towards the second impact point
% from perigee and prints the travel time and distance
[t_sattomoon,w_sattomoonvals,distancefromperigee,timefromperigee] = ...
    minimum(time*4,w_sattomoonperi,w_impmoon(2,1:3));
plot3(...
    w_sattomoonvals(:,1),w_sattomoonvals(:,2),w_sattomoonvals(:,3));
distancefromperigee
timefromperigee

% Sets up inital conditions for satellite leaving from apogee
rp_sattomoon = w_satapogee(1,1:3);
vp_sattomoon = w_satapogee(1,4:6);
[rp_sattomoon,vp_sattomoon] = ...
    sattomoon(dv2,rp_sattomoon,vp_sattomoon);
w_sattomoonapo = [rp_sattomoon,vp_sattomoon];

% Graphs the satellite's trajectory towards the first impact point
% from apogee and prints the travel time and distance
[t_sattomoon,w_sattomoonvals,distancefromapogee,timefromapogee] = ...
    minimum(time*5,w_sattomoonapo,w_impmoon(1,1:3));
plot3(...
    w_sattomoonvals(:,1),w_sattomoonvals(:,2),w_sattomoonvals(:,3));
distancefromapogee
timefromapogee

%{
-------------------------------------------------------------------
ACCOUNTING FOR MOON GRAVITY ON SATELLITE
-------------------------------------------------------------------
%}

fig2 = figure;
hold on;
[e1,e2,e3] = sphere;
surf(e1*re,e2*re,e3*re);
axis equal;
grid on;

options = odeset('RelTol',0.00000001,...
		 'Events',@detect_impact_point);
[t_vals,w_vals,t_impmoon,w_impmoon] = ...
    ode45(@eom2,[0,time*8],[w_sat,w_moon],options);
plot3(w_vals(:,1),w_vals(:,2),w_vals(:,3));
plot3(w_vals(:,7),w_vals(:,8),w_vals(:,9));
timetoimpact = t_impmoon(2);

% Graphs the potential impact point
plot3(w_impmoon(2,7),w_impmoon(2,8),w_impmoon(2,9),'ko', ...
      'MarkerSize',7,'MarkerFaceColor','k')

timetoimpact;
tmooni = timetoimpact-timefromperigee;

for i = 1:length(t_vals)
  if t_vals(i) >= tmooni
    ind = i;
    break;
  end
end

options = odeset('RelTol',0.00000001,...
		 'Events',@stopper);
[t_vals2,w_vals2,t_stop,w_stop] = ...
    ode45(@eom2,[0,time*8],[w_sattomoonperi,w_vals(ind,7:12)],options);
plot3(w_vals2(:,1),w_vals2(:,2),w_vals2(:,3));

timetoimpact-t_stop

fig3 = figure;

for i = 1:length(t_vals2)
  clf;
  hold on;
  axis equal;
  grid on; 
  surf(e1*re,e2*re,e3*re);
  plot3(w_vals(:,1),w_vals(:,2),w_vals(:,3));
  plot3(w_vals2(:,1),w_vals2(:,2),w_vals2(:,3),'k');
  plot3(w_vals(:,7),w_vals(:,8),w_vals(:,9));
  plot3(...
      w_vals2(i,1),w_vals2(i,2),w_vals2(i,3),'ro','MarkerSize',7,'MarkerFaceColor','r');
  plot3(...
      w_vals2(i,7),w_vals2(i,8),w_vals2(i,9),'ko','MarkerSize',5,'MarkerFaceColor','c');
%  moonloc = [w_vals2(i,7),w_vals2(i,8),w_vals2(i,9)];
  pause(0.0001);
end

[m1,m2,m3] = sphere;
surf(m1*moon_radius+w_stop(7),m2*moon_radius+w_stop(8),m3*moon_radius+w_stop(9));

velocity_at_impact = ...
    norm([w_stop(10)-w_stop(4),w_stop(11)-w_stop(5),w_stop(12)-w_stop(6)])


%{
-------------------------------------------------------------------
FOR PREMILINARY REPORT                            
-------------------------------------------------------------------
%}

%{
% Velocity of satellite while in orbit
v_sat = w_satvals(:,4:6);
% Speed of satellite while in orbit
speed_sat = (w_satvals(:,4).^2+w_satvals(:,5).^2+w_satvals(:,6).^2).^.5;
% Position vector of satellite while in orbit
r_sat = w_satvals(:,1:3);
% Distance to center of eartch of satellite while in orbit
d_sat = (w_satvals(:,1).^2+w_satvals(:,2).^2+w_satvals(:,3).^2).^.5;
% Total energy of satellite while in orbit
e_sat = .5 * speed_sat.^2 - mu./d_sat;

% Graphs the total energy per unit mass of the satellite
%figure;
%plot(t_sat,e_sat);
%title('Graph of total energy of satelite per unit mass');
%axis([0 time*2 min(e_sat) max(e_sat)]);

% Finds the angular momentum per unit mass of the satellite while
% in orbit by taking the cross product of its position vector and
% velocity vector at every time index
%for i = 1:length(t_sat);
%  L_sat(i) = norm(cross(r_sat(i,:),v_sat(i,:)));
%end

% Graphs angular momentum per unit mass of the satellite while in orbit
%figure;
%plot(t_sat,L_sat);
%title('Graph of angular momentum of satellite per unit mass');
%axis([0 time*2 min(L_sat) max(L_sat)]);

% Prints the perigee of the moon
peri_moon
% Prints the apogee of the moon
apo_moon

% Displays the number of days that it takes for the moon to go
% around the earth
moontime = (t_impmoon(3) - t_impmoon(1))/60/60/24;
disp(sprintf('Moon orbits in %d days',moontime));
%}

%{
-------------------------------------------------------------------
FUNCTIONS
-------------------------------------------------------------------
%}

% Equation of motion as a functoin of time and a position vector
% and velocity vector
function dwdt = eom(t,w)
  x=w(1); y=w(2); z=w(3);
  vx=w(4); vy=w(5); vz=w(6);
  r = sqrt(x^2+y^2+z^2);
  dwdt = [vx;vy;vz;-mu*x/r^3;-mu*y/r^3;-mu*z/r^3];
end

function [event_val,stopthecalc,direction] = ...
    detect_impact_point1(t,w)
  % Position vector of moon (this assumes w(1)=x,w(2)=y,w(3)=z)
  r = w(1:3);
  % Detect when r.n=0
  event_val = dot(n,r);
  stopthecalc = 0;
  direction = 0;
end

function [event_val,stopthecalc,direction] = satapogee1(t,w);
  event_val = sqrt(w(1)^2+w(2)^2+w(3)^2)-(apo_sat+re);
  stopthecalc = 0;
  direction = 0;
end

function [tvals,wvals,test_d,time_to_reach] = minimum(t,w,r_imp)
  options = odeset('RelTol',0.00000001,'Event',@min_dist);
  [tvals,wvals,tevent,wevent] = ...
      ode45(@eom,[0,t],w,options);
 
  [test_d,time_to_reach] = ...
      mindisttomoon(tevent,wevent,r_imp);

  function [event_val,stopthecalc,direction] = ...
      min_dist(t,w)
    event_val = dot((r_imp - transpose(w(1:3))),w(4:6));
    stopthecalc = 1;
    direction = 0;
  end

  function [d,time_to_reach_min] = ...
      mindisttomoon(t_event,w_event,r_impact)
    rmin = w_event(1,1:3)-r_impact;
    d = sqrt(dot(rmin,rmin)); % min dist to moon.
    time_to_reach_min = tevent(1); % Time to reach the min dist.
  end

end

% The new and improved equation of motion that solves for the
% effect of the moon's gravity on the satellite
function dWdt = eom2(t,w)
  xs=w(1); ys=w(2); zs=w(3);
  vxs=w(4); vys=w(5); vzs=w(6);
  xm=w(7); ym=w(8); zm=w(9);
  vxm=w(10); vym=w(11); vzm=w(12);
  rs = sqrt(xs^2+ys^2+zs^2);
  rm = sqrt(xm^2+ym^2+zm^2);
  rmosat = sqrt((xm-xs)^2+(ym-ys)^2+(zm-zs)^2);
  % kappa is the ratio between the moon's mass and the earth's
  dWdt = [vxs;vys;vzs;...
	  -mu*xs/rs^3+mu*kappa*(xm-xs)/rmosat^3;...
	  -mu*ys/rs^3+mu*kappa*(ym-ys)/rmosat^3;...
	  -mu*zs/rs^3+mu*kappa*(zm-zs)/rmosat^3;...
	  vxm;vym;vzm;...
	  -mu*xm/rm^3;-mu*ym/rm^3;-mu*zm/rm^3];
end

function [event_val,stopthecalc,direction] = ...
    detect_impact_point(t,w)
  % Position vector of moon (this assumes w(1)=x,w(2)=y,w(3)=z)
  r = w(7:9);
  % Detect when r.n=0
  event_val = dot(n,r);
  stopthecalc = 0;
  direction = 0;
end

function [event_val,stopthecalc,direction] = satapogee(t,w);
  event_val = sqrt(w(1)^2+w(2)^2+w(3)^2)-(apo_sat+re);
  stopthecalc = 0;
  direction = 0;
end

function [event_val,stopthecalc,direction] = stopper(t,w)
  r = sqrt((w(7)-w(1))^2+(w(8)-w(2))^2+(w(9)-w(3))^2);
  event_val = r - moon_radius;
  stopthecalc = 1;
  direction = 0;
end

end
