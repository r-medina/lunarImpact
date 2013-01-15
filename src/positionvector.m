% now obsolete
function [t_vals,w_vals,tevent,wevent,indexevent] = ...
    positionvector(t,rp,vp,grav_par,normal,which)

w0 = [rp,vp];

if which == 1
  options = odeset('RelTol',0.0000001,'Event',@detect_impact_point);
  [t_vals,w_vals,tevent,wevent,indexevent] = ...
      ode45(@eom,[0,t],w0,options,normal);
end

% $$$ if which == 2
% $$$   options = odeset('RelTol',0.00000001,'Event',@min_dist);
% $$$   [t_Vail's,w_vals,tevent,wevent,indexevent] = ...
% $$$       ode45(@eom,[0,t],w0,options);
% $$$ end

if which == 0
  options = odeset('RelTol',0.0000001);
  [t_vals,w_vals] = ode45(@eom,[0,t],w0,options);
end

function [event_val,stopthecalc,direction] = ...
    detect_impact_point(t,w,n)
  % Position vector of moon (this assumes w(1)=x,w(2)=y,w(3)=z)
  r = w(1:3);
  % Detect when r.n=0
  event_val = dot(n,r);
  stopthecalc = 0;
  direction = 0;
end

function dwdt = eom(t,w,n)
  x=w(1); y=w(2); z=w(3);
  vx=w(4); vy=w(5); vz=w(6);
  r = sqrt(x^2+y^2+z^2);
  dwdt = [vx;vy;vz;-grav_par*x/r^3;-grav_par*y/r^3;-grav_par*z/r^3];
end

end