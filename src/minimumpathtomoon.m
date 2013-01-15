% now obsolete
function [tvals,wvals,test_d,time_to_reach] = ...
    minimumpathtomoon(t,w,r_imp)

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