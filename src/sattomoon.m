function [rp,vp] = sattomoon(delta_v,r,v)

rp = r;
vp = v + delta_v*v/sqrt(dot(v,v));

end