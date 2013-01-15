function [rp,vp,normalvector] = satellite(inc,peri_alt,apo_alt,arg,long,rb,grav_par)

radians = @ (degrees) degrees*pi/180;

inc = radians(inc);
arg = radians(arg);
long = radians(long);

peri_alt = peri_alt + rb;
apo_alt = apo_alt + rb;

rp = (peri_alt) * ...
     [ (cos(arg)*cos(long)-sin(arg)*sin(long)*cos(inc)), ...
       (cos(arg)*sin(long)+sin(arg)*cos(long)*cos(inc)), ...
       (sin(arg)*sin(inc)) ];

vp = sqrt(2*grav_par*apo_alt/(peri_alt*(apo_alt+peri_alt))) * ...
     [ -(sin(arg)*cos(long)+cos(arg)*sin(long)*cos(inc)), ...
     (cos(arg)*cos(long)*cos(inc)-sin(arg)*sin(long)), ...
     (cos(arg)*sin(inc)) ];

normalvector = [sin(long)*sin(inc), -cos(long)*sin(inc), cos(inc)];

end