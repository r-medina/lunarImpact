function [rp,vp,perigee,apogee] = moondata(file,grav_par,moonrad)

% Converts degrees to radians
radians = @ (degrees) degrees*pi/180;

% Reads a csv and makes it an array
moondata = csvread(file);

% Stores the perigee and apogee and where they happen
[perigee,pind] = min(moondata(:,8));
[apogee,aind] = max(moondata(:,8));

perigee = perigee + moonrad;
apogee = apogee + moonrad;

% Stores the ra and dec at perigee and apogee
[ra_p,dec_p] = radec(pind);
[ra_a,dec_a] = radec(aind);

% Makes the position vectors at perigee and apogee.
rp = rvec(ra_p,dec_p,perigee);
ra = rvec(ra_a,dec_a,apogee);

% Crosses the position vectors at perigee and apogee to find a
% vector normal to the moon's plane
m = cross(rp,ra)./(norm(cross(rp,ra)));

% Finds velocity when moon is at perigee
vp = sqrt(2*grav_par*apogee/(perigee*(apogee+perigee))) * ...
     (cross(m,rp)/perigee);

% Function that outputs a position vector given the distance, ra,
% and dec (in radians)
function r = rvec(ra,dec,dist)
  r = dist * ...
      [cos(ra)*cos(dec), sin(ra)*cos(dec), sin(dec)];
end


% Function that outputs the ra and dec in radians for a given index
% of moondata.csv
function [ra,dec] = radec(index)
  % Stores the ra and dec in the file as some variables
  ra_deg = moondata(index,1:3);
  dec_deg = moondata(index,4:7)*moondata(index,4);  

  % Converts ra_deg first to decimal degrees
  ra = ra_deg(1)*(360/24) + ...
       ra_deg(2)*360/(24*60) + ...
       ra_deg(3)*60/(24*60*60);
  % Then to radians and stores it as the output of the function
  ra = radians(ra);
  
  % Converts dec values in the file to decimal degrees then radians
  % and stores it as the second output of the function
  dec = dec_deg(2) + dec_deg(3)/60 + dec_deg(4)/(60*60);
  dec = radians(dec);
end

end