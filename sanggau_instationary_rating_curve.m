% Fri  3 Mar 11:43:45 CET 2017
% estimate of instationary termes on the rating curve
u = 1;
c = 3/2*u;
% this is similar in the surface level measurement in 12/2016 and srtm data
% 
S = 2.5e-5;
% this is the standard deviation of bed level in 1 day
dh_dt=0.2/86400;

d = 1-sqrt(1+1/(c*S)*dh_dt)
d = 1-sqrt(1-1/(c*S)*dh_dt)
d = 0.5*1/(c*S)*dh_dt

