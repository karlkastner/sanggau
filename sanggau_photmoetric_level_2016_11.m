% Fri 10 Feb 10:40:50 CET 2017

% TODO measurement 7 is too far off
% TODO compare to the bed level offset

meta = sanggau_metadata();                                                      

dt_max = 1;
load ../water-level/mat/water-level.mat
load ../discharge/mat/sanggau-stage-discharge.mat

clear t z

% calibration 1 (installation)
% image IMG_0213
% camera time + 6
t(1) = datenum('2013/12/10 10:32:38');
m_per_pixel = 0.25/863;
l(1) = 10-1424*m_per_pixel;

% read out 2014 01 no picture of gage and no calibration (no boat)

% calibration 2 - before redeployment
% IMG_1155.JPG
% camera time + 0
t(2) = datenum('2014/02/18 07:58:44');
m_per_pixel = (6-4)/2288;
l(2) = 4-1940*m_per_pixel;

% calibration 2 - after redeployment
% IMG_1213 
t(3) = datenum('2014/02/20 08:11:11');
m_per_pixel = (5-3)/1939;
l(3) = 3 - 835*m_per_pixel;

% 
% IMG_1776
t(4) = datenum('2014/04/17 22:20:52');
m_per_pixel = (10-5)/3364;
l(4) = 5 - 76*m_per_pixel;

% read out - no calibration (no large enough boat)
% IMG_4838
% this was slanted (difficult)
t(5) = datenum('2014/07/26 01:28:32');
m_per_pixel = (6-3)/1047;
l(5) = 3-792*m_per_pixel;

% removal of instrument
% 0134
t(6) = datenum('2015/04/26 07:56:58');
l(6) = 7; %m (exact match of marker)

% read out by student
% DSC01733
t(7)        = datenum('2014/03/22 15:14:17');
m_per_pixel = (8.75-5.75)/1938;
l(7)        = 61*m_per_pixel;

% 0666
t(8) = datenum('2014/01/11 09:11:27');
m_per_pixel = 0.25/275;
% elevation of yellow dot in 0216
zp = 10-(462-154)*m_per_pixel; 
% distance between the two yellow dots in 12/2016
% slanted distance
d1 = 315*m_per_pixel;
% distance along the vertical only
d1 = (295/315)*d1;

% same elevation of yellow dot in 0666
d2 = d1;
% the 464 pixels are also slanted, but at a different angle
m_per_pixel = d2/464;
l(8) = zp - 714*m_per_pixel;
l = cvec(l);

% DSCF0158
% local time from gps
t(9) = datenum('2016/11/05 11:40:21')
m_per_pixel = (10-8)/1044
l(9) = 8-529*m_per_pixel;

id      = nid.sanggau_merged;
fdx     = isfinite(K(id).depth);
l_gauge = interp1_limited(K(id).time(fdx),K(id).depth(fdx),t,dt_max,'linear');

l_gauge = l_gauge;
%- meta.z_offset;

dl     = l-l_gauge;
offset = nanmedian(dl);
l      = l-offset;

l       = l - meta.z_offset;
l_gauge = l_gauge - meta.z_offset;
clear Q
% predict discharge for reference campaign
Q(:,1) = powerrc.predict(NaN(size(l_gauge)),l_gauge);
Q(:,2) = powerrc.predict(NaN(size(l)),l);

% TODO offset from bed level profile is 2.44m only (but gps antenna antenna altitude correction maybe too small)
% there might also be an issue with the redeployment
% note that the bed level might actually be lower, consistend with the tide at low flow
fprintf('offset (l0 2016 11): %g\n', offset);

figure(1);
clf
subplot(2,2,1)
plot(dl,'.');
ylabel('l_{optical} - l_{gauge} (m)')
hline(offset);

subplot(2,2,2)
plot(sort(dl),'.')

subplot(2,2,3)
plot(l_gauge,'o');
hold on
plot(l,'*');
ylabel('Surface level (m)')

subplot(2,2,4)
plot(Q(:,1)/1e3,'o')
hold on
plot(Q(:,2)/1e3,'*')
ylabel('Discharge (10^3 m^3/s)');

