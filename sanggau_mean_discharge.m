% 2016-01-22 00:28:16.624541754 +0100
load ../water-level/mat/water-level.mat
id = nid.sanggau_merged;
t=K(id).slot.time;
dt = t(2)-t(1);
d=K(id).slot.depth;
% TODO quick fix (slot erronously extrapolated)
fdx = t>735517 & t<736081;
d(~fdx) = NaN;

dt_max = 1;
d = fixnan(t,d,1);
%fdx = isnan(d);
%d(fdx) = interp1(find(~fdx),d(~fdx),find(fdx));
%d=d-max(d)+10;

figure(1);
clf()
subplot(2,1,1)
plot(t,d)
[nanmedian(d) nanmean(d)]

%fdx = isfinite(d);
a = acf_man(d,round(365/dt));
subplot(2,1,2)
plot(dt*(1:length(a)),a)

[rho T] = ar1delay(a);
T = dt*T
vline(T)

%de = rho
%dt = -1/log(de)*median(diff(K(23).slot.time))
% Note, this is twice as much as estimated by Meyer-Ragu
100*1e9 / (365*86500)

