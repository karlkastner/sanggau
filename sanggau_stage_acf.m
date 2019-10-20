% 2015-07-24 18:06:43.674289728 +0200
% Karl Kastner, Berlin

meta = sanggau_metadata();
% preload water level
load([ROOTFOLDER,'/src/water-level/mat/water-level.mat']);
level.val  = K(nid.sanggau_merged).slot.depth;
level.time = K(nid.sanggau_merged).slot.time;


% preload vadcp calibrarion discharge

% daily average
val = level.val;
val0 = val;
n   = length(val);
dt  = level.time(2)-level.time(1);
m   = round(n*dt);
l   = round(1/dt); 
val = reshape(val(1:m*l),l,m);
t   = reshape(level.time(1:m*l),l,m);
val = cvec(nanmean(val));
t   = cvec(mean(t));
[a] = acf_man(val);
%[a dd] = acf_man(val);

clf
subplot(2,2,1)
plot(t,val)
datetick('x','mmm');
subplot(2,2,2)
plot(a)

% what about instationarity?

k = find(a<exp(-1),1,'first');
printf('Correlation length in days %d\n',k);
rho = exp( -1/k )
k_ = (1*rho)/(1-rho)

dof = sum(isfinite(val))/k;
printf('Degrees of freedom %f\n',dof);

printf('Subtracting seasonal cycle');
fdx = isfinite(val);
A   = [ones(size(t)) sin(2*pi*t/365.25) cos(2*pi*t/365.25)];
c   = A(fdx,:)\val(fdx)
val_ = val-A*c;
a   = acf_man(val_);
subplot(2,2,1);
hold on;
plot(t,A*c);
subplot(2,2,2);
hold on
plot(a)

k = find(a<exp(-1),1,'first');
printf('Correlation length w/o seasonal variation in days %d\n',k);
a = a-exp(-(0:length(a)-1)'./k);
plot(a)
legend('Total','w/o seasonal variation','w/n seasonal nor first order correlation');

subplot(2,2,3)
plot(t,[val-nanmean(val),val_])
nanrms(val0-nanmean(val0))
nanrms(val-nanmean(val))
nanrms(val_)
legend('Water level','Water level w/o seasonal variation')
datetick('x','mmm');



