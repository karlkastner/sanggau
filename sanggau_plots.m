% 2014-07-14 13:06:24.237653350 +0200
% Karl Kastner, Berlin

scale_sintang = 2;
pflag = 1;
d_offset = 2.25;
or = 2;

C = { ...
'2013-12-09-sanggau.mat', ...
'2013-02-20-sanggau.mat', ...
'2013-04-18-sanggau.mat', ...
'2013-06-18-sanggau.mat' };

for idx=1:length(C)
	C{idx}
	load(['mat/' C{idx}]);
	d(idx,1) = discharge;
	u0v(idx,1) = sqrt(sum(discharge.ubar(1:2).^2));
	q0v(idx,1) = sqrt(sum(discharge.discharge(1:2).^2));
	h0v(idx,1) = discharge.area/discharge.width;
	val(idx,1) = u0v(idx);
	val(idx,2) = q0v(idx)/10000;
	val(idx,3) = h0v(idx)/10;
	t0(idx,1) = mean(discharge.ens.time);
end % for idx
w = discharge.width;

'regressing u = exp(c0 + c1*h)'
c = [h0v.^0 h0v] \ log(u0v)
u_est = exp(c(1) + c(2)*h0v);
[u0v u_est (u_est-u0v)./(u0v)]
'regressing q = exp(c0 + c1*h)'
c = [h0v.^0 h0v] \ log(q0v)
q_est = exp(c(1) + c(2)*h0v);
[q0v q_est (q_est-q0v)./(q0v)]

'regressing [u;q-hw] '
one = ones(size(h0v));
c = [[one; one] [h0v*w; h0v*w]] \ [log(u0v); log(q0v)-log(h0v*w)];

u_reg = exp(c(1) + c(2)*h0v*w);
q_reg = (h0v*w).*exp(c(1) + (c(2)+1));
c(1)
c(2)
[ u0v u_reg]
mean( (u_reg-u0v).^2)./var(u0v,1)
[  q0v q_reg]
mean((q_reg-q0v).^2)./var(q,1)

%figure(1);
%clf();
%bar(u0v);
%set(gca,'xticklabel',datestr(t0,'dd/mm/yy'));
%print('img/sanggau-calibration-velocity.eps');
%system('LD_LIBRARY_PATH= epstopdf img/sanggau-calibration-velocity.eps');
%ylabel('m/s');
%
%figure(2);
%clf();
%bar(q0v);
%set(gca,'xticklabel',datestr(t0,'dd/mm/yy'));
%print('img/sanggau-calibration-discharge.eps');
%system('LD_LIBRARY_PATH= epstopdf img/sanggau-calibration-discharge.eps');
%ylabel('m^3/s');
%
%figure(3);
%clf();
%bar(h);
%set(gca,'xticklabel',datestr(t0,'dd/mm/yy'));
%print('img/sanggau-calibration-mean-depth.eps');
%system('LD_LIBRARY_PATH= epstopdf img/sanggau-calibration-mean-depth.eps');
%ylabel('m');
%
figure(1);
clf();
subplot(1.5,1,1);
bar(val);
set(gca,'xticklabel',datestr(t0,'dd/mm/yy'));
grid on
legend('location','eastoutside','mean velocity [m/s]','discharge [1e4 m^3/s]', 'mean depth [10 m]');
print('img/sanggau-calibration.eps');
system('LD_LIBRARY_PATH= epstopdf img/sanggau-calibration.eps');

figure(2);
clf();
h_ = linspace(0,max(h0v),100);
u_ = exp(c(1) + c(2)*h_*w);
plot(h0v,u0v,'ko','markerfacecolor',[0 0 0]);
hold on
plot(h_,u_,'k-');
xlabel('average depth [m]');
ylabel('average velocity [m/s]');
t_str = datestr(t0,'  dd/mm/yy');
text(double(h0v),double(u0v),t_str);
grid on
print('img/sanggau-h-u.eps');
system('LD_LIBRARY_PATH= epstopdf img/sanggau-h-u.eps');

% Mon Jul 14 00:18:15 WIB 2014
% Karl Kastner, Berlin

load('../hadcp/mat/hadcp-sanggau.mat');
adcp = hadcp;
% dummy bottom track data
adcp.btvel = zeros(size(adcp.VEL,2),4,'int16');

% convert raw values
time                 = ADCP.convert_raw_time(adcp.timeV);
%for idx=1:4
%subplot(2,2,idx)
%imagesc(adcp.VEL(:,:,idx))
%end
%pause
[vel btvel]          = HADCP.convert_raw_velocity(adcp.VEL,adcp.btvel);
[heading pitch roll] = HADCP.convert_raw_attitude_angles(adcp.heading,adcp.pitch,adcp.roll); 
[pressure depth]     = HADCP.convert_raw_pressure(adcp.pressure);

% transform to earth velocities
[vel btvel] = HADCP.to_earth(adcp.corstr,adcp.FileNumber,vel,btvel,heading,pitch,roll);
%very hot quick fix


% quick fix for firmware bug:
% invalidate all except every fourth sample
vel(2:4:end,:,:) = NaN;
vel(3:4:end,:,:) = NaN;
vel(4:4:end,:,:) = NaN;
% velocity magnitude
vmag=sqrt(sum(vel(:,:,1:2).^2,3));
fdx = find(vmag > 2.0);
vmag(fdx) = NaN;
depth = depth(:) ;                        
uh = nanmean(vmag(1:50,:));

% skip fast ping samples
time_ = time;
time  = (time(1):1/48:time(end))';
fdx = find(isfinite(depth));
depth = interp1(time_(fdx),depth(fdx),time,'nearest');
uh    = interp1(time_,uh,time,'nearest');

%jdx=1;
%for idx=2:length(time)
%	dt = time(idx)-time(jdx);
%	if (dt > 1/49)
%		t(
%	end
%end

% correct for redeployment
% TODO, this is approximate
%fdx = find(time > t0(2));
%depth(fdx) = depth(fdx)-d_offset;

% water level
figure(3);
clf();
plot(time,depth,'.');
hold on
datetick();
h0h = interp1(time,depth,t0+1/24);
plot(t0,h0v,'ro','markerfacecolor',[1 0 0],'markersize',10);
ylabel('depth as recorded (m)'); %water column above initial instrument installation level [m]');
xlim([time(1)-1/24 time(end)+1/24]);
grid on
if (pflag)
print('img/sanggau-water-level.eps','-depsc');
system('LD_LIBRARY_PATH= epstopdf img/sanggau-water-level.eps');
end

% hadcp velocity magnitude image
figure(4);
clf
%imagesc(uh);
%hist(uh(:),linspace(0,2,100)); xlim([0 1.8])
%figure(2);
fdx=find(isfinite(vel(:,1,1)));
imagesc(time_,fdx,vmag(fdx,:))
caxis([0 2]);
colorbar()

% time series
uh_ = uh;
% hold fix
for idx=2:length(uh)
	if (~isfinite(uh(idx)))
		uh(idx) = uh(idx-1);
	end
	if (~isfinite(uh_(end-idx+1)))
		uh_(end-idx+1) = uh_(end-idx+2);
	end
end
uh = 0.5*(uh + uh_);
u0h = interp1(time,uh,t0+1/24);
nf = 48;
uh_ = [uh(1)*ones(nf,1); uh(:); uh(end)*ones(nf,1)];
uh_ = medfilt1(double(uh_),nf);
uh_ = uh_(nf+1:nf+length(uh));

% vadcp velocity magnitude time series (1d)
figure(5);
clf
plot(time,uh_,'.');
ylabel('m/s');
datetick();
xlim([time(1) time(end)]);
grid on
%if (pflag)
%	print('img/sanggau-hadcp-velocity.eps','-depsc');
%	system('LD_LIBRARY_PATH= epstopdf img/sanggau-hadcp-velocity.eps');
%end

% relation since last flood onset
fdx2 = find(time > 7.356504791666666e+05 & depth > 0.15);
% u = exp(c0 - c1*width*depth)
% TODO joint regression with two constants
uh = uh(:);
c = [ones(size(fdx2)) w*depth(fdx2)] \ log(uh(fdx2));
u_ = exp(c(1)+c(2)*w*depth(fdx2));
hold on
plot(time(fdx2),u_,'r.');
'R^2'
1 - mean((u_ - uh(fdx2)).^2)/var(uh(fdx2),1)
1 - mean((u_ - uh_(fdx2)).^2)/var(uh_(fdx2),1)


% regress first part
fdx1 = find(time <= 7.356356250000000e+05 & isfinite(depth) & (time <= 7.356202083333334e+05 | time > 7.356258125000000e+05));
uh =uh(:);                                                                  
c = [ones(size(fdx1)) w*depth(fdx1)] \ log(uh(fdx1));                            
u_ = exp(c(1)+c(2)*w*depth(fdx1));                                      
hold on                                               
plot(time(fdx1),u_,'.','color',[255 165 0 ]/255);
plot(t0,u0v,'ro','markerfacecolor',[1 0 0],'markersize',10);
'R^2'
1 - mean((u_ - uh(fdx1)).^2)/var(uh(fdx1),1)

% extract hadcp velocity at calibrations
%for idx=1:length(t0)
%	[mv mdx] = min( (time-t0(idx)).^2 );
%	uh0v(idx) = uh(mdx);
%	hh0v(idx) = depth0v(mdx);
%end
%
if (pflag)
	print('img/sanggau-hadcp-velocity.eps','-depsc');
	system('LD_LIBRARY_PATH= epstopdf img/sanggau-hadcp-velocity.eps');
end

figure(6);
clf();
subplot(1.5,1,1);
%h0_ = h0 + (h0v(1)-h0(1));
bar([h0v/10 h0h/10 u0v u0h])
set(gca,'xticklabel',datestr(t0,'dd/mm/yy'));
grid on
legend('location','eastoutside','vadcp depth [10m]','hadcp depth [10m]', 'vadcp velocity [m/s]','hadcp velocity [m/s]');

if (pflag)
	print('img/sanggau-hadcp-vadcp.eps','-depsc');
	system('LD_LIBRARY_PATH= epstopdf img/sanggau-hadcp-vadcp.eps');
end

disp('correspondence between hadcp and vadcp during the first campaign');
fdx = 1;
cu1 = u0h(fdx) \ u0v(fdx);
rr_u1 = 1-mean((u0v(fdx) - cu1*u0h(fdx)).^2)/var(u0v(fdx),1);
fprintf('u_v = %f u_h\nRR_u %f\n',cu1,rr_u1);
ch1 = ones(length(fdx),1) \ (h0v(fdx)-h0h(fdx));
rr_h1 = 1 - mean((h0v(fdx) - h0h(fdx) - ch1).^2)/var(h0v(fdx),1);
fprintf('h_v = h_h + %f\nRR_h %f\n',ch1,rr_h1);

disp('correspondence between hadcp and vadcp during the last three campaigns');
fdx = 2:length(u0v);
cu2 = u0h(fdx) \ u0v(fdx);
rr_u2 = 1-mean((u0v(fdx) - cu2*u0h(fdx)).^2)/var(u0v(fdx));
fprintf('u_v = %f u_h\nRR_u %f\n',cu2,rr_u2);
ch2 = ones(length(fdx),1) \ (h0v(fdx)-h0h(fdx));
rr_h2 = 1 - mean((h0v(fdx) - h0h(fdx) - ch2).^2)/var(h0v(fdx));
fprintf('h_v = h_h + %f\nRR_h %f\n',ch2,rr_h2);

% now that u ad h are parametrised the discharge can be estimated
hh = depth;
hh_ = NaN(size(hh));
uh_ = NaN(size(hh));
hh_(fdx1) = hh(fdx1)+ch1;
hh_(fdx2) = hh(fdx2)+ch2;
uh_(fdx1) = cu1*uh(fdx1);
uh_(fdx2) = cu2*uh(fdx2);

qh = w*uh_.*hh_;

figure(7);
clf();
subplot(3,1,3);
plot(time,qh);
hold on
plot(t0,q0v,'ro','Markerfacecolor',[1 0 0]);
datetick();
xlim([time(1) time(end)]);


xlim([time(1) time(end)]);
ylabel('q (m^3/s)')
grid on
title('Discharge');

subplot(3,1,1);
plot(time,hh_);
hold on
plot(t0,h0v,'ro','Markerfacecolor',[1 0 0]);
datetick();
ylabel('h (m)')
grid on
title('Water level');
subplot(3,1,2);
plot(time,uh_);
hold on
plot(t0,u0v,'ro','Markerfacecolor',[1 0 0]);
datetick();
xlim([time(1) time(end)]);
ylabel('u (m/s)')
grid on
title('Velocity level');

if (pflag)
	print('img/sanggau-hadcp-discharge.eps','-depsc');
	system('LD_LIBRARY_PATH= epstopdf img/sanggau-hadcp-discharge.eps');
end

figure(8);
clf();
subplot(2,1,1);
plot(time,qh);
hold on
plot(t0,q0v,'ro','Markerfacecolor',[1 0 0]);
load('mat/discharge-sintang.mat');
daily30 = [daily30(:); daily30(:)];
%monthly = [monthly(:); monthly(:)];
%plot(datenum('15/01/2013','dd/mm/yyyy')+30.5*(0:length(monthly)-1), scale_sintang*monthly,'-g');
plot(datenum('01/01/2013')+(0:length(daily30)-1),scale_sintang*daily30,'.g');
%datetick();
xlim([time(1) time(1)+365]);
tick = datenum('01/12/2013','dd/mm/yyyy')+30.5*(0:12);
set(gca,'xtick',datenum('01/12/2013','dd/mm/yyyy')+30.5*(0:12));
set(gca,'xticklabel',datestr(tick(:),'mm/yy'));
ylabel('q [m^3/s]');
grid on
if (pflag)
	print('img/sanggau-hadcp-vs-precipitation.eps','-depsc');
	system('LD_LIBRARY_PATH= epstopdf img/sanggau-hadcp-vs-precipitation.eps');
end

% rating curve
figure(9);
clf();
qr_ = NaN(size(time));
%c = [[1; zeros(3,1)] [0;ones(3,1)] w*h0h] \ log(q0v);
%q_(fdx1) = exp(c(1) + c(3)*w*depth(fdx1));
%q_(fdx2) = exp(c(2) + c(3)*w*depth(fdx1));
or_ = 2;
c = [ones(4,1) w*[h0h(1); h0h(2:end)-or_] ] \ log(q0v);
qr_(fdx1,1) = exp(c(1) + c(2)*w*depth(fdx1));
qr_(fdx2,1) = exp(c(1) + c(2)*w*(depth(fdx2)-or_));

%c = [[1; zeros(3,1)] [0;ones(3,1)] w*[h0h(1); h0h(2:end)-or] ] \ log(q0v);
%q_(fdx1) = exp(c(1) + c(3)*w*depth(fdx1));
%q_(fdx2) = exp(c(2) + c(3)*w*(depth(fdx1)-or));

%fdx=1
%c = [ones(size(fdx)) w*h0h(fdx)] \ log(q0v(fdx));
%q_(fdx1) = exp(c(1) + c(2)*w*depth(fdx1));
%fdx = (2:4)';
%c = [ones(size(fdx)) w*h0h(fdx)] \ log(q0v(fdx));
%q_(fdx2) = exp(c(1) + c(2)*w*depth(fdx2));
plot(time,qh);
hold on;
plot(time,qr_,'r');
hold on                                                                         
plot(t0,q0v,'ro','Markerfacecolor',[1 0 0]);
grid on
fdx = [fdx1;fdx2];
R2 = 1 - mean((qr_(fdx) - qh(fdx)).^2)/var(qh(fdx),1);
disp(R2)
datetick();
xlim([time(1) time(end)]);
ylabel('q [m^/s]');
legend('HADCP','rating curve');
if (pflag)
	print('img/sanggau-rating-curve.eps','-depsc');
	system('LD_LIBRARY_PATH= epstopdf img/sanggau-rating-curve.eps');
end

