% Tue Mar 24 19:53:52 CET 2015
% Karl Kastner, Berlin

if (~exist('reload','var') || reload)
meta = sanggau_metadata();
% get HADCP velocity
[hadcp hlevel]= sanggau_load_hadcp(meta,[],1800);
hl.time = hadcp.time;
hl.val  = hlevel;
calib = load_vadcp_discharge(meta.filename.discharge,hl);
% rotate hadcp velocity do the cross section
hadcp.velocity.cs = rotvel(calib.dir, hadcp.velocity.earth);

% load rating curve discharge
rc = load('../discharge/mat/sanggau-stage-discharge.mat');
reload = 0;
end
% TODO quick fix
calib.bottom.median(end) = calib.bottom.median(end-1);

figure(1);
clf();
subplot(2,2,1)
% select time spans of the HADCP to use
plot(hadcp.time,hlevel);
%plot(hadcp.time,hadcp.idepth_m);
hline(0);
hline(-meta.hadcp.redeploy.d);

T = 1e5*[	7.358279430292825,   7.358563805292825;
		7.359471439499999,   7.360159772833333];

% select valid cells
fdx = (hadcp.time > T(1,1) & hadcp.time < T(1,2)) | (hadcp.time > T(2,1) & hadcp.time < T(2,2));

time  = hadcp.time(fdx);
%level = hlevel - calib.lH;
level = hlevel(fdx);
U = hadcp.velocity.cs(:,fdx,1);


load('mat/fgood.mat')
u0 = U(:,fgood);
level0 = level(fgood);

med = nanmedian(u0)';
mea = nanmean(u0)';
order = [1 1];
reg = PolyOLS(order);
s = [];
c = [];
for idx=1:size(u0,1)
	%A = vander_1d(med,1);
	reg.regress(mea,u0(idx,:)');
	c(idx,:) = reg.param;
	s(idx,:) = reg.params;
end
N = 1:size(u0,1);

namedfigure(100,'coefficients');
errorlines(N,c,s)
% predicted
up = (vander_1d(mea,order)*c')';
D = up - u0;
imagesc(D);

n = 1000;
for idx=1:size(D,1)
	acf1(idx,:) = autocorr(D(idx,:),1000);
end

namedfigure(101,'autocorrelation function (time)');
subplot(2,1,1)
imagesc(acf1);
subplot(2,1,2)
ma = nanmean(acf1,1);
plot(ma);
ma(2)

for idx = 1:size(D,2)
	acf2(:,idx) = autocorr(D(:,idx),148);
end

namedfigure(102,'Autocorrelation function (space)')
subplot(2,1,1)
imagesc(acf2)
ma = nanmean(acf2,2);
subplot(2,1,2)
plot(ma,'.');
ma(2)
se = nanstd(D(:));
sdat = nanstd(u0(:));
R2 = 1 - (se/sdat)^2;

fprintf(1,'error s^2 %f\n',se);
fprintf(1,'data  s^2 %f\n',sdat);
fprintf(1,'R2        %f\n',R2);

% print velocity profile for different stages
n = 4;
mea = mean(u0)';
%q = [-inf, quantile(level,(1:n-1)/n), inf];
%q = quantile(level,[1/n,1-1/n]);
q = quantile(mea,[1/n,1-1/n]);
q = [-inf, q(1)+(q(2)-q(1))*(0:n-1)/(n-1), inf]
%q = [-inf, quantile(mea,(1:n-1)/n), inf];
clear u_;
for idx=1:length(q)-1
	%fdx = level0 > q(idx) & level0 < q(idx+1);
	fdx = mea > q(idx) & mea < q(idx+1);
	u_(idx,:) = mean(u0(:,fdx),2);
end
figure(10);
plot(u_);
plot(bsxfun(@times,u_',1./mean(u_')));

