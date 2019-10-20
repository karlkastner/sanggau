% Sun Jan 25 21:15:40 CET 2015
% Karl Kastner, Berlin

%
% time derivative of stage in Sanggau
	if (~exist('pflag','var'))
		pflag = false;
	end
	T = 2;

	meta = sanggau_metadata();	
	% load hadcp data
	if (~exist('reload','var') || reload)
		[hadcp hlevel offset] = sanggau_load_hadcp(meta);
		load([ROOTFOLDER,'/src/water-level/mat/water-level.mat']);
		level.val  = K(nid.sanggau_merged).slot.depth;
		level.time = K(nid.sanggau_merged).slot.time;
		calib = load_vadcp_discharge(meta.filename.discharge,level);
		reload = 0;
	end
	dt_max = 3/24;
	level.val = fixnan(level.time,level.val,dt_max,'linear',NaN);

	level.val = level.val + calib.lH; % meta.lH
	calib.l0 = calib.l0 + calib.lH;

	% compute derivative in time
	dt    = diff(level.time);
	dH    = diff(level.val);
	dH_dt = dH./dt;
	% filter daily
	dn    = ceil(T./median(dt));
	fH    = trifilt1(level.val,dn);
	dH_dt = diff(fH)./dt;
%	dH_dt = meanfilt1(double(dH_dt),dn);
	%dH_dt = medfilt1(double(dH_dt),dn);
	time = level.time(1:end-1);
	fdx  = isfinite(dH_dt);
	dH_dt0 = interp1(time(fdx),dH_dt(fdx),calib.t0,'nearest');

	figure(1);
	clf();
%	subplot(2,1,1)
%	plot(level.time(1:end-1),dH_dt);
%	hold on
%	plot(calib.t0,dH_dt0,'ro','Markerfacecolor','r');
%	subplot(2,1,2)
%	plot(level.time,level.val);
%	hold on;
%	plot(calib.t0,calib.l0,'ro','Markerfacecolor','r')
	[ax h1 h2] = plotyy(time,dH_dt,level.time,level.val);
	linkaxes(ax,'x');
	hold(ax(1),'on');
	hold(ax(2),'on');
	plot(ax(1),calib.t0,dH_dt0,'ro','Markerfacecolor','r');
	violet = hsv2rgb([32 255 196]/255);
	plot(ax(2),calib.t0,calib.l0,'o','color',violet,'Markerfacecolor',violet);
	ylim(ax(1),[-0.8 0.8]);
	ylim(ax(2),[0 16]);
	datetick(ax(1));
	set(ax(2),'xtick',[]);
	set(ax(1),'ytick',-1:0.2:1');
	set(ax(2),'ytick',0:2:16);

	% CDF
	figure(2);
	clf();
	n = sum(fdx);
	cdf_ = sort(dH_dt(fdx));
	% make unique
	cdf_ = cdf_ + (0:n-1)'*1e-7;
	subplot(2,2,1)
	plot(cdf_,(1:n)/n);
	cdf0 = interp1(cdf_,(1:n)/n,dH_dt0,'linear');
	hold on
	plot(dH_dt0,cdf0,'ro','markerfacecolor','r');
	ylabel(ax(1),'dH/dt (m/s)');
	ylabel(ax(2),'H (m)');

	subplot(2,2,2);
	cla
	fdx = isfinite(dH_dt);
	q = quantile(dH_dt,[0.01 0.99]);
	
	x = linspace(q(1),q(2));
	s = 64*std(dH_dt(fdx))/sqrt(sum(fdx));
	sd = nanstd(dH_dt);
	sk = skewness(dH_dt(fdx))
	y = [];
	y(:,1) = kernel1d(x,dH_dt(fdx),s);
	y(:,2) = normpdf(x,0,sd);
	y(:,3) = skewpdf(x,0,sd,sk);
	rms(bsxfun(@minus,y(:,2:end),y(:,1)))
	plot(x,y);
	vline(dH_dt0)


	% TODO exclude banks
	% infer friction slope from vadcp data
	% TODO, move to discharge class
%	ln_z0 = median(calib.ln_z0_A,2);
if (0)
	clear S c sqrtS
	ln_z0 = nanmedian(calib.ln_z0_A(:));
	ln_z0_A = real(calib.ln_z0_A);
	for idx=1:length(calib.dis)
		dis = calib.dis(idx);
		H   = dis.bgrid.val(:,1);	
		U   = dis.U;
		c(idx,1) = 1.5*nanmedian(U);
%		H   = nanmedian(H);
%		U   = nanmedian(U);
		sqrtS(idx,:,1) = Constant.KAPPA*dis.U./(sqrt(Constant.g*H).*(log(H)-ln_z0_A(:,idx)-1));
		sqrtS(idx,:,2) = Constant.KAPPA*dis.U./(sqrt(Constant.g*H).*(log(H)-ln_z0));
	end
	'std_sS'
	squeeze(nanstd(sqrtS,[],2))
	'mean_sS'
	sqrtS = squeeze(nanmean(sqrtS,2))
	
%	sqrtS = nanmedian(sqrtS);
	figure(3);
	clf
	subplot(2,2,2);
	S = sqrtS.^2;
	% expected (intercep unknown)
	Se = median(S(:)) + 1./c.*(dH_dt0/86400);
	s = 1e5;
	plot(dH_dt0,s*[S Se],'o'); hold on
	plot(dH_dt0,s*Se,'r');
	xlabel('dH/dt (m/d)');
	ylabel('S 10^{-5}');
end

	limits(dH_dt0)
	if (pflag)
		pdfprint(1,'img/sanggau-dh-dt-vs-time.pdf');
		pdfprint(2,'img/sanggau-dh-dt-cdf.pdf');
		pdfprint(3,'img/sanggau-energy-slope-vs-dH-dt.pdf');
	end

