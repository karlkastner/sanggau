% 2015-09-27 15:44:39.276071119 +0200
% Karl Kastner, Berlin

	if (~exist('reload','var') || reload)

	addpath ../sanggau

	% preload water level
	meta = sanggau_metadata();
	% preload water level
	level = sanggau_load_water_level(meta);
	% preload vadcp calibrarion discharge
	meta.filename.discharge
	calib = load_vadcp_discharge(meta.filename.discharge,level);

	% level is whith respect to mean bottom elevation
	bottom = smooth(calib.bottom.median,20)-calib.lH;

	% area : fit of polynomial
	afunc_ = @(h) afunc(h,bottom,calib.dw);

	% get HADCP samples during calibration +- 1h
	hadcpname = '/home/pia/large/phd/src/hadcp/mat/hadcp-sanggau-squeezed-1800-preloaded.mat';
	load(hadcpname);

	% from redeployment onward, the instrument was deployed upside down, that swaps beam 1 and 2
	% TODO, do this in preload already
	hadcp.to_beam();
	fdx = hadcp.dat.FileNumber > 5;
	v1 = hadcp.dat.VEL(:,fdx,1);
	hadcp.dat.VEL(:,fdx,1) = hadcp.dat.VEL(:,fdx,2);
	hadcp.dat.VEL(:,fdx,2) = v1;
	hadcp.convert_raw_velocity();

	reload = 0;
	end % reload



	% beam 3 was partially blocked, so use only beam 1 and 2	
	mode = 12;
	hadcp.to_earth(mode);

	% rotate hadcp earth vel to cross section direction
	hadcp.to_cs(calib.dir);


	% the instrument was redeployed and became misagligned, therefore the
	% periods are calibrated individually
	% 

	index = {   [1, 3216], [1];
		 [3216, 5547], [2]
                 [5547, 11000], [3,4];
                 [11000, length(hadcp.time)], [5] };
	C = [];
	for idx=1:size(index,1)
	cdx = index{idx,2};
	range = index{idx,1}(1):index{idx,1}(2);
	
	A0 = afunc_(calib.h0(cdx));
	Q0 = calib.cs.q0(cdx);
	% first 110 are good in 99.4% of all samples
	nmax = 100;
	u0 = interp1(hadcp.time,hadcp.velocity.cs(:,:,2)',calib.t0(cdx),'linear')';
	u0 = u0(1:nmax,:);
	% calibrate the ivm
	c = (A0.*cvec(mean(u0))) \ Q0
	
	% compute discharge time series
	% TODO fixing invalid values should not be done here
	fdx = isfinite(level.val);
	l = interp1(level.time(fdx),level.val(fdx),hadcp.time)+calib.lH;
	A(range,1) = afunc_(l(range));
	u  = hadcp.velocity.cs(1:nmax,:,2);
	%ue = hadcp.velocity.earth(1:nmax,:,2);
	mean_u = mean(u);
	Q(range,1)  = c*A(range,1).*cvec(mean_u(range));
	C(idx) = c;
	end

	% scale fix for period 2 (misalignment)
	range = index{2,1}(1):index{2,1}(2);
	Q(range,1) = C(4)*A(range,1).*cvec(mean_u(range));
%	Q(range,1) = 9.4825/11.3791*C(3)*A(range,1).*cvec(mean_u(range));
	

	% calibrate the bins individually
%	c_i = A0*u0 \ Q0;
	c_i = Q0 ./ (A0*u0);
	% individual discharge estimate
	Q_i = bsxfun(@times,rvec(A),bsxfun(@times,c_i,u));
	% difference to joint calibration
	d_i = bsxfun(@plus, Q_i,-rvec(Q));
%	Qe = -c*A.*cvec(rmse(ue));

	% load rating curve discharge
	rc = load('../discharge/mat/sanggau-stage-discharge.mat');
	qma = interp1(rc.time,rc.qma,hadcp.time);

	namedfigure(1,'Time series of discharge');
	clf();
	plot(hadcp.time,[Q qma]/1e3,'.');
	datetick('x','dd/mm/yy');
	hold on;

if (0)
	namedfigure(2,'HADCP velocity orthogonal to cross section')
%	surface(hadcp.time,1:149,hadcp.velocity.cs(:,:,2),'edgecolor','none')
%	caxis([-2 0]);
end


 	% difference between rating curve and HADCP
	d  = (qma - Q);
	d2 = d.^2;
	d2(1:12000) = NaN;
	d2(end-48:end) = NaN;

	% largest 4 differences with at least 1 week in between
	Ma = [];
	Mi = [];
	% find time slots with largest deviation
	for idx=1:4
		[mv mdx] = max(d2);
		Ma(idx) = mdx;
		d2(mdx-(48*7):mdx+48*7) = NaN;

		[mv mdx] = min(d2);
		Mi(idx) = mdx;
		d2(mdx-(48*7):mdx+48*7) = NaN;
	end
	
	namedfigure(4,'velocity profile during maximum deviation')
	clf
	namedfigure(5,'Backscatter during maximum deviation');
	clf
	namedfigure(6,'velocity profile during minimum deviation')
	clf
	for idx=1:4
		figure(4)
		subplot(2,2,idx)
		plot(hadcp.velocity.cs(:,Ma(idx),2));
		hold on
		plot(hadcp.velocity.earth(:,Ma(idx),4),'r');
		figure(5)
		subplot(2,2,idx)
		E = hadcp.E(:,Ma(idx),1:3);
		plot(squeeze(E))
		figure(1);
		plot(hadcp.time(Ma(idx)),Q(Ma(idx))/1e3,'or');

		figure(6);
		subplot(2,2,idx)
		plot(hadcp.velocity.cs(:,Mi(idx),2));
		hold on
		plot(hadcp.velocity.earth(:,Mi(idx),3),'r');
		figure(7);
		subplot(2,2,idx)
%		plot(hadcp.velocity.earth(:,Mi(idx),4),'r');
		E = hadcp.E(:,Mi(idx),1:3);
		plot(squeeze(E))
		figure(1);
		plot(hadcp.time(Mi(idx)),Q(Mi(idx))/1e3,'og');
		
	end

	namedfigure(101,'Median beam velocities')
	plot([nanmedian(hadcp.velocity.beam(:,12000:end,1),2), ...
             -nanmedian(hadcp.velocity.beam(:,12000:end,2),2), ...
             -nanmedian(hadcp.velocity.beam(:,12000:end,3),2)])

	% TODO do the same for ton's method
	namedfigure(10,'CDF of Q_rc - Q_hadcp');
	sd=sort(d(isfinite(d(12000:end))));
	n=length(sd);
	plot(sd,(1:n)/n);
	q=quantile(sd,normcdf([-2:1:2]));
	hold on;
	p=(1:n)/n;
	sd_=qstd(sd);
	mu_=median(sd);
	plot((sd_*norminv(p)+mu_), p,'r');
	xlim([-3 3]*sd_)

	% internal consistency
	namedfigure(11,'Internal consistency of the HADCP discharge estimate');
	d_i = d_i(:,12000:end);
	subplot(2,2,1)
	plot(qstd(d_i))        
	subplot(2,2,2)
	plot([nanmedian(d_i'); qstd(d_i')]')
	subplot(2,2,3)
	plot(quantile(d_i',[0.16 0.84])')	

	namedfigure(12,'HADCP vs RC')
	% hysteresis is visible, but small (less 500-1000m^3/s)
	% not 1:1
%	plot(l,d,'-')
	clf;
	plot(qma/1e3,Q/1e3,'.-');
	hold on;
	plot([2 10],[2 10],'k');
	qma = double(qma);
	Q   = double(Q);	
	for idx=1:48*7:length(qma)
		t_str = datestr(hadcp.time(idx),'dd/mm');
		text(qma(idx)/1e3,Q(idx)/1e3,t_str);
	end

	range=[14643:15501];
	plot(qma(range)/1e3,Q(range)/1e3,'r');
	range=[11678:13006];
	plot(qma(range)/1e3,Q(range)/1e3,'r');
	
