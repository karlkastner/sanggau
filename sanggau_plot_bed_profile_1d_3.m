% Mon Dec 15 13:00:56 CET 2014
% Karl Kastner, Berlin

% TODO surface level correction
% TODO time of 2016 campaign

%load ~/phd/dat/kapuas/pilot/baty/15_01_2012_Sanggau/Chart_15_1_12[15].mat
%plot3(sonardata.lat,sonardata.long,sonardata.depth);
%zlim([-15 0])

	if (~exist('reload','var') || reload)
		meta = sanggau_metadata();
		load([ROOTFOLDER,'/src/water-level/mat/water-level.mat']);
		level.val  = Kr(nid.sanggau_merged).depth;
		level.time = Kr(nid.sanggau_merged).time;
		hadcp      = sanggau_load_hadcp(meta);
		calib      = load_vadcp_discharge(meta);
		[bpilot sonardata] = sanggau_load_pilot(calib);
		b2016      = sanggau_load_bed_level_2016(true,calib.cs_(1));
		reload     = 0;
	end
	% TODO no magic numbers
	% put this into the HADCP class
	% sea RDI manual
	hadcp_range = 75;

	N      = bpilot.N;

	dhadcp = [0 meta.hadcp.redeploy.dz]; 

	t0     = [nanmedian(sonardata.time); cvec(calib.tstart); b2016.tstart];

	% reference level in 2016
	l0_201611 = nanmedian(b2016.zb - calib.zb.median);
	fprintf('l0 2016/11: %g\n',l0_201611);

	% bed level with reference to the HADCP installation level
	z_b    = [bpilot.zb calib.zb.A b2016.zb];
	% zb = zb - calib.lH;
	% dhadcp = dhadcp - calib.lH;

	% fix nan values
	dn_max = 5;
	z_b = fixnan(N,z_b,dn_max);

	% vertical reference with respect to the sea level of the hadcp
	% only the 2016 measurement has a vertical reference
	% so the bed level of all other campaings is referenced with respect to the 2016 measurement
	z_offset       = nanmedian(nanmedian(z_b(:,1:end-1),2)-z_b(:,end));
	fprintf('z_offset (gauge level - sea level): %f\n',z_offset);
	z_b(:,1:end-1) = z_b(:,1:end-1)-z_offset;
	h0             = calib.zs0-z_offset;
	level_         = level.val-z_offset;

	dhadcp = dhadcp-z_offset;


	s_ = std(z_b,[],2);
	m_ = mean(z_b,2);

	if (0)
		win = lanczoswin(-15:15);
		m_f = wmeanfilt(win,m_);
		fdx = isnan(m_);
		m_f(fdx) = NaN;

		% extrapolate to banks
		% TODO no magic numbers
		lmax = 14; % max(level_);
		m_f(1)=lmax;
		m_f(end)=lmax;
		fdx=isnan(m_f);
		m_f(fdx)=interp1(N(~fdx),m_f(~fdx),N(fdx),'pchip');
	else
		lmax = 0;
		interp = RegularizedInterpolator1();
		lambda = 100;
%		interp.lambda = lambda;
		fdx = isfinite(m_);
		interp.remesh([min(N),max(N)],length(N));
		interp.init(N(fdx),m_(fdx),[],lambda);

		%m_f = NaN(size(N));
		m_f_ = interp.vali.default;


		fdx = isfinite(z_b);
		NN = repmat(cvec(N),1,size(z_b,2));
	%	interp.lambda = size(z_b,2)*lambda;
		interp.init(NN(fdx),z_b(fdx),[],lambda);
		%interp.init(NN(fdx),z_b(fdx),length(N),[N(1) N(end)]);
		m_f = interp.vali.default;
%	figure()
%	plot([m_ m_f_ m_f])
%	pause
	end

	namedfigure(1,'Average profile');
	clf();

	fprintf('mean level %f std level %f \n',nanmedian(m_),sqrt(nanmedian(s_.^2)));
	plot(N,m_,'k');
	hold on
	% plot hadcp range
	plot(N(end)-15-[0 150],[1 1]'*dhadcp,'k','linewidth',2) 
	% plot calibration level
	hline(h0,'--');
	% plot min and max level
	hline([min(level_) max(level_)],'k-');
	%hline([min(level) max(level)]+calib.lH,'k-');
	xlim([N(1) N(end)]);
	ylim([-3 13]);
	xlabel('N (m)');
	ylabel('z (m WGS84)');
	ax1 = gca;
	addx(ax1,'xlim',[-0.5 0.5],'xlabel','\xi','ytick',[]);
	axes(ax1);

	% pruned
	namedfigure(10,'Pruned Average profile');
	clf();
	plot(N,m_f,'k');
	hold on
	% plot hadcp range
	dz  = (hadcp_range+1)*tan(deg2rad(1.3/2));
	patch([hadcp.N(1), hadcp.N(1)-hadcp_range, hadcp.N(1)-hadcp_range], ...
	      [0 -dz dz]+dhadcp(2),[1 1 1],'facecolor','none');

	% plot min and max level
	hline([min(level.val) max(level.val)]+z_offset,'k--');
%	hline([min(level.val) max(level.val)]+calib.lH,'k--');
	xlim([N(1) N(end)]);
	ylim([-3 13]);
	xlabel('N (m)');
	ylabel('z (m WGS84)');
	a1 = gca;
	title({'',''});
	a2 = addx(a1,'xlim',[-0.5 0.5],'xlabel','\xi','ytick',[],'xtick',[-0.4:0.2:0.4],'xticklabel',num2str([-0.4:0.2:0.4]'));
	axes(a1);
	set(a1,'position',get(a2,'position'));

	namedfigure(2,'Bed level profile of each campaign');
	clf();
	plot(N,z_b,'-');
	hold on
	plot(N(end)-15-[0 120],-[1 1]'*dhadcp,'r')
	legend(datestr(t0,'dd/mm/yy'))
	xlim([N(1) N(end)]);
	%ylim([-15 0])
	ylim([-3 13])

	disp('std of zb with respect to min and max level')
	std(m_)/mean(m_)
	std(m_)/(mean(m_)+(min(level.val)-max(level.val)))

	namedfigure(3,'error pdf')	
	clf();
	d = (m_-nanmean(m_));
	d2 = (m_-nanmean(m_)).^2;
%	d2 = d2(isfinite(d2));
%	n=length(d2);
%	plot(sort(d2),(1:n)/n);
	hist_man(d);
	hold on;
	median(d2)
%	s = s(isfinite(s));
%	n = length(s);
%	plot(sort(s),(1:n)/n,'r')
	title('error pdf');

	
	% bed slope
	N = double(cvec(calib.N));
	zb = double(calib.zb.mean);
	w = N(end)-N(1)+1;
	fdx = abs(N)<=w/2;
%	poly = PolyOLS([1 0 0 0 1]);
	poly = PolyOLS(8); %[1 0 0 0 0 0 1]);
	poly.fit(cvec(N(fdx)),cvec(zb(fdx)));
	'parameter'
	poly.param
	poly.serr
	poly.R2
	p = poly.predict(cvec(N));
	A = vander_1d(cvec(N(fdx)),1);
	%A = vander_1d(cvec(N(fdx)),[1 0 0 0 1]);
	%A = vander_1d(cvec(N(fdx)),1);
	c = A\cvec(zb(fdx))
	c(2)*w
	c(2)*w/mean(calib.zs0)

	namedfigure(4,'')
	A = vander_1d(cvec(N),1);
	plot(N,[-(A*c) -p]+mean(zb))

	N1 = N/(max(N)+1); % trick a bit to stay below 1
%	B1 = double(zb+midrange(calib.l0));
	B1 = double(zb+max(calib.zs0));
%	B1 = -B1/mean(B1);
%	B1 = -B1/max(B1);
%	B1 = B1 - B1(end/2);
	B1 = B1 - min(B1);
	B1 = B1/max(B1);
%	B1 = 1 - sqrt(1-N1.^2);
%	B1 = N1.^2;
%	B1 = N1.^8;
%	f = @(p) abs(n1).^p
	%f = @(n,p) 1 - (1 - abs(n).^p).^(1/p); p = 1;
	f1 = @(n,p) p(3)*(abs(n).^p(1) - 1).*(1 - p(2)*n) + 1; p1 = [1 0 1];
	[p1 resn(1)] = lsqnonlin( @(p) f1(N1,p)-B1,p1);
	f2 = @(n,p) 1 - p(2)*((1 - abs(n).^2).^(1/2)).*(1-p(1)*n); p2 = [0 1];
	%f2 = @(n,p) 1 - p(3)*((1 - abs(n).^p(1)).^(1/p(1))).*(1-p(2)*n); p2 = [1 0 1];
	[p2 resn(2)] = lsqnonlin( @(p) f2(N1,p)-B1,p2);
	%f3 = @(n,p) 1 - (1 - (n-p(1)).*(n-p(2))).^1; p3 = [0 200];
	%f3 = @(n,p) 1 - (1 - ((p(1)*n.^2 + p(2)*n + p(3)).^2)).^0.5; p3 = [0 1 0];
	%f3 = @(n,p) 1 - p(4)*(1 - (p(1)*n.^2 + p(2)*n + p(3)).^p(5)).^(1/p(5)); p3 = [0 -0.5 0 1 2];
	%f3 = @(n,p) 1 - (1 - (n.^2 - p(1)*n a(n-p(1)).*(n-p(2))).^1; p3 = [0 200];
%	f3 = @(n,p)  p(3)*n(n < p(1)) + *p(4)*(n(n >= p(1)) & n(n < p(2))) + n(n >= p(2));
	f3 = @(n,p) pwp(n,[-1;p(1); p(2); 1],[1; p(3); p(4); 1]); p3 = [-150/300, 150/300, 0, 0];
%	f3 = @(n,p) pwp(n,[-1;p(1); p(2); p(3); 1],[1; p(4); p(5); p(6); 1]); p3 = [-150/300, 0, 150/300, 0, 0, 0];
	f3 = @(n,p) pwp(n,[-1; p(1); 1],[p(2); p(3); p(4)]); p3 = [-150/300, 0, 150/300, 0, 0, 0];
	[p3 resn(3)] = lsqnonlin( @(p) f3(N1,p)-B1,p3);
	resn
	
	x=linspace(-1,1,100)';
	figure(5)
	%plot(x,bsxfun(@power,abs(x),[1,2,3,4,5]))
	plot(N,[B1,f1(N1,p1),f2(N1,p2),f3(N1,p3)]);
	[p1 p2]
	rad2deg(grad(@(x) f1(x,p1), 0.99)*max(calib.zs0)/330)
	rad2deg(grad(@(x) f2(x,p2), 0.99)*max(calib.zs0)/330)
	legend('original','polynomial','elliptic','trapezoidal')
	
% for circle: 3.7, for quad : 2 (off course)

	namedfigure(6,'one parameter fits');
	clf
	r = [];
	B1 = B1 - 1;
	plot(N1,B1); hold on
	f = @(n,p) p*(n.^2 - 1); p = 1;
	[p rn(1)] = lsqnonlin(@(p) f(N1,p) - B1, p);
	plot(N1,f(N1,p));
	f = @(n,p) p*(1 - n.^2).^(0.5); p = 1;
	[p rn(2)] = lsqnonlin(@(p) f(N1,p) - B1, p);
	plot(N1,f(N1,p));
	f = @(n,p) p*ones(size(n)); p = mean(B1);
	[p rn(3)] = lsqnonlin(@(p) f(N1,p) - B1, p);
	plot(N1,f(N1,p));
	f = @(n,p) p(1) + p(2)*n; p = [1 0];
	[p rn(4)] = lsqnonlin(@(p) f(N1,p) - B1, p);
	{'orig','para','ellipic','rect','linear rect (2 param)'}
	rn
	plot(N1,f(N1,p));
	legend('orig','para','ellipic','rect','linear rect (2 param)');

	% consecutive fits
	figure(7)
	clf
	plot(N1,B1)
	hold on
	resn = [];
	for idx=1:7
		f = @(n,p) pwp(n,[-1 p(1:idx) 1]',[1 p(idx+1:end) 1]');
		p = [2*((1:idx)/(idx+1)-0.5) ones(1,idx)];
		[p resn(idx)] = lsqnonlin(@(p) f(N1,p) - B1, p);
		plot(N1,f(N1,p));
	end	
	resn

	fdx = abs(N) < 250;
	reg = PolyOLS(1,false);
	reg.fit(N(fdx),zb(fdx)+midrange(calib.zs0));

	z=linspace(min(m_f),max(m_f));
	id=[];
	for idx=1:length(z);
		d=max(0,z(idx)-m_f);
		A=sum(d);u=d.^(1/2);
		U=sum(d.^(3/2))/A;
		id(idx,1)=max([NaN find(u>U,1,'FIRST')]);
		id(idx,2)=max([NaN find(u>U,1,'LAST')]);
	end
	fdx=isfinite(id);
	nmin = nan;
	nmin(fdx(:,1),1) =  -(N(1)-N(id(fdx(:,1),1)));
	nmin(fdx(:,2),2) =    N(end)-N(id(fdx(:,2),2));
	plot(z,nmin);
	vline(min(level_));
	xlabel('low flow water level');
	ylabel('minimum required profiling range');
	legend('left bank','right bank (hadcp installed)');

		fprintf('Bankful area: %f\n',csarea(lmax,m_f,1));
	fprintf('Bankful perimeter: %f\n',csperimeter(lmax-sqrt(eps),m_f,1));
	fprintf('Bankful hydraulic radius: %f\n',csradius(lmax-sqrt(eps),m_f,1));
	fprintf('Thalweg level: %f\n',min(m_f));
	fprintf('HADCP level before and after deployment: %f %f\n',dhadcp-min(m_f));
	fprintf('HADCP depth below top of banks: %f %f\n',lmax-dhadcp);

	res = bsxfun(@minus,z_b,mean(z_b,2));
	inner = (N>-250 & N < 250);
	rmse = rms(flat(res(inner,:)))
	k = sum(inner)-1;
	[ac av] = autocorr_man4(res(inner,:),k,[],true,true);
	'sac'
	sum(ac)
	figure(1e3)
	clf
	subplot(2,2,1)
	plot(av)
	subplot(2,2,2)
	m = mean(av,2);
	s = std(av,[],2);
	s = s/sqrt(size(av,2));
	errorlines(1:length(m),m,m-s,m+s);
	hold on
	%$serr(av')';
	%errorlines(1:length(m),NaN*m,m-s,m+s);
	%plot(median(av,2))
	%plot(hodges_lehman(av,2))
	c = lsqnonlin(@(c) w.*flat(bsxfun(@minus,c(1)*acfar1(c(2),k,(0:k-1)',true),av)),[av(1),ac(2)]);
	hold on
	c
	plot(c(1)*acfar1(c(2),k,(0:k-1)',true));
	av=mean(av,2);
	c = lsqnonlin(@(c) w.*(c(1)*acfar1(c(2),k,(0:k-1)',true)-av),[av(1),ac(2)]);
	hold on;
	c
	plot(c(1)*acfar1(c(2),k,(0:k-1)',true));

	subplot(2,2,3)
	plot(mean(ac,2))
	w = 1;
	ac = av/av(1);
	rho = lsqnonlin(@(rho) w.*(acfar1(rho,k,(0:k-1)',true)-ac),ac(2))
	hold on
	plot(acfar1(rho,k,(0:k-1),true))
	hold on
	%plot(acfar1(rho,k,(0:k-1),false))
	plot(acfar1(rho^0.5,k,(0:k-1),true))
	


	L_A = -1/log(rho)
	width = range(N)
	s_A = width*[std(mean(z_b(inner,:))), std(mean(z_b))]./mean(calib.area)
	      width*[rmse, rms(res(:))]*sqrt(L_A/width).*1./mean(calib.area)

	if (exist('pflag','var') && pflag)
		ps = [2 2.5 3];
		for idx=1:length(ps)
			pdfprint(1,'img/sanggau-bed-level-1d-averaged',ps(idx));
			pdfprint(10,'img/sanggau-bed-level-1d-averaged-pruned',ps(idx));
		end

		%pdfprint('img/sanggau-zb-1d'); %-' method '-' num2str(worder) '.eps']);
		fid = fopen('mat/sanggau-bottom-1d.csv','w');
		fprintf(fid,'%7s; %7s; %7s;\n','N (m)','D (m)','std(D) (m)');
		fprintf(fid,'%7.2f; %7.2f; %7.2f;\n',[cvec(N), cvec(m_), cvec(s_)]');
		fclose(fid);
		pflag=0;
	end

