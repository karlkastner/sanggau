% 2015-01-09 17:53:41.747242217 +0100
% Karl Kastner, Berlin

	hudim = 1;
	dt_plot = 7;

	method = {'empirical', 'empirical','theoretical'};
	worder = [0,1,1];
	
	% load rating curve
	if (reload)
		clear a
		rating = load('../discharge/mat/sanggau-stage-discharge-2014-12-16.mat');
		for mdx=1:length(method)
			a(mdx) = load(['mat/sanggau-hadcp-calibration-' method{mdx} '-' num2str(worder(mdx)) '.mat']);
		end
		reload = 0;
	end

	fdx = find(isfinite(rating.discharge),1,'first');
	ldx = find(isfinite(rating.discharge),1,'last');
%	tlim = [a(1).hadcp.time(1) a(1).hadcp.time(end)];
	tlim = [rating.time(fdx) rating.time(ldx)];

	% joint cw
	figure(100);
	clf()
	cw =a(2).bin.cw;
	cw(:,2) = -cw(:,2);
	% a21 is ident with a11
	fdx = find(cw(:,1),1);
	cw(1:fdx-1,:) = NaN;
	plot([cw]);
	ylim([0 1500]);
	xlabel('m');
	ylabel('m');	
	legend('c_{q0}','-c_{q1}');
	pdfprint('img/sanggau-specific-discharge-weights.eps');

	figure(101)
	clf();
	%subplot(2,1,1);
	%plot(hadcp.time,[qh qh_all]);
	%hold on;
	%plot(hadcp.time,[qh-qherr qh+qherr],'b--');
	%datetick()
	%subplot(2,2,3);
	subplot(2,1,1);
	%[ax h1 h2] = 
	plot([NaN NaN NaN]); hold on
	plot(a(1).hadcp.time,[a(1).qh a(2).qh]/1e3,'.','Markersize',2);
	hold on
	skip = round(1/(4*median(diff(rating.time))));
	plot(rating.time(1:skip:end),rating.discharge(1:skip:end)/1e3,'.r','Markersize',2);
	plot(a(1).calib.t0,a(1).calib.q0/1e3,'ok','markerfacecolor','k');
	ylim([0 14]);
	ylabel('q (10^3 m^3/s)');
	%ylabel(ax(2),'s_q (10^3 m^3/s)');
	set(gca,'ytick',(0:2:14));
	datetick('x','mmm');
	set(gca,'xtick',monthspace(tlim(1),tlim(2)));
	xlim(tlim);
	legend('HADCP constant','HADCP linear','rating curve');
	title('Discharge')
	subplot(2,1,2);
	title('standard error')
	plot(a(1).hadcp.time,[a(1).qherr a(2).qherr]/1e3,'.','markersize',2);
	set(gca,'ytick',(0:0.2:1.4));
	datetick('x','mmm');
	set(gca,'xtick',monthspace(tlim(1),tlim(2)));
	ylabel('q (10^3 m^3/s)');
	xlim(tlim);
	ylim([0 1]);
	title('Standard error \sqrt{(\bar Q_i - Q_i)^2/n}')
	pdfprint('img/sanggau-discharge-time-series.eps');
	
	for mdx=1:length(method)

	%load(['mat/sanggau-hadcp-calibration-' method{mdx} '-' num2str(worder(mdx)) '.mat']);
	unpack_struct(a(mdx));

	Nz = calib.dis(1).z0grid.cX1();
	Nb = calib.dis(1).bgrid.cX1();
	Nq = calib.dis(1).vgrid.cX1();

	% level to depth average conversion
	lH = nanmean(nanmean(calib.bottom_A));

%{
	% plot roughness length and bottom profile
	namedfigure(1,'Roughness length and bottom');
	clf();
	subplot(2,1,1)
	plot(Nz,ln_z_0_A_);
	xlim([min(Nz) max(Nz)]);
	title('Roughness length');
	subplot(2,1,2);
	hold on
	plot(Nb,-calib.bottom_A);
	title('Bottom');
	xlim([min(Nb) max(Nb)]);
%}
	namedfigure(1,'Stage and discharge time series');
	clf();
	[ax h1 h2] = plotyy(hadcp.time,qh/1e3,hadcp.time,hlevel+lH);
%	[ax h1 h2] = plotyy(hadcp.time,qh/1e3,hadcp.time,level.val+lH);
	ylim(ax(1),[0 14]);
	ylim(ax(2),[0 14]);
	set(h1,'LineStyle','.');
	set(ax(1),'ytick',(0:2:14));
	set(ax(2),'ytick',(0:2:14));
	linkaxes(ax,'x');
	set(ax(2),'xtick',[])
	datetick();
	xlim(tlim)
	ylabel(ax(2),'H (m)');
	ylabel(ax(1),'q (10^3 m^3/s)');
	hold on
	plot(ax(1),calib.t0,calib.q0/1e3,'or','markerfacecolor','r');
	hold on
	% mark bad samples
	qf = medfilt1(double(qh),48);
	qout = 4*quantile(abs(qh-qf),0.5);
	odx = find(abs(qh-qf) > qout);
	plot(hadcp.time(odx),qh(odx)/1e3,'.g');
	grid on

	namedfigure(2,'HADCP discharge vs. rating curve');
	clf;
	%subplot(2,1,1);
	%plot(hadcp.time,[qh qh_all]);
	%hold on;
	%plot(hadcp.time,[qh-qherr qh+qherr],'b--');
	%datetick()
	%subplot(2,2,3);
	[ax h1 h2] = plotyy(hadcp.time,qh/1e3, hadcp.time,qherr/1e3);
	hold on
	plot(ax(1),rating.time,rating.discharge/1e3,'r');
	plot(ax(1),calib.t0,calib.q0/1e3,'or','markerfacecolor','r');
	ylim(ax(1),[0 14]);
	ylim(ax(2),[0 1.4]);
	set(h1,'LineStyle','.');
	set(h2,'LineStyle','.');
	ylabel(ax(1),'q (10^3 m^3/s)');
	ylabel(ax(2),'s_q (10^3 m^3/s)');
	set(ax(1),'ytick',(0:2:14));
	set(ax(2),'ytick',(0:0.2:1.4));
	linkaxes(ax,'x');
	set(ax(2),'xtick',[])
	datetick('x','mmm');
	set(gca,'xtick',monthspace(tlim(1),tlim(2)));
	xlim(tlim);
	'bias median'
%	nanmedian(qh-qh_all)
	% add rating curve

	
	% plot 2D hadcp velocity timeseries
	namedfigure(3,'Depth averaged streamwise velocity 2D');
	imagesc(bin.Ubar);
	caxis([-2 2]);
	title('depth averaged hadcp velocity');

	% plot of depth averaged HADCP and VADCP velocity 
	namedfigure(4,'Streamwise HADCP and VADCP depth average velocity vs range');
	% extract hadcp velocity during calibration
	% TODO, provide this by the calibration script
	Uh0 = interp1(hadcp.time,bin.Ubar',calib.t0,'linear')';
	sig = sign(nanmedian(Uh0(:)));
	% extract vadcp velocity in hadcp range
	Uv = interp1(Nq(:),calib.U,Nh,'linear');
	for idx=1:size(Uv,2)
		subplot(2,2,idx);
		%plot([u(1:149,idx) U0(:,idx)]);
		plot(Nh,sig*[Uv(:,idx) Uh0(:,idx)]);
		title(datestr(calib.t0(idx),'dd/mm/yyyy'));
		ylim([0 1.3]);
		xlim([min(Nh) max(Nh)]);
		ylabel('u (m/s)');
		xlabel('N (m)');
		grid on
	end
	
	namedfigure(5,'Streamwise HADCP and VADCP velocity at HADCP depth');
	% TODO provide d0 by the calibration scrip
	hd  = bin.h - bin.z;
	lh  = bin.h - bin.z - bin.level;
	lh0 = interp1(hadcp.time,lh',calib.t0,'linear')';
	dh0 = interp1(hadcp.time,hd',calib.t0,'linear')';
	uh0 = interp1(hadcp.time,hvel(:,:,hudim)',calib.t0,'linear')';
	sig = sign(nanmedian(uh0(:)));

	for idx=1:length(calib.dis)
		subplot(2,2,idx);
%		NN = discharge(idx).vgrid.cXX1()*discharge(1).cs.width;
%		ZZ = discharge(idx).vgrid.cXX2();

		% streamwise hadcp velocity
		uu = double(flipud(calib.dis(idx).vgrid.val(:,:,1)'));
		% n-coordinate of VADCP vel mesh bin centres
		Nv = double(calib.dis(idx).vgrid.cX1());
		% d-coordinate of VADCP vel mesh bin cenctres
		Dv = double(flipud(-calib.dis(idx).vgrid.cX2'));
%		Zv = dischargedischarge(idx).vgrid.cX2();
%		uv = interp2(NN,double(ZZ),double(uu),Nh,double(hz0(:,idx)));

%		uv = interp1(Nv(:),uu',Nh,'linear')';
%		uv = interp1(Dv(:),uv,dh0(1,idx),'linear')';
		uv(:,idx) = interp2(Nv,Dv,uu,Nh,dh0(:,idx),'linear');

		% interpolate vadcp velocity to hadcp depth
		plot(Nh,sig*[uv(:,idx) uh0(:,idx)]);
		title(datestr(calib.t0(idx),'dd/mm/yyyy'));
		ylim([0 1.3]);
		xlim([min(Nh) max(Nh)]);
		ylabel('u (m/s)');
		xlabel('N (m)');
		grid on
	end

	namedfigure(6,'CDF of u_h - u_v');
	clf()
%	subplot(2,2,1);
%	plot(u);
%	subplot(2,2,3);
%	plot(U0);
%	subplot(2,2,4)
%	n = length(uh0);
	cc = colormap('lines');
	for idx=1:size(uh0,2)
		d = uh0(:,idx) - uv(:,idx);
		fdx = find(isfinite(d));
		n = length(fdx);
		plot(sort(d(fdx)), (1:n)/n,'color',cc(idx,:));
		hold on
	end
%	hist(uh0 - uv);
	title('CDF of u_hadcp - u_vadcp')
	%plot(u(1:149,:),U0,'.');

	namedfigure(7,'Discharge weigths');
%	subplot(2,1,1);
	imagesc(bin.c_q);
%	title('Discharge weights');

	namedfigure(71,'Weigths');
	clf();
	c_u = nanmean(bin.c_u,2);
	c_q = nanmean(bin.c_q,2);
	c_qu = c_u.*c_q;
	s_u = nanstd(bin.c_u,[],2);
	s_q = nanstd(bin.c_q,[],2);
	s_qu = nanstd(bin.c_u.*bin.c_q,[],2);
%	imagesc( );
	subplot(3,1,1)
	plot(c_u);  hold on;
	plot(c_u*[1 1]+s_u*[-1 1],'b--');
	title('c_u');
	subplot(3,1,2)
	plot(c_q);  hold on;
	plot(c_q*[1 1]+s_q*[-1 1],'b--');
	title('c_q');
	subplot(3,1,3)
	plot(c_u);  hold on;
	plot(c_qu*[1 1]+s_qu*[-1 1],'b--');
	title('c_u c_q');
%	namedfigure(72,'Velocity weigths');
%	imagesc(nanmean(bin.c_q.*bin.c_u,2));
%	namedfigure(73,'Velocity weigths');
%	imagesc(nanmean(bin.c_q.*bin.c_u,2));

	namedfigure(8,'Discharge weight coefficients');
%	subplot(2,1,2);
	plot(Nh,bin.cw);

%	namedfigure(10,'Earth velocity 2D');
	namedfigure(9,'Discharge vs. stage');
	clf();
	% interpolate to 30min interval
	ti = hadcp.time(1):1/48:hadcp.time(end);
	ti_ = interp1(hadcp.time,hadcp.time,ti,'constant');
	gap = abs(ti-ti_);
	fdx = find(gap > 1/24);
	li = interp1(hadcp.time,hlevel,ti,'linear') + lH;
	%li = interp1(hadcp.time,level.val,ti,'linear') + lH;
	qi = interp1(hadcp.time,qh,ti,'linear');
	qf = medfilt1(double(qi),48);
	qf(fdx) = NaN;
	qi(fdx) = NaN;
%	t=round(hadcp.time(1):7:hadcp.time(end));
	t=round(hadcp.time(1):dt_plot:hadcp.time(end));
	lt=interp1(ti,li,t);
	qt=interp1(ti,qf,t);
	plot(li,qf/1e3,'.-');
	hold on
 	plot(calib.l0+lH,calib.q0/1e3,'ro','markerfacecolor','r')
	text(double(lt),qt/1e3,datestr(t));
	title('stage discharge relation');
	ylabel('q (10^3 m/s)');
	xlabel('H (m)');
	xlim([0 13]);
	ylim([0 10]);
	grid on

	for idx=1:4
		subplot(4,1,idx)
		imagesc(hadcp.velocity.earth(:,:,idx))
		caxis(1.5*[-1 1]);
	end

	namedfigure(11,'Earth velocity along range');
	clf();
	subplot(2,1,1)
	fdx=1:1500;
	%fdx=11150:11960;
	%plot([mean(vel0(:,fdx,1),2) mean(hadcp.velocity.earth(:,fdx,1),2)],'.-');
	plot([mean(hadcp.velocity.earth(:,fdx,1),2)],'.-');
	subplot(2,1,2)
	plot([mean(hadcp.velocity.earth(:,fdx,2),2)],'.-');
	%plot([mean(vel0(:,fdx,2),2) mean(hadcp.velocity.earth(:,fdx,2),2)],'.-');

%	quant = 0.98;
%	hist(qherr./qh,linspace(0,quantile(qherr./qh,quant),50)); xlim([0 quantile(qherr./qh,quant)])
	namedfigure(13,'Error vs. Discharge');
%	subplot(2,2,3);
	clf();
%	datetick()
	%n = length(qh);	for idx=1:length(qh); plot(qh(idx),qherr(idx),'.','color',[1-idx/n idx/n 0]); hold on; end
	n = length(qh);
	x = (0:n-1)'/n;
	scatter(qh/1e3,qherr/1e3,2*x.^0,[x 1-x 0.*x]);
	xlim(quantile(qh/1e3,[0.01 0.99]));
	ylim(quantile(qherr/1e3,[0.01,0.99]))
	xlabel('q (1e3 m^3/s)')
	ylabel('s_q (1e3 m^3/s)');

	namedfigure(14,'Water depth and instrument depth time series')
	plot([nanmedian(bin.h)', nanmedian(bin.z)'])

	namedfigure(15,'Specific discharge 2D');
	imagesc(bin.qs);
	caxis([0 20]);

	namedfigure(16,'Rating curve vs. HADCP discharge');
	clf
	qr = interp1(rating.time,rating.discharge,hadcp.time);
	scatter(qh/1e3,qr/1e3,[],'.');
	hold on;
	plot([0 10],[0 10]);
	axis equal
	xlabel('q_{hadcp} 1e3 m^3/s');
	ylabel('q_{rating} 1e3 m^3/s');
	
	% plotting
	figure(1);
	pdfprint(['img/sanggau-q-qerr-vs-time' method{mdx} '-' num2str(worder(mdx)) '.eps']);

%	figure(2);
%	pdfprint(['img/sanggau-q-h-vs-time-' method{mdx} '-' num2str(worder(mdx)) '.eps']);

	figure(2);
	pdfprint(['img/sanggau-q_hadcp-q_rating-vs-time-' method{mdx} '-' num2str(worder(mdx)) '.eps']);

	figure(5);
	pdfprint(['img/sanggau-vel-hadcp-vadcp-' method{mdx} '-' num2str(worder(mdx)) '.eps']);

	figure(71);
	pdfprint(['img/sanggau-discharge-weights-' method{mdx} '-' num2str(worder(mdx)) '.eps']);

	figure(8);
	pdfprint(['img/sanggau-discharge-weight-coefficients-' method{mdx} '-' num2str(worder(mdx)) '.eps']);

	figure(9);
	pdfprint(['img/sanggau-q-vs-h-' method{mdx} '-' num2str(worder(mdx)) '.eps'])

	figure(13);
	pdfprint(['img/sanggau-q-vs-qerr-' method{mdx} '-' num2str(worder(mdx)) '.eps']);

	figure(16)
	pdfprint(['img/sanggau-q_rat-vs-q_hadcp-' method{mdx} '-' num2str(worder(mdx)) '.eps']);

	end

