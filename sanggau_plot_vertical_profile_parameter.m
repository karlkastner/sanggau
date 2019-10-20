% Thu Oct 23 08:44:38 CEST 2014
% Fr 31. Jul 11:04:11 CEST 2015
% plot histograms of u_*, ln_z_0 and standard errors

	% TODO compare different filter length (1,4,16,64m)
	% TODO no magic file names
	% TODO the error bars are bogus, have to be corrected for by correlation

	ps  = 2.5;
	n0  = 80;
	nf  = 0;
	ax2flag = 0;

	% number of kenel density bins
	nh    = 100;

	% quantiles of density axis limits
	qplot = 0.98;

	order = 1;

	% dependence of roughness on depth
	mode = {'direct', 'roughness'};

	if (~exist('pflag','var'))
		pflag = false;
	end
	flag = pflag;

	field_C = {'R_h','Q'}
	abc     = 'ABCDE';

	% load data
	if (~exist('reload','var') || reload)

		% load water level
		meta = sanggau_metadata();
		load([ROOTFOLDER,'/src/water-level/mat/water-level.mat']);
		level.val  = Kr(nid.sanggau_merged).depth;
		level.time = Kr(nid.sanggau_merged).time;

		% load bed form and grain size data
		gravel  = load([ROOTFOLDER,'/src/bed-material/mat/sanggau-gravel-data.mat']);
		bedform = load('mat/sanggau-bedform-data.mat');

		% load processed vadco data (joint)
		calib = load_vadcp_discharge(meta);

 		level.val = level.val  - meta.hadcp.z0;
		calib.zs0  = calib.zs0 - meta.hadcp.z0;

		% TODO no magic file names
		%f=dir([meta.ofolder, filesep, '*sanggau-dw-1-vmethod-1d-errmode-standard-discharge.mat']);
		f=dir([meta.ofolder, filesep, '*sanggau-dn-1.mat']);
	
		cs_A = calib.cs_;
	
		% load unprocessed vadcp data
		adcp = VADCP();
		for idx=1:length(meta.filename.vadcp)
			load(meta.filename.vadcp{idx},'vadcp')
			vadcp = VADCP(vadcp);
			adcp(idx) = vadcp;
		end
	
		% order processed vadcp data in time
		[t tdx] = sort(arrayfun(@(x) x.t0(1), cs_A));
		cs_A = cs_A(tdx);
		
		% order unprocessed vadcp data in time
		[t tdx]=sort(arrayfun(@(x) x.time(1),adcp));
		adcp = adcp(tdx);

		% bring into same order as calib
		t0        = calib.tstart;
		[t0_ sdx] = sort(calib.tstart);
		adcp(sdx) = adcp;
		cs_A(sdx) = cs_A;
		abc_(sdx) = abc;
		for idx=1:5;
			leg{idx} = sprintf('%s  %s',abc_(idx), datestr(cvec(t0(idx)),'dd/mm/yyyy'));
			leg__{idx} = sprintf('%s',datestr(cvec(t0(idx)),'dd/mm/yyyy'));
		end

		%%% split

		% inner region
		inner = abs(calib.N) < meta.nmax;

		% low, mid and high flow
		zs_p = linspace(min(calib.zs0),max(calib.zs0),3);
		%cvec(linspace(min(calib.h0),max(calib.h0),3));
	
		Q      = calib.Q;
		Q      = sign(Q(1)).*Q;
		R_h     = calib.radius;

		reload = 0;
	end % load data

	

	clear h
	R_h_mid = midrange(R_h);
	R_h0 = [min(R_h),R_h_mid,max(R_h)];


	N      = cs_A(1).N;
	us.A    = calib.u_s.A.val;
	ln_z0.A = calib.ln_z0.A.val;
	s_ln_z0 = NaN(size(ln_z0.A));
	s_us    = NaN(size(us.A));
	pp      = calib.perturbation.A.val;
	zb.A      = calib.zb.A;
	depth   = bsxfun(@plus,-zb.A,+(calib.zs0+meta.hadcp.z0)');
	zs      = calib.zs0;
	Q = sign(calib.Q(1))*calib.Q;
	A = calib.area;
	
	% instrument elevation above bottom
	zi      = -zb.A - meta.hadcp.z0 - meta.hadcp.redeploy.d;

	% spatial filtering of profiles, not necessary anymore, since
	% regularisation is done during discharge processing
	if (nf > 0)
		win = kaiserwin(1:nf)';
		depth  = wmeanfilt(win,depth,1);
		%us = wmeanfilt(win,us);
		ln_z0.A = wmedfilt(win,ln_z0.A);
		pp      = wmedfilt(win,pp);
	end
	
	wake.A = pp./depth;

	% cross sectional averaged quantities	
	[us.m s us.l us.u] = median_man(us.A(inner,:));
	[ln_z0.m s ln_z0.l ln_z0.u] = median_man(ln_z0.A(inner,:));
	[wake.m s wake.l wake.u] = median_man(wake.A(inner,:));

	%
	% fit profile parameters locally
	%

	% fit shear velocity
	fit.us = PowerLS();
	fit.us.fit(R_h,us.A');

	% fit roughness length
	fit.ln_z0  = PolyOLS(order);
	[void res(:,:,idx)] = fit.ln_z0.fit(R_h-R_h_mid,ln_z0.A');
	
	% fit wake parameter
	fit.wake = PolyOLS(order);
	fit.wake.fit(R_h-R_h_mid,wake.A');

	%
	% fit profile parameters cross sectionally averaged
	%

	x.R_h = R_h;
	x.Q  = Q/1e3;
	xlim_.Q = [3 13];
	xlim_.R_h = [4 14];
	xlab.R_h = 'R_h (m)';
	%xlab_.R_h = 'R_h';
	xlab_.R_h   = '(R_h-R_{h,mid})';
	xlab.Q  = 'Q (10^3 m^3/s)'
	%xlab_.Q = 'Q';
	xlab_.Q = '(Q-Q_{mid})';
	
	for idx=1:length(field_C)
		field = field_C{idx};
	
		x0.(field)           = midrange(x.(field));
	
		fit.(field).us       = PowerLS();
		fit.(field).us.fit(x.(field),cvec(us.m));
	
		fit.(field).ln_z0    = PolyOLS(1);
		fit.(field).ln_z0.fit(x.(field)-x0.(field),cvec(ln_z0.m));
	
		fit.(field).wake          = PolyOLS(1);
		fit.(field).wake.fit(x.(field)-x0.(field),cvec(wake.m));

	end
	%	
	% histograms of vertical profile parameters
	% 
	namedfigure(1,'Parameter distribution');
	clf();
	subplot(2,2,1);
	q = quantile(us.A(:),[1-qplot qplot]);
	h = linspace(q(1),q(2),nh);
	z = kernel1d(h,us.A(isfinite(us.A)));
	plot(h,z,'k');
	xlabel('u_*');
	hold on
	xlim(q);

	subplot(2,2,2);
	q = quantile(ln_z0.A(:),[1-qplot qplot]);
	h = linspace(q(1),q(2),nh);
	z = kernel1d(h,ln_z0.A(isfinite(ln_z0.A)));
	plot(h,z,'k');
	xlim(q);
	hold on
	xlabel('log10_{z0}');

	% relative error of shear velocity
	subplot(2,2,3);
	s   = s_us./us.A;
	q   = quantile(s,[1-qplot qplot]);
	h   = linspace(q(1),q(2),nh);
	z = kernel1d(h,s(isfinite(s)));
	plot(h,z,'k');
	hold on
	xlabel('\sigma/u_*')

	% relative error of roughness length
	subplot(2,2,4);
	s = s_ln_z0./ln_z0.A;
	q = quantile(s(:),[1-qplot, qplot]);
	h = linspace(q(1),q(2),nh);
	z = kernel1d(h,s(isfinite(s)));
	plot(h,z,'k');
	hold on
	xlabel('\sigma (log10_{z0})');

	for idx=1:calib.nc
		subplot(2,2,1);
		q = quantile(us.A(:),[1-qplot qplot]);
		h = linspace(q(1),q(2),nh);
		fdx = isfinite(us.A(:,idx));
		kernel1d(h,us.A(fdx,idx));
		hold on

		subplot(2,2,2);
		q = quantile(ln_z0.A(:),[1-qplot qplot]);
		h = linspace(q(1),q(2),nh);
		fdx = isfinite(ln_z0.A(:,idx));
		kernel1d(h,ln_z0.A(fdx,idx));
		hold on
		xlim(q);

		% relative error of shear velocity
		subplot(2,2,3);
		s   = s_us./us.A;
		q   = quantile(s,[1-qplot qplot]);
		h   = linspace(q(1),q(2),nh);
		fdx = isfinite(s(:,idx));
		kernel1d(h,s(fdx,idx));
		%n = size(us,1);
		%xlabel('\sigma/u_*')

		% relative error of roughness length
		subplot(2,2,4);
		s = s_ln_z0./ln_z0.A;
		q = quantile(s(:),[1-qplot, qplot]);
		h = linspace(q(1),q(2),nh);
		fdx = isfinite(s(:,idx));
		kernel1d(h,s(fdx,idx));
		%hold on
		%xlabel('\sigma (log10_{z0})');
	end % for idx
	subplot(2,2,3)
	leg_ = {'combined', arrayfun(@num2str,round(Q'/1e3),'uniformoutput',false)};
	legend(leg_{:});

	%
	% profile parameters across section
	%

	splitfigure([2 3],[1 1],flag,'Shear velocity');
	cla();
	plot(N,double(us.A),'-');
	ylabel({'$\;\;u_s$';'(m/s)$\;$'}, ...
	       'rot',0,'interpreter','latex');
	xlabel('N (m)');
	if (ax2flag)
	ax  = gca();
	ax(2) = addx(gca,'xlim',[-1 1]*max(calib.N), ...
		'xlabel','\xi','ytick',0, 'yticklabel',' ',...
		'xtick',[-0.8:0.4:0.8]*max(calib.N), ...
		'xticklabel',num2str([-0.4:0.2:0.4]'));
	linkaxes(ax,'x')
	end

	splitfigure([2 3],[1 2],flag,'Roughness length');
	cla();
	N = calib.N;
	plot(N,ln_z0.A); % f
	n = 0.5*(N(end)-N(1))*(-1:1)'/2;
	z = interp1(N,ln_z0.A,n); % f
	if (0)
	for idx=1:size(ln_z0,2)	 % f
		text(n,z(:,idx),abc_(idx)); 
		%line_fewer_markers(N,ln_z0f(:,idx),5,');
	end
	end
	hold on
	plot(N(N<n0),bedform.ln_z0*ones(size(N(N<n0))),'-k');
	plot(N(N>n0),gravel.ln_z0*ones(size(N(N>n0))),'-k');
	ylim([-10 -3.5]);
	yt = -4:0.5:1;
	set(gca,'ytick',log(10)*yt,'yticklabel',num2str(cvec(yt),'10^{%g}'));
	ylabel({'$\;z_0$';'(m)'},'rot',0,'interpreter','latex');
	xlabel('N (m)');
	xlim(limits(calib.N));
	%legend(datestr(t0));
	%C = cellfun(@(x) datestr(x,'dd/mm/yyyy'),num2cell(calib.tstart),'uniformoutput',false);
	legend('location','southwest',leg__{:},'from bed form and grain size');
	% dummy for spacing
%	ylabel('$      $','interpreter','latex','rot',0);

	if (ax2flag)
	ax    = gca;
	lim   = ylim();
	ax(2) = addx(gca,'xlim',[-0.5 0.5],'xlabel','\xi', ...
		'ytick',[], ...
		'xtick',[-0.4:0.2:0.4],'xticklabel',num2str([-0.4:0.2:0.4]'));
	linkaxes(ax,'y');
	y = [1e-4 1e-3 1e-2 1e-1]';
	%y = [1e-4 5e-4 1e-3 0.005 0.01 0.05 0.1 0.5];
	set(ax(2),'ylim',lim,'ytick',log(y),'yticklabel',num2str(log10(y),'10^{%d}'));
	box off
	if (0)
		a2 = axes('YAxisLocation', 'Right');
		set(a2, 'color', 'none');
		set(a2, 'XTick', []);
		%dlim = round(log10(exp(lim))); set(a2,'ylim',exp(lim),'ytick',[1 5 1].*10.^([dlim(1),dlim(1),dlim(2)])); ylabel(a2,'$z_0 \; (m)$','interpreter','latex','rot',0)
		dlim = round(log10(exp(lim)));
		set(a2,'ylim',lim,'ytick',log([0.005 0.01 0.05 0.1 0.5]),'yticklabel',[0.005 0.01 0.05 0.1 0.5]);
		ylabel('$z_0 \; (m)$','interpreter','latex','rot',0)
		set(a,'position',get(a2,'position'))
	end
	end

	% wake across section
	splitfigure([2 3],[1 3],flag);
	cla();
	plot(N,wake.A);
	ylim(1.5*[-1/6,0.25])
	ylabel('c/h','rot',0);
	hline([1/4, -1/6],'k--')
	xlim(limits(N));
	xlabel('N (m)');
	if (ax2flag)
	ax = gca();
	%addx(gca,'xlim',[-0.5 0.5],'xlabel','\xi','ytick',[]);
	ax(2) = addx(gca,'xlim',[-0.5 0.5],'xlabel','\xi','ytick',[], ...
		'xtick',[-0.4:0.2:0.4],'xticklabel',num2str([-0.4:0.2:0.4]'));
	linkaxes(ax,'x');
	box off
	legend(datestr(floor(calib.tstart)));
	end

	% predicted shear stress
	ax2flag = 0;
	splitfigure([2 3],[1 4],flag);
	cla();
	plot(calib.N,fit.us.predict(R_h0));
	ylabel('$u_s (m/s)$','rot',0,'interpreter','latex');
	xlabel('N (m)');
	xlim(limits(calib.N));
	legend('location','northeast',num2str(cvec(R_h0),'R_h=%2.1f m'));
	if (ax2flag)
	ax1 = gca;
	ax2 = addx(gca,'xlim',[-1 1]*max(calib.N), ...
		'xlabel','\xi','ytick',[], ...
		'xtick',[-0.8:0.4:0.8]*max(calib.N), ...
		'xticklabel',num2str([-0.4:0.2:0.4]'));
	linkaxes([ax1 ax2],'x');
	end

	% predicted roughness length
	splitfigure([2 3],[1 5],flag);
	cla();
	plot(calib.N,fit.ln_z0.predict(R_h0-R_h_mid));
	hold on
	plot(N(N<n0),bedform.ln_z0*ones(size(N(N<n0))),'-k');
	plot(N(N>n0),gravel.ln_z0*ones(size(N(N>n0))),'-k');
	%errorlines(calib.N,rpred(2).ln_z0,rpred(2).serr);
	%ylabel('$ln(z_0)$','rot',0,'interpreter','latex');
	yt = -4:0.5:1;
	set(gca,'ytick',log(10)*yt,'yticklabel',num2str(cvec(yt),'10^{%g}'));
	ylabel({'$\;z_0$';'(m)'},'rot',0,'interpreter','latex');
	xlabel('N (m)');
	xlim(limits(calib.N));
	leg___ = arrayfun(@(x) num2str(x,'R_h=%2.1f m'),R_h0,'uniformoutput',false);
	handle = legend('location','southwest',leg___{:},'from bed form and grain size');
	%handle.title = 'R_h';
	if (ax2flag)
	ax1 = gca;
	ax2 = addx(gca,'xlim',[-1 1]*max(calib.N), ...
		'xlabel','\xi','ytick',[], ...
		'xtick',[-0.8:0.4:0.8]*max(calib.N), ...
		'xticklabel',num2str([-0.4:0.2:0.4]'));
	linkaxes([ax1 ax2],'x');
	end


	% predicted wake parameter
	splitfigure([2 3],[1 6],flag);
	cla();
	plot(calib.N,fit.wake.predict(R_h0-R_h_mid));
	%errorlines(calib.N,rpred(2).ln_z0,rpred(2).serr);
	ylabel('$c/h$','rot',0,'interpreter','latex');
	xlabel('N (m)');
	xlim(limits(calib.N));
	legend('location','south',num2str(cvec(R_h0),'R_h=%2.1f m'));
	ylim(0.4*[-1 1])
	hline([1/4, -1/6],'k--')
	if (ax2flag)
	ax1 = gca;
	ax2 = addx(gca,'xlim',[-1 1]*max(calib.N), ...
		'xlabel','\xi','ytick',[], ...
		'xtick',[-0.8:0.4:0.8]*max(calib.N), ...
		'xticklabel',num2str([-0.4:0.2:0.4]'));
	linkaxes([ax1 ax2],'x');
	end

	%
	% profile parameter averaged across section
	%
	[t0_ sdx] = sort(calib.tstart);
	for jdx=1:length(field_C)

	field = field_C{jdx};
	h   = x.(field);
	h0  = x0.(field);
	h__           = linspace(0,20,100)';

	splitfigure([2 3],[2 1+3*(jdx-1)],flag,'Shear velocity');
	cla();

	cc = colormap('lines');
	for idx=1:calib.nc
		set(gca,'ColorOrderIndex',sdx(idx));
		hold on
		errorbar(h(sdx(idx)),us.m(sdx(idx)),us.m(sdx(idx))-us.l(sdx(idx)), ...
			 us.u(sdx(idx))-us.m(sdx(idx)),'ko','markerfacecolor',cc(sdx(idx),:));
	end
	[m_ ]       = fit.(field).us.predict(h__);
	plot(h__,m_,'k');

	xlim(xlim_.(field));
	xlabel(xlab.(field));
	ylabel({'$\;\;u_s$';'(m/s)$\;$'},'rot',0,'interpreter','latex');
	text(4,0.075, sprintf(['u_s = %4.2g ',field,'^{%4.2g}\nRMSE = %4.2g\n{R^{2}} = %5.3g\n'], ...
			exp(fit.(field).us.param(1)), ...
			fit.(field).us.param(2), ...
			fit.(field).us.rmse, ...
			fit.(field).us.r2));
	legend('location','southeast',leg(sdx));
	
	splitfigure([2 3],[2 2+3*(jdx-1)],flag,'Roughness length');
	cla();
	for idx=1:calib.nc
		set(gca,'ColorOrderIndex',sdx(idx));
		hold on
		%errorbar(h(sdx(idx)),ln_z0.m(sdx(idx)),ln_z0.m(sdx(idx))-ln_z0.l(sdx(idx)), ...
		%	 ln_z0.u(sdx(idx))-ln_z0.m(sdx(idx)),'ko','markerfacecolor',cc(sdx(idx),:));
		plot(h(sdx(idx)),ln_z0.m(sdx(idx)),'ko','markerfacecolor',cc(sdx(idx),:));
		hold on
		text(double(h(sdx(idx))), double(ln_z0.m(sdx(idx))), ['  ',abc_(sdx(idx))]);
	end
	[m_ s_]       = fit.(field).ln_z0.predict(h__-h0);
	plot(h__,m_,'k');
	ylabel({'$\;z_0$';'(m)'},'rot',0,'interpreter','latex');
	xlabel(xlab.(field));
	xlim(xlim_.(field));
	%ylim([-7.5 -3.5]);
	ylim(log([10^-4.001 10^-0.999]))
	if (flag)
		legend('location','southeast',leg(sdx));
	end
	yt = (-4:0);
	set(gca,'ytick',log(10)*yt,'yticklabel',num2str(cvec(yt),'10^{%g}'));
	%legend(leg{:});
	%fprintf(1,'ln_z0 = (%f + %f (h - \\bar h), RMSE = %f, {R^{2}} = %f\n',poly.param,poly.serr,poly.R2);
	text(5,-3.5,sprintf(['ln(z_0)=%1.2f%+1.2f',xlab_.(field), '\nRMSE = %0.2f\n{R^2} = %0.2f'], ...
		fit.(field).ln_z0.param, fit.(field).ln_z0.serr, fit.(field).ln_z0.R2))

	if (0)
	lim=ylim();
	ax(1) = gca;
	ax(2)=addy(ax(1));
	%a2 = axes('YAxisLocation', 'Right');
	linkaxes(ax,'x');
	set(ax(2), 'color', 'none');
	set(ax(2), 'XTick', []);
	dlim = round(log10(exp(lim)));
	%yt = logspace(-3,-1,5)';
	yt = [1e-3 3e-3 1e-2 3e-2 1e-1];
	ytl = arrayfun(@(x) sprintf('%0.3f',x),yt,'uniformoutput',false);
	set(ax(2),'ylim',lim,'ytick',log(yt),'yticklabel',vertcat(ytl));
	ylabel(ax(2),'$z_0 \; (m)$','interpreter','latex','rot',0)
	%log([0.001 0.005 0.01 0.05 0.1]), 'yticklabel', [0.001 0.005 0.01 0.05 0.1]);
	end
	
	%wake parameter
	splitfigure([2 3],[2 3+3*(jdx-1)],flag);
	cla();

	for idx=1:calib.nc
		set(gca,'ColorOrderIndex',sdx(idx));
		hold on
		errorbar(h(sdx(idx)),wake.m(sdx(idx)),wake.m(sdx(idx))-wake.l(sdx(idx)), ...
			 wake.u(sdx(idx))-wake.m(sdx(idx)),'ko','markerfacecolor',cc(sdx(idx),:));
		%text(double(h(sdx(idx))), double(m(sdx(idx))), ['  ',abc_(sdx(idx))]);
	end
	[m_ s_]       = fit.(field).wake.predict(h__-h0);
	errorlines(h__,m_,NaN*m_+s_,NaN*m_-s_,'k'); 
	text(3,0.04,sprintf(['c/h=%1.2f%+1.3f',xlab_.(field),'\nRMSE = %0.2f\n{R^2} = %0.2f'], ...
			fit.(field).wake.param, fit.(field).wake.serr, fit.(field).wake.R2));
	legend('location','southeast',leg(sdx));

	ylabel('Wake parameter');
	xlabel(xlab.(field));
	xlim(xlim_.(field));

	end % for jdx

	if (exist('pflag','var') && pflag)
		for idx=1:length(ps)
			pdfprint(11,'img/sanggau-shear-velocity-vs-n-meas',ps(idx));
			pdfprint(12,'img/sanggau-roughness-length-vs-n-meas',ps(idx));
			pdfprint(13,'img/sanggau-wake-parameter-vs-n-meas',ps(idx));

			pdfprint(14,'img/sanggau-shear-velocity-vs-n-predicted',ps(idx));
			pdfprint(15,'img/sanggau-roughness-length-vs-n-predicted',ps(idx));
			pdfprint(16,'img/sanggau-wake-parameter-vs-n-predicted',ps(idx));

			pdfprint(21,'img/sanggau-shear-velocity-vs-rh',ps(idx));
			pdfprint(22,'img/sanggau-roughness-length-vs-rh',ps(idx));
			pdfprint(23,'img/sanggau-wake-parameter-vs-rh',ps(idx));

			pdfprint(24,'img/sanggau-shear-velocity-vs-q',ps(idx));
			pdfprint(25,'img/sanggau-roughness-length-vs-q',ps(idx));
			pdfprint(26,'img/sanggau-wake-parameter-vs-q',ps(idx));
		
		end % for idx
		pflag = 0;
	end % if pflag

