% 2017-05-15 00:05:30.081212053 +0200
function cs = sanggau_plot_backscatter_flux(pflag)

	if (nargin()<1)
		pflag = 0;
	end
	
	fflag = pflag;
	
	% TODO use URc
	nwin1 = 17; % 1
	nwin2 = 0; %151; % 101
	
	meta  = sanggau_metadata();
	calib = load_vadcp_discharge(meta)

	% load ../discharge/mat/sanggau-stage-discharge.mat
	% Qrange=[min(qm
	%qma        qmean      qmedfilt1  qmf        qmoments   qmr        
	%>> Qrange=[min(qma),max(qma)]


%	load(opt.filename.water_level);
%	level      = struct;	
%	level.time = Kr(nid.sanggau_merged).time;
%	level.val  = Kr(nid.sanggau_merged).depth;
	
	N = calib.N;
	Q = calib.Q;
	s = sign(Q(1));
	Q     = s*Q;
	q_n   = s*calib.q_tn;
	qbs_n = s*calib.sediment_discharge_tn;
	UN    = s*calib.U_tn;
	U     = s*calib.U;
	
	%bsc_n = calib.bsc_tn;
	bsc_n = qbs_n./q_n;
	H   = calib.radius;
	cs  = calib.cs_;
	
	% bc
	qbs_n(1,:)   = 0;
	qbs_n(end,:) = 0;
	qbs_n = fixnan(N,qbs_n,inf);
	bsc_n([1 end],:)=0;
	
	% intecrage
	dn  = N(2)-N(1);
	Qbs = dn*sum(qbs_n)'
	Bsc = Qbs./Q;
	
	%Qbs = nansum(qbs_n)';
	
	if (nwin1>1)
	%win   = hanwin(1:nwin1);
	%UN    = wmeanfilt(win,UN,1);
	%qbs_n = wmeanfilt(win,qbs_n,1);
	UN    = qmedfilt1(UN,nwin1);
	q_n   = qmedfilt1(q_n,nwin1);
	qbs_n = qmedfilt1(qbs_n,nwin1);
	%qbs_n = wmeanfilt(win,qbs_n,1);
	%qbs_n = wmeanfilt(win,qbs_n,1);
	end
	
	% calibration (constant coefficient)
	if (0)
	win = hanwin(1:151);
	bsc_n = wmeanfilt(win,bsc_n,1);
	sanggau_ssc_sample();
	c=cs(1).centre;
	[Ni T] = xy2nt(ssc(:,1), ssc(:,2)+1e7, c(1), c(2), cs(1).dir);
	Ni     = mean(Ni);
	bsc = interp1(N,bsc_n(:,1),Ni)
	ssc_ = mean(ssc(:,10))
	% mg/l to kg/m^3
	ssc_ = 1e-3*ssc_
	cal = ssc_/bsc
	%pause
	else
		cal = 1
	end
	
	% apply calibration (no compensation for attenuation)
	Qbs   = cal*Qbs
	qbs_n = cal*qbs_n;
	bsc_n = cal*bsc_n;
	
	'theoretic eh'
	d = [  5.1762    4.3589    3.6606    3.0627    2.5706    2.0030    1.4163    0.9153    0.6527    0.5477    0.4610    0.3884    0.3263    0.2739    0.2302    0.1953    0.1500    0.1061    0.0753    0.0532    0.0380];
	h = [ 0    0.0284    0.0227    0.0295    0.0302    0.0728    0.1314    0.2113    0.0856    0.0988    0.0667    0.0630    0.0548    0.0294    0.0162    0.0185    0.0257    0.0082    0.0039    0.0020    0.0009];
	% C = cs(1).Chezy;
	% TODO
	C = 50;
	[Qs_i, Qs] = fractional_transport_engelund_hansen(h,C,d,U(1),H(1),cs(1).width);
	%total_transport_engelund_hansen(C,sum(d.*h),U(1),H(1),cs(1).width)
	%pause
	
	% engelund-hansen (5th-power)
	if (0)
		c0       =  (U.^5) \ Qbs(:,1);
		Qbs(:,3) = c0*U.^5;
	end
	
	%c = [nanmean(Qbs)./nanmean(Q);3];
	if (0)
		[c sd] = freg(U,Qbs(:,1))
		% rms
		%f = @(Q,c) c(1)*
		jn = Jackknife(@freg);
		jn.estimate(U,Qbs(:,1));
		Qbs(:,2) = c(1)*U.^c(2);
		'parameter (linearised, jn)'
		[c jn.val]
		'parameter se'
		[sd jn.serr]
		'95% (2-sigma) confidence intervals'
	t = tinv([0.025 0.975],3)
	c + t*sd
	jn.val + t*jn.serr
end
	
	% regression along cs individual
	n1 = length(cs(1).N); %cs(1).grid_n.-1;
	cc = [];
	cc = NaN(2,n1);
	for idx=1:n1
		try
		cc(:,idx) = freg(UN(idx,:)',qbs_n(idx,:)');
		catch
		end
	end
	cc = cc';
	
	if (nwin2>1)
		win = hanwin(1:nwin2);
		cc  = wmeanfilt(win,cc,1);
	%qbs_n = wmeanfilt(win,qbs_n,1);
	end
	
	
	nlim = limits(N);
	
	
	%Q0 = [min(Q) midrange(Q) max(Q)]';
	Q0 = [1e3; 5.5e3; 1e4];
	Q0_ = linspace(0,1e4);
	x = Q-midrange(Q);
	x = x/max(x);
	x0 = [-1,0,1];
	l  = arrayfun(@(x) num2str(x,'%3.2f'),round(Q(:,1)),'uniformoutput',false);
	l0  = arrayfun(@(x) num2str(x,'%3.2f'),round(Q0(:,1)),'uniformoutput',false);
	%fit.bsc_n = PolyOLS(1);
	%fit.bsc_n.fit(x,bsc_n');
	%pred.bsc_n = fit.bsc_n.predict(x0)';
	%fit.bsc_n.fit(x,log(bsc_n'));
	%pred.bsc_n = exp(fit.bsc_n.predict(x0)');
	%fit.bsc_n.fit(x,(bsc_n'));
	%pred.bsc_n = (fit.bsc_n.predict(x0)');

if (1)	
	% bootstrap
	id_ = resample_with_replacement(1:5)';
	%id_ = nchoosek(1:5,4)'; % this is the jackknife
	%id_ = nchoosek(1:5,3)';
	%id_ = randi(5,100,1e3);
	for idx=1:size(id_,2)
		id = id_(:,idx);
		%randi(5,5,1);
		fit.qbs_n = PowerLS();
		fit.qbs_n.nonlinear = false;
		fit.qbs_n.fit(Q(id)/1e3,(qbs_n(:,id)'));
		pred.qbs_n(:,:,idx) = (fit.qbs_n.predict(Q0/1e3))';
		pred.qbs_n_(:,:,idx) = (fit.qbs_n.predict(Q0_/1e3))';

		fit.q_n = PowerLS();
		fit.q_n.nonlinear = false;
		fit.q_n.fit(Q(id)/1e3,(q_n(:,id)'));
		pred.q_n(:,:,idx) = (fit.q_n.predict(Q0/1e3))';
		pred.q_n_(:,:,idx) = (fit.q_n.predict(Q0_/1e3))';

		fit.Qbs = PowerLS();
		fit.Qbs.nonlinear = false;
		fit.Qbs.fit(Q(id)/1e3,(Qbs(id)));
		pred.Qbs(:,idx) = (fit.Qbs.predict(Q0/1e3))';
		pred.Qbs_(:,idx) = (fit.Qbs.predict(Q0_/1e3))';
	end
	pred.qbs_n_ = real(pred.qbs_n_);
	fdx = ~all(isfinite(qbs_n),2);
	pred.qbs_n_(fdx,:,:) = NaN;
	%pred.q_n = real(prep.q_n);
	pred.qbs_n  = nanmedian(pred.qbs_n,3);
	pred.q_n    = nanmedian(pred.q_n,3);
	pred.qbs_n_ = quantile(pred.qbs_n_,[0.16,0.5,0.84],3);
	%pred.q_n_ = nanmedian(pred.q_n,3);
	pred.q_n_   = quantile(pred.q_n_,[0.16,0.5,0.84],3);
	pred.Qbs_  = quantile(pred.Qbs_,[0.16,0.5,0.84],2);
	pred.Qbs  = quantile(pred.Qbs,[0.16,0.5,0.84],2);
	
	pred.q_n(isnan(pred.q_n)) = 0;
	pred.Q = dn*sum(pred.q_n);
	pred.qbs_n(isnan(pred.qbs_n)) = 0;
	%pred.Qbs = dn*nansum(pred.qbs_n);
	%pred.Qbs_ = squeeze(dn*nansum(pred.qbs_n_,1));
%	fit.qbs_n = PowerLS();
%	fit.qbs_n.fit(Q/1e3,(qbs_n'));
%	pred.qbs_n = (fit.qbs_n.predict(Q0/1e3))';
%	pred.qbs_n(isnan(pred.qbs_n)) = 0;
%	pred.Qbs = dn*sum(pred.qbs_n);
%else
%	fit.qbs_n = PolyOLS(1);
%	fit.qbs_n.fit(Q/1e3,(qbs_n'));
%	pred.qbs_n = (fit.qbs_n.predict(Q0/1e3))';
%	pred.qbs_n(isnan(pred.qbs_n)) = 0;
%	pred.Qbs = dn*sum(pred.qbs_n);
end
	
	
	pred.qbs_n_prof = bsxfun(@times,pred.qbs_n,1./pred.Qbs(:,2)');
	pred.bsc_n = pred.qbs_n./pred.q_n;
	
	%qbs_prof = bsxfun(@times,qbs_n,1./rvec(Qbs(:,1)));
	
	%prof_ = mean(qbs_n,2)/mean(Qbs(:,1));
	%fit.prof = Theil();
	%fit.prof = PolyOLS(1);
	%fit.prof.fit(x,(prof'));
	%pred.prof = fit.prof.predict(x0)';
	
	
	splitfigure([2 3],[1 1],fflag);
	cla();
	plot(N,qbs_n);
	ylabel('Suspended sand transport profile');
	xlabel('N');
	xlim(nlim);
	legend(l{:});

	splitfigure([2, 3],[4, 1],fflag);
	cla();
	qbs_profile = bsxfun(@times,qbs_n,1./nanmean(qbs_n));
	q_profile  = bsxfun(@times,q_n,1./nanmean(q_n));
	'deviation'
	median(std(qbs_profile,[],2))
	sqrt(mean(std(qbs_profile,[],2).^2))
	median(std(q_profile,[],2))
	sqrt(mean(std(q_profile,[],2).^2))
	plot(N,qbs_profile);
	ylabel('Transport proxy profile');
	xlabel('N');
	xlim(nlim);
	lh= legend(l{:}); title(lh,'Q (m^3/s)');
		
	splitfigure([2 3],[1 4],fflag);
	cla
	plot(N,bsc_n)
	ylabel('Concentration proxy')
	xlabel('N (m)');
	
	splitfigure([2 3],[1 2],fflag);
	cla();
	plot(N,pred.qbs_n)
	ylabel('Transport proxy (predicted)');
	xlabel('N (m)');

	splitfigure([2, 3],[4, 2],fflag);
	cla();
	qbs_profile = bsxfun(@times,pred.qbs_n,1./nanmean(pred.qbs_n));
	q_profile  = bsxfun(@times,pred.q_n,1./nanmean(pred.q_n));
	'deviation'
	median(std(qbs_profile,[],2))
	sqrt(mean(std(qbs_profile,[],2).^2))
	median(std(q_profile,[],2))
	sqrt(mean(std(q_profile,[],2).^2))
	plot(N,qbs_profile,'linewidth',1.5);
	ylabel('w q_s / Q_s');
	%Transport proxy profile (predicted)');
	xlabel('N (m)');
	xlim(nlim);
	lh= legend(l0{:}); title(lh,'Q (m^3/s)');
%	set(gca,'defaultaxiscolororder',[0,0,0;
 %                 0.9,0,0;
  %                0,0.7,0]);
	ylim([0,4.5]);
	
	splitfigure([2 3],[1 5],fflag);
	cla();
	plot(N,pred.bsc_n)
	ylabel('Concentration proxy (predicted)')
	xlabel('N (m)');
	
	splitfigure([2 3],[1 3],fflag);
	cla();
	w = range(N);
	plot(N,w*pred.qbs_n_prof)
	ylabel('Transport profile')
	xlabel('N (m)');
	ax    = gca();
	ax(2) =addx();
	linkaxes(ax,'x');
	xlim(ax(1),limits(N));
	set(ax(2),'xtick',-300:100:300,'xticklabel',num2str((-300:100:300)'/w));
	
	
	%splitfigure([2 2],[10 4],fflag);
	%cla();
	%bsc_n./Bsc;
	%plot(N,pred.qbs_n_prof)
	%hold on
	%plot(N,pred.prof_,'--k');
	%ylabel('Flux profile')
	%xlabel('N (m)');
	
	
	splitfigure([2 2],[2 2],fflag,'Sediment transport profile');
	cla();
	xlabel('N (m)')
	ylabel('$\bar u (m/s)$','interpreter','latex')
	xlim(nlim);
	
	%xlim(ax(2),'
	
	legend(l{:});
	
	splitfigure([2 2],[2 3],fflag,'Sediment transport profile');
	cla();
	plot(N,cc(:,2));
	%hline(nanmedian(cc(:,2)))
	hline(nanmean(cc(:,2)));
%	hline(c(2,1))
	ylabel('power');
	xlabel('N (m)')
	xlim(nlim);
	
	splitfigure([2 2],[2 4],fflag,'Sediment transport profile');
	cla();
	plot(N,cc(:,1));
	hline(nanmedian(cc(:,1)))
	hline(nanmean(cc(:,1)));
	%hline(c(1,1));
	xlabel('N (m)');
	ylabel('scale');
	ylim([0 1e-4])
	xlim(nlim);
	
	
	splitfigure([2 2],[3 1],fflag,'Sediment transport through cross section');
	cla();
	%subplot(2,2,1);
	plot(Q/1e3,Qbs/1e3,'or','markerfacecolor','r');
	xlabel('Q_w (10^3 m/s)');
	ylabel('Q_s (10^3 kg/s)');
	if (exist('c','var'))
		str=sprintf('Q_{s} = %0.3f U^{%1.1f}',c); %,'interpreter','latex')
		text(1,0.8*max(Qbs(:))/1e3,str);
	end
	%l = axis();
	xlim([0 10]);
	hold on
	plot(Q0_/1e3,pred.Qbs_/1e3,'k');
	
	splitfigure([2 2],[3 2],fflag,'Sediment transport through cross section');
	cla();
	plot(U,Qbs,'.');
	xlabel('U_w (m/s)');
	ylabel('Q_s (BS-proxy)');
	
	splitfigure([2 2],[2 3],fflag);
	cla
	plot(Q/1e3,Bsc,'o');
	ylabel('Bsc (concentration proxy)');
	xlabel('Q (10^3 m^3/s)');
	
	%cc
	if (0)
	'test'
	c = [1 5];
	n1 = 1e2;
	n2 = 1e2;
	U = linspace(0.5,1.5,n1)';
	U = repmat(U,1,n2);
	%Qbs = c(1)*U.^c(2);
	%Qbs = Qbs + 0.01*rand(size(Qbs));
	c1 = c(1);
	c2 = c(2); %+randn(100,1);
	Qbs = c1.*U.^c2 + 1*randn(n1,n2);
	
	%jn = Jackknife(@freg);
	%jn.estimate(U,Qbs(:,1));
	
	
	'lin'
	[c sd] = freg(U,Qbs);
	mean(c,2)
	sqrt(mean(sd.^2,2))
	std(c,[],2)
	
	%'jn'
	%jn.val
	%jn.serr
	
	%f = @(c) c(1)*U.^c(2) - Qbs;
	%grad(f,c)
	%A = [U.^c(2),c(1)*log(U).*U.^c(2)]
	end
	
	if (pflag>0)
		ps = 2.5;
	%	pdfprint(21,'img/sanggau-qbs-vs-qw.pdf',ps);
	%	pdfprint(11,'img/sanggau-qbs-vs-N.pdf',ps);

		pdfprint(31,'img/sanggau-Qs-vs-Qw.pdf',ps);
		pdfprint(31,'img/sanggau-Qs-vs-Qw.pdf',ps);
	end
	
end % backscatter_sanggau

function [c_ sd_] = freg(U_,Qbs_)
	for idx=1:size(U_,2)
		U = U_(:,idx);
		Qbs = Qbs_(:,idx);
	% log space as start value
	A = vander_1d(log(U),1);
	c = A \ log(Qbs);
	c(1) = exp(c(1));
	c = lsqnonlin(@(c) c(1)*U.^c(2) - Qbs,c,[0 0]);
	c_(:,idx) = c;
	if (nargout()>1)
		Qbs__ = c(1)*U.^c(2);
		ss2 = sum((Qbs__-Qbs).^2);
		% covariance matrix
		A = [U.^c(2),c(1)*log(U).*U.^c(2)];
		C = ss2/(length(Qbs)-2)*inv(A'*A);
		sd = sqrt(diag(C));
		sd_(:,idx) = sd;
	end
	end
end

