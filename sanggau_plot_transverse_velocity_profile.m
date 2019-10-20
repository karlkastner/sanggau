% Di 2. Sep 16:25:08 CEST 2014
% 2015-01-28 17:37:40.388461825 +0100
% Karl Kastner, Berlin

if (~exist('pflag','var'))
	pflag = 0;
end
	fflag = pflag;

if (~exist('reload','var') || reload)
	abc = 'ABCDE';
	meta = sanggau_metadata();	
	% load hadcp data
	[hadcp hlevel offset] = sanggau_load_hadcp(meta);

	% must come before level
	load([ROOTFOLDER,'/src/discharge/mat/sanggau-stage-discharge.mat']);

	load([ROOTFOLDER,'/src/water-level/mat/water-level.mat']);
	level.val  = K(nid.sanggau_merged).slot.depth;
	level.time = K(nid.sanggau_merged).slot.time;

	calib = load_vadcp_discharge(meta.filename.discharge,level);
	N     = calib.N;
	width = max(calib.cs.width);
	inner = abs(N) < meta.nmax; 
	ln_z0_A = calib.ln_z0.A.val;
	ln_z0 = nanmedian(ln_z0_A(:));
	% TODO, this is already the averaged variant, so better to take the individual data
	meas.h   = bsxfun(@plus, calib.zb.mean, rvec(calib.l0));

	level.val = level.val - meta.z_offset;
	calib.l0 = calib.l0-calib.lH - meta.z_offset;	

	reload = 0;
end
	% transformation of adcp coordinate
	%N_adcp = ([meta.utm.x meta.utm.y]-meta.centre')*meta.dir;
	hadcp_range = 75;

	[t0_ sdx0] = sort(calib.t0);
	abc(sdx0) = abc;
	for idx=1:5;
		leg{idx} = sprintf('%s %s',abc(idx), datestr(cvec(calib.t0(idx)),'dd/mm/yyyy'));
		leg_{idx} = sprintf('%s',datestr(cvec(calib.t0(idx)),'dd/mm/yyyy'));
	end

	cc = colormap('lines');
	for idx=1:length(calib.cs_)
		% TODO u, q should be provided and interpolated by discharge or load function
		meas.U(1,idx) = calib.cs_(idx).U(1);
		meas.u(:,idx) = calib.cs_(idx).u1t(calib.t0(idx));
		% interpolate missing values
		meas.u(1,:)   = 0;
		meas.u(end,:) = 0;
		for jdx=1:size(meas.u,2)
			fdx = isnan(meas.u(:,jdx));
			meas.u(fdx,jdx) = interp1(N(~fdx),meas.u(~fdx,jdx),N(fdx),'linear');
		end
	end
	meas.q    = meas.u.*meas.h;
	meas.ubar = meas.U;
	meas.ubar_verify = sum(meas.q)./sum(meas.h);
	meas.qbar = mean(meas.q);
	meas.urel = bsxfun(@times,meas.u,1./meas.ubar);
	meas.qrel = bsxfun(@times,meas.q,1./mean(meas.q));
	fprintf('ubar %f verified: %f\n',[meas.ubar;meas.ubar_verify]);

	namedfigure(1,'Measured transversal profile of streamwise velocity');
	clf();
	plot(N,filtfunc(double(meas.urel),meta.nf));
	hold on
	urel_mean = mean(meas.urel,2);
	plot(N,filtfunc(urel_mean,meta.nf),'k')
	ylim([0.6 1.2]);
	legend(datestr(calib.t0,'dd/mm/yy'));
	legend('location','south',num2str(calib.h0,'%2.1f m'));
	ylim([0.8 1.2]);
	xlim([-350 350]);
	xlabel('N (m)');
	ylabel('$u/\overline{u}$','interpreter','latex');
	
	
	%figure(101);
	

	namedfigure(2,'Average measured velocity distribution and variance');
	clf();
	subplot(2,1,1)
%	meas.urelf = filtfunc(meas.urel,meta.nf);
	meas.urelf = smooth2(meas.urel,100,'lowess');
	m = mean(meas.urelf, 2);
	u_mean = m;
	s = std(meas.urelf, [], 2);
	errorlines(N,m,m-s,m+s,'k');
%	plot(N,m,'b'); hold on
%	plot(N,m*[1 1] + s*[-1 1],'b--');
	subplot(2,1,2)
	plot(N,s);
	s_ = s;
%	s_ = meanfilt1(s,50);
%	s_ = smooth(s,100,'loess');
	[n0 sd0] = fminsearch(@(n0) interp1(N,s_,n0),0);
	[n0(2,1) sd0(2,1)] = fminsearch(@(n0) -interp1(N,s_,n0),-150);
	[n0(3,1) sd0(3,1)] = fminsearch(@(n0) -interp1(N,s_,n0),150);
	[n0(4,1) sd0(4,1)] = fminsearch(@(n0) interp1(N,s_,n0),200);
	[n0(5,1) sd0(5,1)] = fminsearch(@(n0) interp1(N,s_,n0),-200);
	disp('extrema of transversal velocity profile variation');
	[n0 sd0*100]
	ndx = round(interp1(N,1:length(N),n0));
	hold on
	plot(N(ndx),s(ndx),'o');
	% inflection points
	s_ = smooth(s,200,'loess');                                             
	[n0_(1,1)] = fzero(@(n0) interp1(N,cdiff(cdiff(s_)),n0),+200);
	[n0_(2,1)] = fzero(@(n0) interp1(N,cdiff(cdiff(s_)),n0),-200);
	ndx = round(interp1(N,1:length(N),n0_));
	plot(N(ndx),s(ndx),'o');
	

	% depending on the channel curvature
	r   = meta.Rc + N;
	u   = bsxfun(@times, meas.h, 1./cvec(r));
	% normalise
	u   = bsxfun(@times, u, 1./nanmean(u));
	namedfigure(3,'Ratio depth to radius');
	plot(N,u);

	% regression of transversal water level slope
	gpoly = PolyOLS(1,false);
	gpoly.regress(N(inner),meas.u(inner,:));
	% note : blanckaert defines wududn = (R/u)*(du/dn)
	%        but here (w/u)*(du/dn) is used, because the gradient
	%        has an effect over the width, not over the radius
	wududn.val  = cvec(width).*cvec(gpoly.param(2,:)./gpoly.param(1,:));
	wududn.serr = cvec(width).*cvec(gpoly.params(2,:)./gpoly.param(1,:));
	save('mat/sanggau-wududn.mat','wududn');

	% local model of velocity distribution
	nf = 1;
	%nf = meta.nf;
	u_in     = meanfilt1(meas.urel,nf);
	q_in     = meanfilt1(meas.qrel,nf);
	upoly    = PolyOLS(1);
	qpoly    = PolyOLS(1);
	iqpoly   = PolyOLS(1);
	fnpoly0  = PolyOLS(0);
	fnpoly1  = PolyOLS(1);
	% fit
	upoly.regress(calib.h0,u_in');
	qpoly.regress(calib.h0,q_in');
	iqpoly.regress(calib.h0,bsxfun(@times,cvec(meas.qbar),1./meanfilt1(meas.q,nf)'));
	[void res0] = fnpoly0.regress(calib.h0,bsxfun(@times,meas.u',1./meas.ubar'));
	[void res1] = fnpoly1.regress(calib.h0,bsxfun(@times,meas.u',1./meas.ubar'));
%	[void res0] = fnpoly0.regress(calib.h0,bsxfun(@times,cvec(meas.ubar),1./meanfilt1(meas.u,nf)'));
%	[void res1] = fnpoly1.regress(calib.h0,bsxfun(@times,cvec(meas.ubar),1./meanfilt1(meas.u,nf)'));

	% spatial coorelation
	namedfigure(100,'Autocorrelation of the prediction residual');
	clf();
	subplot(2,2,1);
	n = sum(inner);
	acf0 = acf_man(res0(:,inner)');
%	[r0 acf0] = xar1(res0(:,inner)');
	r0 = ar1delay(acf0);
	plot(acf0,'.');
	hold on
	plot(r0(1).^(0:n-1),'k--');
%	plot(r0(2).^(0:n-1),'k--');
	title('Constant prediction');

	subplot(2,2,2);
	[r1 acf1] = xar1(res1(:,inner)');
	%[r1 acf1] = xar1(res1(:,inner)');
	plot(acf1,'.');
	hold on
	plot(r1(1).^(0:n-1),'k--');
%	plot(r1(2).^(0:n-1),'k--');
	title('Linear prediction');

	% predict
	hp = cvec(linspace(min(calib.h0),max(calib.h0),3));
	[up.val up.serr] = upoly.predict(hp);
	up.val  = up.val';
	up.serr = up.serr';
	[qp.val qp.serr] = qpoly.predict(hp);
	qp.val  = qp.val';
	qp.serr = qp.serr';
	[up.val0 up.serr0] = upoly.predict(calib.h0);
	up.val0 = up.val0';
	up.serr0 = up.serr0';

	asN.val  = cvec(upoly.param(2,:));
	asN.serr = cvec(upoly.serr);
	asN.R2   = cvec(upoly.R2);

	gdx = ~inner;
	fdx = find(inner);
	gdx = find(gdx);
	asN.val(fdx,2) = interp1([2*N(1)-N(2); N(gdx); 2*N(end-1)-N(end) ],[0; asN.val(gdx,1); 0],N(fdx),'cubic');

	u_ = 1./fnpoly1.predict(min(calib.h0) + [0 0.5 1]'*range(calib.h0));
	figure(41);
	plot(N,u_);
	xlim(limits(N));


	namedfigure(4,'Predicted velocity profile');
	clf
%	w = diff(limits(N));
	errorlines(N/width,up.val,up.serr);
	xlim(limits(N)/width);
	ylabel('$\overline u/U$','interpreter','latex','rot',0);
	xlabel('n/w');
	legend('location','south',num2str(cvec(hp),'%2.1f m'));

	% momentum distribution
	% TODO regress on its own
	namedfigure(5,'Momentum distribution');
	clf();
	%up.m = up.val.*up.h;
	up.h = bsxfun(@plus,calib.zb.mean,rvec(hp));	
	m = up.val.*up.h;
	m = bsxfun(@times,m,1./nanmean(m));
	b = bsxfun(@plus,calib.zb.mean,rvec(hp));
	b = bsxfun(@times,b,1./mean(b));

	splitfigure([2 2],[5 1],fflag);
	plot(N,up.val);
	legend(num2str(round(hp,1),'%1.1fm')); hold on;

	splitfigure([2 2],[5 2],fflag);
	plot(N,m);

	splitfigure([2 2],[5 3],fflag);
	plot(N,b);
	legend(num2str(round(hp,1),'%1.1fm')) 

	splitfigure([2 2],[5 4],fflag);
	plot(N,qp.val);

	figure(1);
	hold on
	N = cvec(N);
	u = gpoly.predict(N);
	u = bsxfun(@times,u,1./meas.ubar); %mean(meas.u));
	set(gca, 'ColorOrderIndex', 1)
	plot(N,u);

	namedfigure(6,'w/u du/dn (alpha_s)');
	clf();
	[h sdx] = sort(calib.t0);
	cc      = colormap('lines');
	for idx=1:length(calib.h0)
		errorbar(calib.h0(sdx(idx)),wududn.val(sdx(idx)),wududn.serr(sdx(idx)),'ko','markerfacecolor',cc(sdx(idx),:));
		text(double(calib.h0(sdx(idx))), double(wududn.val(sdx(idx))), ['  ',abc(sdx(idx))]);
		hold on
	end
	poly  = PolyOLS(1);
	xfunc = @(x) x; %@sqrt;
	poly.regress(xfunc(calib.h0),wududn.val);
	h_    = cvec(linspace(4,13,100));
%	min(calib.h0),max(calib.h0),100));
	[wududn_ wududn_s] = poly.predict(xfunc(h_));
	errorlines(h_,wududn_,NaN*wududn_s,[],'k');
	xlabel('z_s (m WGS84)'); %,'interpreter','latex');
	%xlabel('$H (m)$','interpreter','latex');
	ylabel('$\frac{w}{u} \frac{\partial{u}}{\partial{n}}$','interpreter','latex','rot',0);
	fprintf('w/u du/dn = (%f + %f h)\nRMSE = %f\nR^2 = %f\n',poly.param,poly.serr,poly.R2);
	legend(leg(sdx));
	str = {sprintf('$\\frac{w}{u}\\frac{du}{dn} = %1.2f%+1.2f z_s$',poly.param),
	       sprintf('RMSE$=%0.2f$',poly.serr),
	       sprintf('${R^2}=%0.2f$',poly.R2)};
	text(5,0/100,str,'Interpreter','latex');
	xlim([4 13]);
	ylim([-0.05 0.35]);
	%ylim([0 0.25]);
	
	%ldx = find(N > 0 & N<300); rdx = find(N<0&N>-300); 
	%plot(calib.h0,2*w*[mean(meas.u(ldx,:)) - mean(meas.u(rdx,:))]./mean(meas.u([ldx; rdx],:)),'r');
	
	urel_mean = mean(meas.urel,2);
	urel_mean = urel_mean/mean(urel_mean);

	% predict cs-averaged velocity
	inner = all(isfinite(meas.u),2);
	[scale u_pred] = homogenize_profile(meas.u(inner,:));

	% note: a previous version used an infinite wide cross section
	clear s_rel s_err s_dat rho res m2 fdx
	[rmse_n r2_n] = profile_prediction_error(u_pred);
	s_rel = rmse_n/s_dat;
	%[s_rel s_err s_dat rho res m2 u_pred fdx] = velocity_variation(meas.u);
	%u_pred = urel_mean*mean(meas.u) + N*(poly.predict(calib.h0)-poly.predict(mean(calib.h0)))';
	%res_ = u_pred - meas.u;

	% constant model
	fu	  = upoly.predict(midrange(calib.l0));
	pred0.u   = bsxfun(@times,cvec(fu),mean(meas.u));
	pred0.res = pred0.u - meas.u;

	pred1.u   = bsxfun(@times,up.val0,mean(meas.u));
	pred1.res = pred1.u - meas.u;
 	%bsxfun(@times,(up_.val+1),mean(meas.u_))

	namedfigure(9,'Explained variance vs. profiling range');
	clf();
%	plot(s_rel);
%	hold on;
%	plot([0 600], (1-sqrt(manningrc.lin.R2))*[1 1],'g')
%	semilogy(s_rel.^2,'k');
	semilogy(rmse_n);
	hold on;
	plot([0 600], (1-manningrc.lin.R2)*[1 1],'k--')
	ylim([1e-3 1]);
	xlim([1 length(rmse_n)]);
	xlabel('range (m)');
	ylabel('\sigma_{err}^2/\sigma_{Q}^2');

	namedfigure(10,'Velocity variation');
	clf();
	plot(N(fdx),pred0.res);
	hold on
	rmse0 = rmse(pred0.res')';
	rmse1 = rmse(pred1.res')';
	plot(N(fdx),rmse0,'k');
	plot(N(fdx),rmse1,'k--');
%	plot(N,[ones(size(N)) N]*C -1 );
	sqrt(mean(rmse1(inner).^2))/sqrt(mean(rmse0(inner).^2))
	%plot(hadcp.N,0.1);
	vline(hadcp.N(1))
	vline(hadcp.N(1)-hadcp_range); %end))
	xlabel('N (m)');
	ylabel('u_{meas} - u_{pred}')
	legend(cellfun(@(x) datestr(x,'dd/mm/yy'),num2cell(calib.t0),'uniformoutput',false),'RMSE');

	'relative error'
	(wududn.val - midrange(wududn.val))*width ./ (calib.cs.q0./calib.cs.area)
	range(wududn.val)*width
	figure(11)
	cres = cummean(res(inner));
	plot([ rmse([cummean(res(inner,:)) cummean(flipud(res(inner,:)))]')' rmse([cummean(pred1.res(inner,:)) cummean(flipud(pred1.res(inner,:)))]')'])
%	vline((N(1)+250))
	vline((N(1)+250)+120)

	% better error
	rho0  = sqrt(mean(autocorr1(pred0.res(inner,:),1)'.^2));
	rho0 = rho0(2);
	serr0 = sqrt(mean(mean(pred0.res(inner,:).^2,1)));
	rho1  = sqrt(mean(autocorr1(pred1.res(inner,:),1)'.^2));
	rho1 = rho1(2);
	serr1 = sqrt(mean(mean(pred1.res(inner,:).^2,1)));
	idx = (1:sum(inner))';
	nn = length(N);
	s0   = serr0*f_correlation(rho0,idx).*f_finite(nn,idx).*sqrt(1./(idx-1));
	s1  = serr1*f_correlation(rho1,idx).*f_finite(nn,idx).*sqrt(1./(idx-1));
	figure(12);
	clf
	%plot([s0 s1]);
	%R2 = 1-s1(2)^2/s0(2)^2
	sm = std(meas.U);
	plot([ (1-(s0./sm).^2) (1-(s1./sm).^2)] );

	namedfigure(12,'Goodness of fit');
	plot(N,[cvec(fnpoly1.R2) cvec(iqpoly.R2)]);
	legend('fn iu','fn iq')

	namedfigure(13,'Measured and predicted discharge ratio');
	clf();
	qrat = bsxfun(@times,cvec(meas.qbar),1./meanfilt1(meas.q,meta.nf)');
	plot(N,qrat);
	hold on;
	set(gca,'colororderindex',1);
	plot(N,iqpoly.predict(calib.h0),'--');


	namedfigure(14,'Measured and predicted velocity ratio');
	clf();
	urat = bsxfun(@times,cvec(meas.ubar),1./meanfilt1(meas.u,meta.nf)');
	plot(N,urat);
	hold on;
	set(gca,'colororderindex',1);
	plot(N,fnpoly1.predict(calib.h0),'--');

	namedfigure(15,'Transversal velocity profile function');
	clf();
	mu = meanfilt1(mean(meas.u,2),meta.nf);
	% note that u/U is not a distribution function
	% because of varying depth
	% TODO 1/ubar is wrong, this has to be the momentum average
	% TODO this should also be normalised before averaging (average ratio, not ratio of averages)
	mu = mu/mean(meas.ubar);
	plot(N,mu,'k');
	ylim([0.4 1.2]);
	hline(1,'--k');
	xlim([N(1) N(end)]);
	xlabel('N (m)');
	ylabel('$\frac{\bar u}{U}$ ','interpreter','latex','rot',0);
	title({'',''});
	a1 = gca;
	addx(gca,'xlim',[-0.5 0.5]*0.99,'ytick',[],'xlabel','\xi');
	axes(a1);
	ylim([0.4 1.2]);

	namedfigure(16,'Transversal velocity profile function');
	clf();
	% note that u/U is not a distribution function
	% because of varying depth
	% TODO 1/ubar is wrong, this has to be the momentum average
	% TODO use a better filter
	win = hanwin(1:round(3/2*meta.nf))';
	fn = bsxfun(@times,wmeanfilt(win,meas.u,true),1./meas.ubar);
	%fn = bsxfun(@times,meanfilt1(meas.u,meta.nf),1./meas.ubar);
	plot(N,fn(:,sdx0));
	% hadcp range
	vline(hadcp.N(1),'k--')
	vline(hadcp.N(1)-hadcp_range,'k--')
	n = 0.5*(N(end)-N(1))*([-1.1 1.1])'/2;
	fn_n = interp1(N,fn,n);
	for idx=1:size(fn_n,2)	
%		text(double(n),double(fn_n(:,sdx(idx))),abc(sdx(idx))); 
		%line_fewer_markers(N,ln_z0f(:,idx),5,');
	end

	hline(1,'--k');
	xlim([N(1) N(end)]);
	xlabel('N (m)');
	ylabel('$\frac{\bar u}{U}$ ','interpreter','latex','rot',0);
	title({'',''});
	legend('location','south',leg_(sdx));
	ylim([0.6 1.2]);
	%datestr(cvec(calib.t0),'dd/mm/yyyy'));
	a1 = gca;
	a2 = addx(gca,'xlim',[-0.5 0.5]*0.99,'ytick',[],'xlabel','\xi');
	%axes(a1);
	set(a2,'xtick',[-0.4:0.2:0.4],'XTickMode','manual')
	set(a1,'Color', 'None');
	set(a1,'position',get(a2,'position'));

	figure(102);
	S1 = rmse(smooth2(res1',100,'lowess')')';
	S0 = rmse(smooth2(res0',100,'lowess')')';
	plot(S0,S1);
	%rmse(smooth2(res1',100,'lowess')')' rmse(smooth2(res0',100,'lowess')')'])

	% R2 of the velocity distribution
	urel = bsxfun(@times,meas.u',1./meas.ubar');
	urelf = smooth2(urel',100,'lowess');
	R20 = 1 - mean(flat(S0(inner,:).^2))/var(flat(urelf(inner,:)))
	R21 = 1 - mean(flat(S1(inner,:).^2))/var(flat(urelf(inner,:)))
	srat = sqrt(mean(flat(S1(inner,:).^2))/mean(flat(S0(inner,:).^2)))
	fprintf('R2 of constant model %g\n',R20);
	fprintf('R2 of liner model %g\n',R21);
	fprintf('ratio of standard errors serr1/serr0 %g\n',srat);


	D = bsxfun(@minus,calib.zb.A,-calib.l0'); c=[]; for idx=1:size(D,1); c(idx,1) = corr((max(D(idx,:)',0)).^0.5,fn(idx,:)');  end
	fprintf('R2 of predicted velocity profile by depth profile: %f\n',nanmean(c.^2))
	%nanmedian(c.^2)

	if (exist('pflag','var') && pflag)
		pdfprint(1,'img/sanggau-Un-rel');
		pdfprint(6,'img/sanggau-w-du-dn-vs-h');
		ps = [2 2.5 3];
		for idx=1:length(ps)
			pdfprint(6,'img/sanggau-w-du-dn-vs-h',ps(idx));
			pdfprint(15,'img/sanggau-u-div-U-vs-n-joint',ps(idx))
			pdfprint(16,'img/sanggau-u-div-U-vs-n',ps(idx))
		end
%		pdfprint(100,'img/sanggau-u-pred-vs-n')
		pdfprint(9,'img/sanggau-serr-vs-hadcp-range-theoretic')
		pflag = 0;

	end

