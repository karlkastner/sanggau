% Wed Jan 21 14:42:48 CET 2015
% Karl Kastner, Berlin

function compare_models(s_rel,reload)
	addpath ./hadcp

	% load calibration data

	order = 2;
	persistent calib meta
	load ../discharge/mat/sanggau-stage-discharge.mat

	if (isempty(calib) || (nargin() > 1 && reload))
	meta = sanggau_metadata();	
	% load hadcp data
	[hadcp hlevel offset] = sanggau_load_hadcp(meta);
	load([ROOTFOLDER,'/src/water-level/mat/water-level.mat']);
	level.val  = K(nid.sanggau_merged).slot.depth;
	level.time = K(nid.sanggau_merged).slot.time;
	calib = load_vadcp_discharge(meta.filename.discharge,level);
	end

	% arithmetic average
	func = @(c,Q) mean(Q);

	%
	% extract pseudo hadcp velocity (vadcp velocity over entire depth)
	%
	zi = -meta.hadcp.redeploy.d;
	% surface(calib.dis(1).vgrid.cX1,calib.dis(1).vgrid.cX2+calib.l0(1),double(calib.dis(1).vgrid.val(:,:,1))','edgecolor','none')	
	u_A = [];
	for idx=1:length(calib.l0)
		u_A(:,idx) = interp1(calib.dis(idx).gridNZ.cX2+calib.l0(idx), ...
			    double(calib.dis(idx).gridNZ.val.U)', ...
			    -meta.hadcp.redeploy.d,'linear');
	end

	%
	% filter
	%
	bottom = calib.bottom;
	qs_A   = calib.Qn_A;
	ln_z0  = median(real(calib.ln_z0_A),2);
	
	% This is the most critical part - whys should the data be filtered at all?
if(1)
	bottom = smooth(bottom,64);
	% TODO, the z0 estimation is still not satisfactory -> choose lower resolution
	ln_z0 = smooth(ln_z0,4);

	% TODO, there seems to be a filter problem at the first calibration campaign
	d=[smooth(qs_A,8)-qs_A];
	fdx=isfinite(d); l=length(isfinite(d(fdx)));
	q=quantile(d(fdx).^2,0.98);
	fdx=d.^2>q;
	qs_A(fdx) = NaN;
	qs_A = smooth(qs_A,64);

	d=[smooth(u_A,8)-u_A];
	fdx=isfinite(d); l=length(isfinite(d(fdx)));
	q=quantile(d(fdx).^2,0.98);
	fdx=d.^2>q;
	u_A(fdx)=NaN;
	u_A = smooth(u_A,64);
end
	 % plot(smooth(u0,64))
	%plot(d);
	%plot(sort(d(fdx).^2),(1:l)/l);

	% only use continuously wetted region
%	m = repmat(bottom,1,length(calib.l0)) + repmat(rvec(calib.l0),length(bottom),1);
%	m = calib.bottom_A + repmat(rvec(calib.l0),length(bottom),1);
%	m = min(m,[],2);
%	fdx = find(m>0,1,'first'):find(m>0,1,'last');
% 	m = bottom+min(calib.l0);
%	lmin = 1;
%	fdx = find(m>lmin,1,'first'):find(m>lmin,1,'last');
	bottom_A = calib.bottom_A;
	nc  = size(bottom_A,2);
	ns  = size(bottom_A,1);
	h   = min(repmat(bottom,1,nc) + repmat(rvec(calib.l0),ns,1),[],2);
	
	fdx = find( cvec((calib.dis(1).gridN.cX1 > -306) & (calib.dis(1).gridN.cX1 < 312)) & (h > 0) );
	c = [];

%	ln_z0 = median(ln_z0)*ones(size(ln_z0));

%	ivm_joint(c,zi,calib.l0,calib.cs.q0,calib.cs.area,qs_A(fdx,:),u_A(fdx,:),ln_z0(fdx),bottom(fdx),func);
%	Qp_(:,1) = hU_u_est(c,zi,calib.l0,calib.cs.q0,calib.cs.area,qs_A(fdx,:),u_A(fdx,:),ln_z0(fdx),bottom(fdx),func);
%	Qp_(:,2) = Q_u_est(c,zi,calib.l0,calib.cs.q0,calib.cs.area,qs_A(fdx,:),u_A(fdx,:),ln_z0(fdx),bottom(fdx),func);
%	Qp_(:,3) = ivm_poly(c,zi,calib.l0,calib.cs.q0,calib.cs.area,qs_A(fdx,:),u_A(fdx,:),ln_z0(fdx),bottom(fdx),func,order);
%	Qp_(:,4) = vpm_poly(c,zi,calib.l0,calib.cs.q0,calib.cs.area,qs_A(fdx,:),u_A(fdx,:),ln_z0(fdx),bottom(fdx),@wavg,order);
%	figure(1)
%	plot(calib.cs.q0,Qp_);

	% limited range
	% flip, because the ADCP is installed at the end
	fdx = fliplr(fdx);
	err = NaN(length(fdx),1);
	func = @wavg;
	mode = {'linear','affine','h'}
	for mdx=1:3
	for idx=1:length(fdx)
		h0      = calib.l0;
		q0      = calib.cs.q0;
		u0_adcp = u_A(fdx(1:idx),:)';
		bottom_  = bottom(fdx);
		bottomi = bottom(fdx(1:idx));
		ln_z0i  = ln_z0(fdx(1:idx),:);
		ln_z0_  = ln_z0(fdx);

		ivm = IVM(mode{mdx});
		ivm.fit(h0, u0_adcp, q0, zi, bottom_, bottomi, ln_z0_, ln_z0i);
		R2(idx,1,mdx) = ivm.R2; 
		serr(idx,1,mdx) = ivm.serr0;
		rat(idx,1,mdx)  = ivm.rat;
	
		vpm = VPM(mode{mdx});
		vpm.fit(h0, u0_adcp, q0, zi, bottom_, bottomi, ln_z0_, ln_z0i);
		R2(idx,2,mdx) = vpm.R2; 
		serr(idx,2,mdx) = vpm.serr0;
		rat(idx,2,mdx)  = vpm.rat;

		sdm = SDM(mode{mdx});
		sdm.fit(h0, u0_adcp, q0, zi, bottom_, bottomi, ln_z0_, ln_z0i);
		R2(idx,3,mdx) = sdm.R2; 
		serr(idx,3,mdx) = sdm.serr0;
		rat(idx,3,mdx)  = sdm.rat;

		esm = ESM(mode{mdx});
		esm.fit(h0, u0_adcp, q0, zi, bottom_, bottomi, ln_z0_, ln_z0i);
		R2(idx,4,mdx) = esm.R2; 
		serr(idx,4,mdx) = esm.serr0;
		rat(idx,4,mdx)  = esm.rat;
	end % for idx

	figure(40+mdx);
	clf();
	X = calib.N(fdx)-calib.N(fdx(1))+1;
	semilogy(X,1-R2(:,:,mdx));
	legend('location','southeast','ivm','vpm','sdm','esm','rating curve');
	ylabel('s_{err}^2/s_{Q_0}^2');
	xlabel('profiling range (m)');
	xlim([X(1) X(end)])
%	ylim(10.^[-3.5   -0.5]);
	hold on; plot([0 600], (1-powerrc.R2)*[1 1],'g')
	pdfprint(['img/sanggau-serr-vs-hadcp-range-' mode{mdx}]);
	end % for mdx

	figure(100);
	clf();
	% run sp discharge first
%	semilogy(X,[1-R2(:,1,1) s_rel(1:sum(fdx)).^2 (1-powerrc.R2)*X.^0]);
	n = min(length(s_rel),length(X));
	semilogy(X(1:n),[1-R2(1:n,1,1)  s_rel(1:n)'.^2 (1-powerrc.R2)*X(1:n)'.^0]);

	%hold on; plot([0 600], (1-powerrc.R2)*[1 1],'g')
	ylim([10^-2.5 1]);
	legend('theoretic','IVM','rating curve');
	ylabel('s_{err}^2/s_{Q_0}^2');
	xlabel('profiling range (m)');
	xlim([X(1) X(end)])
	pdfprint(['img/sanggau-serr-vs-hadcp-range']);

	figure(5);
	plot(rat)
	legend('location','southeast','ivm','vpm','sdm','esm');
	pdfprint('img/sanggau-affine-contribution-vs-hadcp-range');

	figure(2)
	clf();
	subplot(3,1,1)
	plot(err)
	subplot(3,1,2)
	plot([mean(w,2) mean(ww,2)])
	subplot(3,1,3)
	plot(s2)
	figure(3)
	clf
	plot(w.*u_A(fdx,:) - ones(length(fdx),1)*calib.cs.q0'); hold on;
end % compare models

% ivm individual
function [Qp Q_i res resn] = ivm_individual(c,zi,l0,Q0,area,q0,u0,ln_z0,bottom,func)
	nc = length(Q0);
	ns = size(u0,1);
	for idx=1:ns
		A = [area area.*u0(idx,:)'];
		% estimate calibration parameter
		c = A \ Q0;
		% predict
		Q_i(idx,:) = (A*c)';
	end
	Qp = func(Q_i./u0,u0,zeros(ns,1),ones(ns,1));
	Qp = cvec(Qp);
	res = Qp - Q0;
	resn = norm(res);
end


% VPM individual
function [Qp Q_i res resn] = vpm_individual(c,zi,l0,Q0,area,q0,u0,ln_z0,bottom,func)
	nc = length(Q0);
	ns = size(u0,1);
	h = repmat(rvec(l0),ns,1) + repmat(bottom,1,nc);
	%u = u0.*( (log(h)-1-repmat(ln_z0,1,nc))./(log(zi)-repmat(ln_z0,1,nc)) );
	z = bottom + zi;
	u = u0.*( (log(h)-1-repmat(ln_z0,1,nc))./(repmat(log(z),1,nc)-repmat(ln_z0,1,nc)) );
	for idx=1:ns
		A = [area area.*u(idx,:)'];
		% estimate calibration parameter
		c = A \ Q0;
		% predict
		Q_i(idx,:) = (A*c)';
	end
	Qp = func(Q_i./u,u,zeros(ns,1),ones(ns,1));
	Qp = cvec(Qp);
	res  = Qp - Q0;
	resn = norm(res);
end

% VPM individual
function [Qp Q_i res resn] = sdm_individual(c,zi,l0,Q0,area,q0,u0,ln_z0,bottom,func)
	nc = length(Q0);
	ns = size(u0,1);
	h = repmat(rvec(l0),ns,1) + repmat(bottom,1,nc);
	z = bottom + zi;
	q = h.*u0.*( (log(h)-1-repmat(ln_z0,1,nc))./(repmat(log(z),1,nc)-repmat(ln_z0,1,nc)) );
	%q = h.*u0.*( (log(h)-1-repmat(ln_z0,1,nc))./(log(zi)-repmat(ln_z0,1,nc)) );
	for idx=1:ns
		A = [ones(nc,1) q(idx,:)'];
		% estimate calibration parameter
		c = A \ Q0;
		% predict
		Q_i(idx,:) = (A*c)';
	end
	Qp = func(Q_i./q,q,zeros(ns,1),ones(ns,1));
	Qp = cvec(Qp);
	res  = Qp - Q0;
	resn = norm(res);
end


% Q_q fix
function [Qp_ Qp] = hU_u_est(c,zi,l0,Q0,area,q0,u0,ln_z0,bottom,func)
	% fixed part
	Q_q = mean(repmat(rvec(Q0),size(u0,1),1)./q0,2);
	z = bottom+zi;
	for idx=1:length(l0)
		% C ignored (cancels if assumed constant)
		C = ones(size(u0,1),1);
		% S ignored (cancels if assumed constant)
		S = ones(size(u0,1),1);
		% variable part
		h = (bottom + l0(idx));
		Ubar = C.*(S.^0.5).*(h.^0.5);
		u = Ubar.*(log(z)-ln_z0)./(log(h)-1-ln_z0);
		hU_u = h.*Ubar./u;
		% predict
		Qp(:,idx) = Q_q.*hU_u.*u;
		% average
		Qp_(idx,1) = func(Q_q.*hU_u, Qp(:,idx));
	end
end	

% 1 fix (all modelled)
function [Qp_ Qp] = Q_u_est(c,zi,l0,Q0,area,q0,u0,ln_z0,bottom,func)
	z = bottom + zi;
	% for each calibration
	for idx=1:length(l0)
		% S ignored (cancels if assumed constant)
		S = ones(size(u0,1),1);
		% C ignored (cancels if assumed constant)
		C = ones(size(u0,1),1);
		% depth
		h = (bottom + l0(idx));
		% model Qeren = int C S
		Q = sum( C.*(S.^0.5).*(h.^1.5));
		% model u	
		Ubar = C.*(S.^0.5).*(h.^0.5);
		u = Ubar.*(log(z)-ln_z0)./(log(h)-1-ln_z0);
		% predict 
		Qp(:,idx) = (Q./u).*u0(:,idx);
		% average
		Qp_(idx,1) = func((Q./u), Qp(:,idx));
	end
end

% index velocity method with rating curve
function [Qp_ Qp] = ivm_power(c,zi,l0,Q0,area,q0,u0,ln_z0,bottom,func)
	u0_ = cvec(mean(u0));
	l0 = double(l0);
	Q0 = double(Q0);
	f = @(c) (c(1)*abs(l0 - c(2)).^c(3)).*u0_.*sign(l0 - c(2));
	c = [600 min(l0)-0.1 1];
	opt = optimset();
	% opt.TolX = 1e-12;
	opt.MaxFunEvals = 1000;
	opt.MaxIter = 1000;
	c = lsqnonlin(@(c) Q0 - f(c),double(c),[],[],opt);
	Qp_ = f(c);
end

% index velocity method with polynomial rating curve
function [Qp_ Qp] = ivm_poly(c,zi,l0,Q0,area,q0,u0,ln_z0,bottom,func,order)
	u0_ = cvec(mean(u0));
	l0 = double(l0);
	Q0 = double(Q0);
%	f = @(c) (c(1)*abs(l0 - c(2)).^c(3)).*u0_.*sign(l0 - c(2));
%	c = [600 min(l0)-0.1 1];
%	opt = optimset();
	A = vander_1d(l0,order);
	c = A \ (Q0./u0_);
	% opt.TolX = 1e-12;
	%opt.MaxFunEvals = 1000;
	%opt.MaxIter = 1000;
	%c = lsqnonlin(@(c) Q0 - f(c),double(c),[],[],opt);
	Qp_ = (A*c).*u0_;
end

% vpm (individual for every bin)
function [Qp_ Qp s2w w ww res] = vpm_poly(c,zi,l0,Q0,area,q0,u0,ln_z0,bottom,func,np)
	n  = size(q0,1);
	l0  = double(l0);
	Q0  = double(Q0);
	A   = vander_1d(l0,np);
	Q0_ = repmat(rvec(Q0),size(u0,1),1);
	c   = (A \ (Q0_./u0)');
	w   = (A*c)';
	Qp  = w.*u0;
	% standard error of regression
	res = Qp - Q0_;
	s2w  = sqrt(sum(res.*res,2)/(size(res,2)-np-1));
	s2w = median(s2w)*ones(size(s2w));
	s2u = 0.01*ones(n,1);
	[Qp_ ww] = func(w,u0,s2w,s2u);
	ww = ww*size(ww,1);
end

