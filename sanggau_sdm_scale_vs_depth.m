% Fr 19. Jun 09:44:10 CEST 2015
% Karl Kastner, Berlin

% TODO uncertainty bands
% TODO separate contribution by momentum redistribution and roughness
% TODO rename calib.lH

nplot  = 100;
N0     = -250:100:250;
direct = false;

nf     = 100;
win    = hanwin(1:nf)';

nf_    = 40;
win_   = hanwin(1:nf_);

if (~exist('pflag','var'))
	pflag = 0;
end

% load calibration level
if (~exist('reload','var') || reload)
	meta  = sanggau_metadata();
	load([ROOTFOLDER,'/src/water-level/mat/water-level.mat']);
	level.val  = K(nid.sanggau_merged).slot.depth;
	level.time = K(nid.sanggau_merged).slot.time;

	calib = load_vadcp_discharge(meta.filename.discharge,level);

	% load transversal and vertical profile parameter
	load('mat/sanggau-calibration-error.mat');
	vert  = vertical.roughness;
	trans = transversal;

	% subtract measurement error from systematic error
	% TODO, this is a quick fix
	for pdx=1:2
		trans(pdx).poly.serr = sqrt(max(0,trans(pdx).poly.serr.^2 - rmse(trans(1).camp.res').^2));
%		vert(pdx).poly.serr  = sqrt(max(0,vert(pdx).poly.serr.^2 - rmse(vert(1).camp.res').^2));
		vert(pdx).rpoly.serr  = sqrt(max(0,rmse(vert(pdx).predict.res_r,2).^2 - rmse(vert(1).camp.res_r,2).^2))';
	end

	% corelation of residual in fz and fr
	rho_nz = rho_nz;

	% cross section grid coordinates
	N  = calib.N;

	% depth across section during individual calibrations
	meas.h   = bsxfun(@plus, calib.zb.mean, rvec(calib.l0));

	% level above mean channel depth (initial instrument level)
	lH0 = calib.lH;

	% surface elevation above mean bottom for diagram
	zs = cvec(linspace(min(calib.h0),max(calib.h0),nplot));

	% midrange of calibration surface level
	zm = midrange(calib.h0);

	% instrument elevation level after redeployment
	za = calib.lH - meta.hadcp.redeploy.d;

	% bottom elevation level
	zb	 = -calib.zb.mean;

	cc = colormap('lines');
	for idx=1:calib.nc
		meas.u(:,idx) = calib.cs_(idx).u1t(0);
		% interpolate missing values
		meas.u(1,:)   = 0;
		meas.u(end,:) = 0;
		for jdx=1:size(meas.u,2)
			fdx = isnan(meas.u(:,jdx));
			meas.u(fdx,jdx) = interp1(N(~fdx),meas.u(~fdx,jdx),N(fdx),'linear');
		end
		% extract velocity at HADCP depth
if (0)
		interp     = TriScatteredInterp(flat(calib.cXX1{idx}),double(calib.h0(idx))+flat(calib.cXX2{idx}),flat(double(calib.U{idx})));
		%meas.uadcp(:,idx) = interp(N,double(meas.h(:,idx)*(calib.h0(idx)-lH)./calib.h0(idx))); %/meas.h(mdx,idx)));
		meas.uadcp(:,idx) = interp(N,double(za)*ones(size(N))); %/meas.h(mdx,idx)));
end
	end

	% specific discharge
	meas.q      = meas.u.*meas.h;
	%
	inner	    = abs(N) < meta.nmax;

	meas.ubar   = sum(meas.q)./sum(meas.h);

	meas.u      = wmeanfilt(win,meas.u);
if (0)
	meas.uadcp  = wmeanfilt(win,meas.uadcp);
end

	% scale of velocity at instrument depth to depth average velocity
	meas.fz     = vert(1).camp.ln_z0;

	% scale of depth average to cross section average velocity
	meas.fn     = trans(1).f0;

	% filter for plot
	meas.fz     = wmeanfilt(win,meas.fz);
	meas.fn     = wmeanfilt(win,meas.fn);

if (0)
	% combined scale
	meas.ff     = meas.fz.*meas.fn;
end

	reload      = 0;
end

%
% predict velocity scale at calibration
%
p0 = [];
p  = []
for pdx=1:2
 % at calibration
 if (direct)
	[pc(pdx).fz.val pc(pdx).fz.serr pc(pdx).fz.sp] = vert(pdx).fpoly.predict(calib.h0);
	 pc(pdx).fz = structfun(@transpose,pc(pdx).fz,'uniformoutput',false);
 else
	 [pc(pdx).ln_z0.val pc(pdx).ln_z0.serr pc(pdx).ln_z0.sp] = vert(pdx).rpoly.predict(calib.h0);
	  pc(pdx).ln_z0 = structfun(@transpose,pc(pdx).ln_z0,'uniformoutput',false);
	 [pc(pdx).fz.val df pc(pdx).fz.serr pc(pdx).fz.sp void pc(pdx).fz.res] = ...
 		log_profile(rvec(calib.h0),zb,za,pc(pdx).ln_z0.val, pc(pdx).ln_z0.serr, pc(pdx).ln_z0.sp, []); %ln_z0_resm');
	 %[pc(pdx).fz.val df pc(pdx).fz.serr pc(pdx).fz.sp void pc(pdx).fz.res] = ...
 	%	vertical_profile(rvec(calib.h0),zb,za,pc(pdx).ln_z0.val, pc(pdx).ln_z0.serr, pc(pdx).ln_z0.sp, []); %ln_z0_resm');
	 [pc(pdx).fz.val df pc(pdx).fz.serr pc(pdx).fz.sp void pc(pdx).fz.res] = ...
 		log_profile(rvec(calib.h0),zb,za,pc(pdx).ln_z0.val, pc(pdx).ln_z0.serr, pc(pdx).ln_z0.sp, []); %ln_z0_resm');
 end
 [pc(pdx).fn.val pc(pdx).fn.serr pc(pdx).fn.sp]       = trans(pdx).poly.predict(calib.h0);
 pc(pdx).fn = structfun(@transpose,pc(pdx).fn,'uniformoutput',false);

 n = length(inner);
 m = 1;
% g1              = ar1_var_factor(trans(pdx).model.rho,n,m);
 g1              = 1;
 pc(pdx).fn.serr = g1*pc(pdx).fn.serr;
 pc(pdx).fn.sp   = g1*pc(pdx).fn.sp;


 % combine the scale
 pc(pdx).ff	                   = pc(pdx).fz.val.*pc(pdx).fn.val;
 pc(pdx).s                         = sqrt( (pc(pdx).fz.serr./pc(pdx).fz.val).^2 ... 
					   + 2*rho_nz(pdx)*(pc(pdx).fz.serr./pc(pdx).fz.val).*(pc(pdx).fn.serr./pc(pdx).fn.val) ...
				           + (pc(pdx).fn.serr./pc(pdx).fn.val).^2 );
 pc(pdx).sp                        = sqrt( (pc(pdx).fz.sp./pc(pdx).fz.val).^2 ...
					  + 2*rho_nz(pdx)*(pc(pdx).fz.sp./pc(pdx).fz.val).*(pc(pdx).fn.sp./pc(pdx).fn.val) ...
					  + (pc(pdx).fn.sp./pc(pdx).fn.val).^2 );
% [fz_rm acfm]                      = ar1(pc(pdx).fz.res);
% fz_rho(pdx) = vert(pdx).model.rho;

 % predict velocity scale over entire water level range

 if (direct)
	[p(pdx).fz.val p(pdx).fz.serr p(pdx).fz.sp] = vert(pdx).fpoly.predict(zs);
         p(pdx).fz = structfun(@transpose,p(pdx).fz,'uniformoutput',false);
 else
	 [p(pdx).ln_z0.val p(pdx).ln_z0.serr p(pdx).ln_z0.sp] = vert(pdx).rpoly.predict(zs);
	 p(pdx).ln_z0 = structfun(@transpose,p(pdx).ln_z0,'uniformoutput',false);
	 %[p(pdx).fz.val df p(pdx).fz.serr p(pdx).fz.sp void]     = vertical_profile(rvec(zs),zb,za,p(pdx).ln_z0.val, p(pdx).ln_z0.serr, p(pdx).ln_z0.sp);
	 [p(pdx).fz.val df p(pdx).fz.serr p(pdx).fz.sp void]     = log_profile(rvec(zs),zb,za,p(pdx).ln_z0.val, p(pdx).ln_z0.serr, p(pdx).ln_z0.sp);
 end
 [p(pdx).fn.val p(pdx).fn.serr p(pdx).fn.sp]          = trans(pdx).poly.predict(zs);
 p(pdx).fn = structfun(@transpose,p(pdx).fn,'uniformoutput',false);

 p(pdx).fn.serr = g1*p(pdx).fn.serr;
 p(pdx).fn.sp   = g1*p(pdx).fn.sp;

 % filter for better plotting
 p(pdx).fn.val  = wmeanfilt(win_,(p(pdx).fn.val);
 p(pdx).fn.serr = wmeanfilt(win_,(p(pdx).fn.serr);
 p(pdx).fn.sp   = wmeanfilt(win_,(p(pdx).fn.sp);
 p(pdx).fz.val  = wmeanfilt(win_,(p(pdx).fz.val);
 p(pdx).fz.serr = wmeanfilt(win_,(p(pdx).fz.serr);
 p(pdx).fz.sp   = wmeanfilt(win_,(p(pdx).fz.sp);


 % correlation coefficients
 fn_r(pdx) = trans(pdx).model.rho;
 fz_r(pdx) = vert(pdx).model.rho;


 % combine the scale

 p(pdx).ff	                     = bsxfun(@times,p(pdx).fz.val,p(pdx).fn.val);
 p(pdx).s                            = sqrt( (p(pdx).fn.serr./p(pdx).fn.val).^2  ...
					     + 2*rho_nz(pdx)*(p(pdx).fn.serr./p(pdx).fn.val).*(p(pdx).fz.serr./p(pdx).fz.val) ...
		                             + (p(pdx).fz.serr./p(pdx).fz.val).^2 );
 p(pdx).sp                           = sqrt( (p(pdx).fn.sp./p(pdx).fn.val).^2 ...
					     + 2*rho_nz(pdx)*(p(pdx).fn.sp./p(pdx).fn.val).*(p(pdx).fz.sp./p(pdx).fz.val) ...
                                               + (p(pdx).fz.sp./p(pdx).fz.val).^2 );



end % for pdx

% predict the scale parameter for selected points in the cross section
for idx=1:length(N0)

	% index of point in the cross section
	[mv mdx] = min(abs(N-N0(idx)));

	if (0)
		fprintf(1,'N %f za %f zb %f\n',calib.N(mdx),za,zb(mdx));
		fprintf(1,'midrange ln_z0 %f Q/q %f \n',p(1).ln_z0(mdx),p(1).fn.val(mdx));
	end

% linear fit
if (0)
	A = vander_1d(zs,1);
	c = A\f(:,2);
	f(:,5) = A*c;
end

	scale = 1./interp1(zs,p(1).ff(mdx,:),zm);

	namedfigure(1,'HADCP scale parameter vs. stage');
	subplot(2,3,idx); cla();
	% predicted scale
	errorlines(zs,scale*[p(1).ff(mdx,:)' p(2).ff(mdx,:)'], scale*[NaN(100,1), p(2).s(mdx,:)']);
	% measured scale (ratio)
	hold on;
%	plot(calib.h0,meas.ff(mdx,:)'/midrange(meas.ff(mdx,:)),'ko','markerfacecolor','k');
	xlabel('H (m)');
	ylabel('$\frac{f_q f_u}{f_q f_u |_0}$','interpreter','latex','rot',0);
	legend('location','southeast','constant','linear','f_q-linear','f_u-linear','appr. linear')
	ylim([0.75 1.2]);
	xlim(limits(calib.h0))
	title(sprintf('n = %dm',round(N0(idx))));

	namedfigure(2,'Parameter estimates');
	subplot(2,3,idx);
	cla();
	plot(zs,scale*[p(2).fn.val(mdx,:)' p(2).fz.val(mdx,:)'])
%	plot(zs,scale*[p(2).fn.val(mdx,:)' p(2).fz.val(mdx,:)' p(2).fz.bias(mdx,:)'])

%	namedfigure(10+idx,'Standard error vs stange');
	namedfigure(3,'Standard error vs. stage');
	subplot(2,3,idx); cla
	plot(zs,[(p(2).fn.serr(mdx,:)./p(2).fn.val(mdx,:))' (p(2).fz.serr(mdx,:)./p(2).fz.val(mdx,:))']);
	leg = legend('$\frac{\hat \sigma_{f_t}}{f_t}$', '$\frac{\hat \sigma_{f_v}}{f_v}$');
	set(leg,'interpreter','latex');
	leg.Color = 'none';

end % for idx (selected points in the cross section

set(0,'defaultAxesColorOrder',[0 0 1; 1 0 0; 0 0.6 0]);

namedfigure(30,'Standard error vs. stage');
clf();
% mean reseponse variable
if (pflag) figure(31); clf(); else subplot(2,2,1); end
subplot(2,2,1)
plot_style(zs,[  sqrt(mean((p(1).fn.serr(inner,:)./p(1).fn.val(inner,:)).^2))' ...
	  ,sqrt(mean((p(1).fz.serr(inner,:)./p(1).fz.val(inner,:)).^2))' ...
          ,sqrt(-2*rho_nz(1)*mean((p(1).fn.serr(inner,:)./p(1).fn.val(inner,:).*(p(1).fz.serr(inner,:)./p(1).fz.val(inner,:)))))' ] );
leg = legend('$\frac{\hat \sigma_{f_t}}{f_t}$', '$\frac{\hat \sigma_{f_v}}{f_v}$' ...
             ,'$\sqrt{- 2 \rho_{t,v} \frac{\hat \sigma_{f_v}}{f_v}\frac{\hat \sigma_{f_t}}{f_t}}$');
set(leg,'interpreter','latex');
	leg.Color = 'none';
 %ylim([0 0.05]);
xlim(limits(zs));
xlabel('H (m)');
ylabel('%','rot',0);
percenttick('y');

% prediction response constant parameter
if (pflag) figure(32); clf(); else subplot(2,2,2); end
plot_style(zs,[  sqrt(mean((p(1).fn.sp(inner,:)./p(1).fn.val(inner,:)).^2))' ...
	  ,sqrt(mean((p(1).fz.sp(inner,:)./p(1).fz.val(inner,:)).^2))' ...
          ,sqrt(-2*rho_nz(1)*mean((p(1).fn.sp(inner,:)./p(1).fn.val(inner,:).*(p(1).fz.sp(inner,:)./p(1).fz.val(inner,:)))))' ] );
leg = legend('$\frac{\hat \sigma_{f_t}}{f_t}$', '$\frac{\hat \sigma_{f_v}}{f_v}$' ...
             ,'$\sqrt{- 2 \rho_{t,v} \frac{\hat \sigma_{f_v}}{f_v}\frac{\hat \sigma_{f_t}}{f_t}}$');
set(leg,'interpreter','latex');
	leg.Color = 'none';
%set(leg,'FontSize',1.5*get(leg,'Fontsize'));
 %ylim([0 0.05])
xlim(limits(zs));
xlabel('H (m)');
ylabel('%','rot',0);
ylim([0 6]/100+1e-4*[-1 1]);
percenttick('y');

% mean response variable
if (pflag) figure(33); clf(); else subplot(2,2,3); end
plot_style(zs,[  sqrt(mean((p(2).fn.serr(inner,:)./p(2).fn.val(inner,:)).^2))' ...
	  ,sqrt(mean((p(2).fz.serr(inner,:)./p(2).fz.val(inner,:)).^2))' ...
          ,sqrt(-2*rho_nz(2)*mean((p(2).fn.serr(inner,:)./p(2).fn.val(inner,:).*(p(2).fz.serr(inner,:)./p(2).fz.val(inner,:)))))' ] );
hold on
set(gca,'colororderindex',1);
%plot(zs,[ sqtudt(mean((p(1).fn.serr(inner,:)./p(1).fn.val(inner,:)).^2))', ...
%	  sqrt(mean((p(1).fz.serr(inner,:)./p(1).fz.val(inner,:)).^2))'],'--');
leg = legend( '$\frac{\hat \sigma_{f_t}}{f_t}$', '$\frac{\hat \sigma_{f_v}}{f_v}$' ...
             ,'$\sqrt{- 2 \rho_{v,t} \frac{\hat \sigma_{f_v}}{f_v}\frac{\hat \sigma_{f_t}}{f_t}}$');
set(leg,'interpreter','latex');
	leg.Color = 'none';
%set(leg,'FontSize',1.5*get(leg,'Fontsize'));
ylim([0 0.05])
xlim(limits(zs));
ylabel('%','rot',0);
ylim([0 6]/100+1e-4*[-1 1]);
xlabel('H (m)');
percenttick('y');

% prediction error variable
if (pflag) figure(34); clf(); else subplot(2,2,4); end
plot_style(zs,[  sqrt(mean((p(2).fn.sp(inner,:)./p(2).fn.val(inner,:)).^2))' ...
	  ,sqrt(mean((p(2).fz.sp(inner,:)./p(2).fz.val(inner,:)).^2))' ...
          ,sqrt(-2*rho_nz(2)*mean((p(2).fn.sp(inner,:)./p(2).fn.val(inner,:).*(p(2).fz.sp(inner,:)./p(2).fz.val(inner,:)))))' ] );
hold on
set(gca,'colororderindex',1);
%plot(zs,[ sqrt(mean((p(1).fn.sp(inner,:)./p(1).fn.val(inner,:)).^2))', ...
%	  sqrt(mean((p(1).fz.sp(inner,:)./p(1).fz.val(inner,:)).^2))'],'--');
leg = legend('$\frac{\hat \sigma_{f_t}}{f_t}$', '$\frac{\hat \sigma_{f_v}}{f_v}$' ...
             ,'$\sqrt{- 2 \rho_{v,t} \frac{\hat \sigma_{f_v}}{f_v}\frac{\hat \sigma_{f_t}}{f_t}}$');
set(leg,'interpreter','latex');
	leg.Color = 'none';
%set(leg,'FontSize',1.5*get(leg,'Fontsize'));
ylim([0 0.05]);
xlim(limits(zs));
xlabel('H (m)');
ylabel('%','rot',0);
ylim([0 6]/100+1e-4*[-1 1]);
percenttick('y');

	namedfigure(4,'Predicted and measured scale parameter during calibration');
	clf();
	subplot(2,2,1);
	plot(N,pc(2).fz.val);
	hold on;
	set(gca,'colororderindex',1);
	plot(N,meas.fz,'--');
	title('f_v');
	ylim([0.5 1.5]);

	subplot(2,2,2);
	plot(N,pc(2).fn.val);
	hold on;
	set(gca,'colororderindex',1);
	plot(N,meas.fn,'--');
	title('f_t')
	ylim([0.5 1.5]);

	subplot(2,2,3);
	plot(N,pc(2).ff);
	hold on;
	set(gca,'colororderindex',1);
%	plot(N,meas.ff,'--');
	title('f_t f_v');
	ylim([0.5 1.5]);

	subplot(2,2,4)
	plot(N,wmeanfilt(win,calib.ln_z0.A.val);
	hold on;
	set(gca,'colororderindex',1);
%	switch (vdx)
%	case {3,4}
	plot(N,pc(2).ln_z0.val,'--');
%	end

	namedfigure(5,'');
	clf();
	if (pflag) figure(51); clf(); else subplot(2,2,1); end
	contourf(N,zs,p(2).s',0.01*(0:10));
	ylabel('H (m)');
	xlabel('N (m)');
	title('Variable parameter mean response error');
	if (pflag) figure(52); clf(); else subplot(2,2,2); end
	contourf(N,zs,p(1).s',0.01*(0:10));
	ylabel('H (m)');
	xlabel('N (m)');
	title('Constant parameter mean response error');

	if (pflag) figure(53); clf(); else subplot(2,2,3);
		hold on;
		title('Variable parameter prediction error');
	end
	contourf(N,zs,p(2).sp',0.01*(0:10));
	ylabel('H (m)');
	xlabel('N (m)');
	h = percenttick('c');
	ylabel(h,'%','rot',0);

	if (pflag) figure(54); clf(); else
		subplot(2,2,4);
		title('Constant parameter prediction error');
		hold on
	end
	contourf(N,zs,p(1).sp',0.01*(0:10));
%	imagesc(N,zs,p(1).sp');
	ylabel('H (m)');
	xlabel('N (m)');
	h = percenttick('c');
	ylabel(h,'%','rot',0);

	namedfigure(8,'Contribution of fz and fn to the error')
	subplot(2,2,1)
	contourf( (p(2).fn.serr./p(2).fn.val)' ,0.01*(0:10));
	ylabel('H (m)'); xlabel('N (m)');
	subplot(2,2,2)
	contourf( (p(2).fz.serr./p(2).fz.val)' ,0.01*(0:10));
	ylabel('H (m)'); xlabel('N (m)');
	subplot(2,2,3)
	contourf( (p(1).fn.sp./p(1).fn.val)',0.01*(0:10));
	ylabel('H (m)'); xlabel('N (m)');
	subplot(2,2,4)
	contourf( (p(1).fz.sp./p(1).fz.val)',0.01*(0:10));
	ylabel('H (m)'); xlabel('N (m)');

	% error reduction by spatial averaging at mid flow
	p(1).fn.sp_  = p(1).fn.sp./p(1).fn.val;
	p(1).fz.sp_  = p(1).fz.sp./p(1).fz.val;
	p(1).fz.sp_E = sqrt(mean(p(1).fz.sp_(inner,nplot/2).^2));
	p(1).fn.sp_E = sqrt(mean(p(1).fn.sp_(inner,nplot/2).^2));

	p(2).fn.sp_  = p(2).fn.sp./p(2).fn.val;
	p(2).fz.sp_  = p(2).fz.sp./p(2).fz.val;
	p(2).fz.sp_E = sqrt(mean(p(2).fz.sp_(inner,nplot/2).^2));
	p(2).fn.sp_E = sqrt(mean(p(2).fn.sp_(inner,nplot/2).^2));
	n = sum(inner);
	m = (1:n)';
	for pdx=1:2
		% TODO function name : this is the var, not the factor!
		s2_fn(:,pdx) = trans(pdx).model.s2*ar1_var_factor(trans(pdx).model.rho,n,m);
		s2_fz(:,pdx) = vert(pdx).model.s2*ar1_var_factor(vert(pdx).model.rho,1e12,m);
		s2_fzfn(:,pdx) = abs(2*rho_nz(pdx)*sqrt(s2_fn(:,pdx)).*sqrt(s2_fz(:,pdx)));
	end
	n=m;
%	s2_fn(:,1) = f_finite(660,n).^2.*(sumpower2(fn_r(1),n)./n)*p(1).fn.sp_E.^2;
%	s2_fz(:,1) =                 (1./n).*(sumpower2(fz_r(1),n)./n)*p(1).fz.sp_E.^2;
%	s2_fn(:,2) = f_finite(660,n).^2.*(sumpower2(fn_r(2),n)./n)*p(2).fn.sp_E.^2;
%	s2_fz(:,2) =                 (1./n).*(sumpower2(fz_r(2),n)./n)*p(2).fz.sp_E.^2;
	figure(9);
	subplot(2,2,1)
%	plot(acfm);
	subplot(2,2,2)
%	plot(acf);

	namedfigure(10,'HADCP error estimate with increasing range');
	clf();
	if (pflag) figure(101); clf(); else subplot(2,2,1); end
	plot_style(n,sqrt([s2_fn(:,1) s2_fz(:,1) s2_fzfn(:,1)]));
	ylim([0 0.05]);
	xlim(limits(n));
	xlabel('Range (m)');
	ylabel('%');
	title({'',''});
	ylim([0 6]/100+1e-4*[-1 1]);
	percenttick('y');
	%leg = legend('$\frac{\hat \sigma_{f_t}}{f_t}$', '$\frac{\hat \sigma_{f_v}}{f_v}$');
	leg = legend( '$\frac{g_t}{m} \frac{\hat \sigma_{f_t}}{f_t}$' ...
                     ,'$\frac{g_v}{m} \frac{\hat \sigma_{f_v}}{f_v}$' ...
		     ,'$\sqrt{-2 \rho_{v,t} \frac{g_{t,v}}{m} \frac{\hat \sigma_{f_v}}{f_v} \frac{\hat \sigma_{f_t}}{f_t} }$' ...
			);
	set(leg,'interpreter','latex');
	leg.Color = 'none';
%	set(leg,'FontSize',1.5*get(leg,'Fontsize'));
%	addx(gca,'xlabel','Range (%)','xlim',100*[0 n(end)]/660,'ytick',[]);

	if (pflag) figure(102); clf(); else subplot(2,2,2); end
	plot_style(n,sqrt([s2_fn(:,2) s2_fz(:,2) s2_fzfn(:,2)])); hold on
%	set(gca,'colororderindex',1);
%	plot(n,sqrt([s2_fn(:,1) s2_fz(:,1) s2_fzfn(:,1)]),'--');
	ylim([0 0.05]);
	xlim(limits(n));
	xlabel('Range (m)');
	ylim([0 6]/100+1e-4*[-1 1]);
	ylabel('%');
	percenttick('y');
	title({'',''});
	leg = legend( '$\frac{g_t}{m} \frac{\hat \sigma_{f_t}}{f_t}$' ...
                     ,'$\frac{g_v \hat g_v}{m} \frac{\hat \sigma_{f_v}}{f_v}$' ...
		     ,'$\sqrt{-2 \rho_{v,t} \frac{g_{v,t}}{m} \frac{\hat \sigma_{f_v}}{f_v} \frac{\hat \sigma_{f_t}}{f_t} }$' ...
		     ... %,'$-2 \rho_{v,t} \sqrt{ \frac{g_v}{m}\frac{g_t}{m} \frac{\hat \sigma_{f_v}}{f_v} \frac{\hat \sigma_{f_t}}{f_t} }$' ...
		);
	set(leg,'interpreter','latex');
	leg.Color = 'none';
%	set(leg,'FontSize',1.5*get(leg,'Fontsize'));
%	addx(gca,'xlabel','Range (%)','xlim',100*[0 n(end)]/660,'ytick',[]);

	namedfigure(1000,'rmse over the entire water level range for single point velocity measurements along the cs');
	subplot(2,2,1);
	plot(N,[rmse(p(1).sp'); rmse(p(2).sp')]'); xlim(limits(N))
	lh =legend('linear','constant');
	title('ovr hydrograph')
	ylim([0 0.3]);

	subplot(2,2,2);
	plot(N,[rmse(pc(1).sp'); rmse(pc(2).sp')]'); xlim(limits(N))
	title('during calibration')
	ylim([0 0.3]);

	subplot(2,2,3);
	plot(N,[ rmse(transversal(1).model.res')', ...
                 rmse(vert(1).model.res_f')', ...
                 rmse(vert(1).model.res_flin')', ...
 	         rmse(pc(1).fz.sp')', ...
		 rmse(pc(1).fn.sp')' ]);
% 	         rmse(pc(1).fz.serr')', ...
%		 rmse(pc(1).fn.serr')' ]);
	ylim([0 0.3]);
	legend('transversal','vertical','vertical linearised')

	figure(1001);
	clf
%	subplot(2,2,2);
	y = wmeanfilt(win,[rmse(pc(1).sp'); rmse(pc(2).sp')]');
	plot(N,y);

	abc = 'ABCDE';
	n = 0.5*(N(end)-N(1))*(0.9*[-1 1])'/2;
	z = interp1(N,y,n);
	for idx=1:size(y,2)	
		text(n,z(:,idx),abc(idx)); 
		%line_fewer_markers(N,ln_z0f(:,idx),5,');
	end

	xlim(limits(N))
	ylim([0 0.1]);
	percenttick('y')
	xlabel('N (m)');
	ylabel('$\frac{\varepsilon_{\hat Q}}{{\hat Q}}$','interpreter','latex','rot',0);
	legend('location','north','A constant model','B linear model');
	a1 = gca;
	a2 = addx(gca,'xlim',[-0.5 0.5]*0.99,'ytick',[],'xlabel','t');
	set(a2,'xtick',[-0.4:0.2:0.4],'XTickMode','manual')
	set(a1,'Color', 'None');
	set(a1,'position',get(a2,'position'));

if (exist('pflag','var') && pflag)
%	pdfprint('img/sanggau-hadcp-scale-parameter-vs-stage');
	yl = [0 0];
	for idx=[32 34 101 102]
		figure(idx);
		yl_ = ylim();
		yl(1) = min(yl(1),yl_(1));
		yl(2) = max(yl(2),yl_(2));
	end
	for idx=[32 34 101 102]
		figure(idx);
		ylim(yl);
	end

	ps = [2 2.5];
	for idx=1:length(ps)
		pdfprint(32,'img/hadcp-prediction-error-1d-constant-parameter',ps(idx));
		pdfprint(34,'img/hadcp-prediction-error-1d-variable-parameter',ps(idx));
% 2d plots cost time
%		pdfprint(53,'img/hadcp-prediction-error-2d-variable-parameter',ps(idx));
%		pdfprint(54,'img/hadcp-prediction-error-2d-constant-parameter',ps(idx));
		pdfprint(101,'img/hadcp-prediction-error-vs-range-constant-parameter',ps(idx));
		pdfprint(102,'img/hadcp-prediction-error-vs-range-variable-parameter',ps(idx));
		pdfprint(1001,'img/hadcp-prediction-error-vs-n',ps(idx));
	end

	pflag = 0;
end % if pflag

