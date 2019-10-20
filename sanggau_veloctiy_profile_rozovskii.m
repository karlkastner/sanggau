% Sat Feb  7 14:43:23 CET 2015
% Karl Kastner, Berlin

wf = 30;
if (~exist('reload','var') || reload)
	meta = sanggau_metadata();	
	% load hadcp data
	[hadcp hlevel offset] = sanggau_load_hadcp(meta);
	load([ROOTFOLDER,'/src/water-level/mat/water-level.mat']);
	level.val  = K(nid.sanggau_merged).slot.depth;
	level.time = K(nid.sanggau_merged).slot.time;

	calib     = load_vadcp_discharge(meta.filename.discharge,level);
	reload    = 0;
end
	calib.l0  = calib.l0;
	Un        = calib.Un_A;
	Un	  = meanfilt1(Un,wf);
	Un_rel    = Un ./ repmat(nanmean(Un),length(calib.N),1);
	bottom = mean(calib.bottom.A,2);

	% acceleration by gravity [m/s]
	g = 9.81;

	Tlim = [0 0.5*pi];
	% steps along cross section (17 in rozovskii)
	nn     = 1000;
	% step width along bend
	dtheta = 0.5*pi/100;
	% rozovski uses a different van Karman constatn in bends
	kappa = 0.5;

	% Chezy's coefficient (constant)
	C = 50;
	% radius of curvature
	Rc = 1500;
	% width of the bend
	B  = 600;
	% grid along cross section
	y = B*((0:nn-1)'/(nn-1)-0.5);
	r  = Rc+y;

	% longitudinal water level slope
	% seems to have no influence on the result
%	It   = 0.0015;
	It = 1e-5;

	% local depth
	yt = 2*y/B;	

	figure(1);
	clf
	figure(2);
	clf
	cc = colormap('lines');

	% sanggau calibration level
	N = calib.N;
        H = double(calib.l0 + calib.lH);
	vavg = double(calib.cs.u0);

	figure(1); clf
	figure(10); clf
	V = NaN(nn,length(H));

	Un_sim = NaN(size(Un_rel));
	Un_scheme = NaN(size(Un_rel));
	for idx=1:length(H)

	% simulate for rectangular cross section
	hmax = H(idx);
	h = hmax.*ones(size(yt));
	v = (h/max(H)).^0.4;
	v = vavg(idx)*v/mean(v);
	fdx = h > 0;
	% solve pde
        [theta v] = rozovskii(v(fdx),h(fdx),r(fdx),It,C,kappa,g,Tlim,dtheta);
	v = real(v);
	Un_scheme(fdx,idx) = v(:,end)/mean(v(:,end));

	% repeat for real cs
	h = double(bottom+calib.l0(idx));
	v = ones(size(h));
%	v = (h/max(h)).^0.4;
	v = vavg(idx)*v/nanmean(v);
	fdx = h > 0;
	% solve pde
        [theta v] = rozovskii(v(fdx),h(fdx),r(fdx),It,C,kappa,g,Tlim,dtheta);
	v = real(v);

	Un_sim(fdx,idx) = v(:,end); %/mean(v(:,end));

	% regress slope
	reg = Theil();
	reg.regress(cvec(N),calib.Un_A(:,idx));
%	[c pserr] = medianslope(cvec(N),calib.Un_A(:,idx));
	c = reg.param;
	pserr = reg.params;
	%c = robustlinreg(cvec(N),u_A(:,idx));
	% TODO, this is not alpha_s but alpha_s*R
	as(idx,1) = c(2)/c(1);

	%[c pserr] = medianslope(cvec(N(fdx)),v(:,end));
	reg.regress(cvec(N),Un_sim(:,idx));
%v(:,idx));
	c = reg.param;
	pserr = reg.params;
	%c = robustlinreg(cvec(N),u_A(:,idx));
	% TODO, this is not alpha_s but alpha_s*R
	as_sim(idx,1) = c(2)/c(1);
	

	end
	Un_sim_rel    = Un_sim ./ repmat(nanmean(Un_sim),length(calib.N),1);

	namedfigure(1,'Observed velocity distribution');
	clf();
	plot(calib.N,Un_rel);

%	namedfigure(2,'Schematic vel dist');
%	clf();
%	plot(y,Un_scheme);

	namedfigure(3,'Simulated velocity distribution');
	clf();
	plot(calib.N,Un_sim_rel);

	namedfigure(11,'Simulated - Observed');
	plot(calib.N,Un_sim_rel-Un_rel);

	namedfigure(4,'Slope of velocity distribution');
	clf();
	bar(H,B*as)
%	ylabel('w \alpha_s R'
	ylabel('w \partial{u}/\partial{n}');
	xlabel('$\overline h$','interpreter','latex');

	figure(40);
	bar(H,[B*as B*as_sim])
	legend('observed','predicted');
	
%	fdx   = find(y > 0);
%	vsim  = nansum(V)';
%	vmeas = nansum(V(fdx,:))';
%	A = [vmeas.^0 vmeas];
%	c = A\v;
%	res = A*c - v;
%	serr2 = res'*res/sqrt(nn-2)
%	R2 = 1 - serr2/var(v)

	d = Un_sim_rel - Un_rel;
	serr2 = nanmean(d(:).^2)
	R2 = 1 - serr2./nanvar(Un_sim_rel(:))

	if (exist('pflag','var') && pflag)
%		figure(10);
%		pdfprint('img/sanggau-vel-depth-averaged-rosovsky.eps');
		figure(4);
		pdfprint('img/sanggau-slope-of-transversal-velocity-distribution.eps');
	end
	fprintf('as %f 10^-4\n',as*10^4);

