% Tue 20 Dec 14:52:46 CET 2016
% Karl Kastner, Berlin

function [out] = sanggau_load_bed_level_2016(recompute,cs)
	meta = sanggau_metadata();

	ofile = 'mat/sanggau-bed-level-2016.mat';

	if (nargin() < 1 || (~isempty(recompute) && ~recompute))
		load(ofile);
	else

	% TODO no magic file names
	load([ROOTFOLDER,'/src/bathymetry/mat/sonar-2016.mat']);
	load([ROOTFOLDER,'/src/bathymetry/mat/gps-2016.mat']);
	addpath([ROOTFOLDER,'/src/sanggau/']);

	% TODO no magic numbers
	dw    = 1;
	d_max = 350;
	y0 = 1e7;
	% nf    = 25;

	% cross section centre	
	pc = [0.5*(meta.cs.xpmax+meta.cs.xpmin);
	      0.5*(meta.cs.ypmax+meta.cs.ypmin)];

	% cross section direction
	dir = [(meta.cs.xpmax-meta.cs.xpmin);
	       (meta.cs.ypmax-meta.cs.ypmin)];

	% cross section width
	dwidth = norm(dir);
	dir   = dir/dwidth;

	X  = cvec(double(slg.utm.X));
	Y  = cvec(double(slg.utm.Y))-y0;
	D  = cvec(slg.Depth());
	t  = cvec(slg.time());

	% TODO interpolate missing x, but not d
	fdx = isfinite(X) & isfinite(D);
	X   = X(fdx);
	Y   = Y(fdx);
	D   = D(fdx);
	t   = t(fdx);

	% crop to sanggau cross section
	d = hypot(X - pc(1), Y - pc(2));
	fdx = d < d_max;

	X = X(fdx);
	Y = Y(fdx);
	D = D(fdx);
	t = t(fdx);
	tlim = limits(t);

	simple = 0;
	if (simple)
		load ../water-level/mat/kapuas-surface-level-2016.mat
		gridN  = Grid1(0.5*dwidth*[-1 1], dw);
		gridN.val.z_s = interp1(ws.S,ws.alt,0)*ones(gridN.nX1-1,1);
	else
		

	% fetch GPS coordinates near the cross section
	gt = gps.time;
	gX = gps.utm.X;
	gY = gps.utm.Y-y0;
	gZ = gps.alt;

	d  = hypot(gX - pc(1), gY - pc(2));
	fdx = d < d_max;
	gt = gt(fdx);
	gX = gX(fdx);
	gY = gY(fdx);
	gZ = gZ(fdx);

	% interpolate to same time span and sampling interval
	tlim(1) = min(t(1),gt(1));
	tlim(2) = max(t(end),gt(end));
	t_ = tlim(1):1/86400:tlim(2);
	X_ = interp1(t,X,t_,'linear')-pc(1);
	Y_ = interp1(t,Y,t_,'linear')-pc(2);
	gX_ = interp1(gt,gX,t_,'linear')-pc(1);
	gY_ = interp1(gt,gY,t_,'linear')-pc(2);
	X_(isnan(X_))   = 0;
	Y_(isnan(Y_))   = 0;
	gX_(isnan(gX_)) = 0;
	gY_(isnan(gY_)) = 0;

	% TODO antenna alt + transducer depth
	% assign altitude

	% the start time of the echo sounder is set manually and only an approximate reference
	% determine the exact time lag by cross correlation
	xc      = cvec(xcorr(gX_,X_,'unbiased'));
	xc(:,2) = cvec(xcorr(gY_,Y_,'unbiased'));
	tc = ((1:length(xc))' - length(xc)/2)/86400;
	% discard junk at series start
	xc(1:500,:) = 0;
	% best fit lag in taps
	[mv mdx] = max(xc);
	% best fit lag in time [days]
	dt = mean(tc(mdx))

	% interpolate gps values to echo sounder samples
	% approximate distance from antenna to transducer
	z_offset = 0.8; 
	t0 = median(gt);
	t0(2) = median(t);

	c.X      = interp1(gt,gX,t+dt,'linear');
	c.Y      = interp1(gt,gY,t+dt,'linear');
	c.z_s     = interp1(gt,gZ,t+dt,'linear')-z_offset;

	if (0)
		% test
		mdx = mean(mdx);
		gX_ = circshift(cvec(gX_),[-mdx,0]);
		gY_ = circshift(cvec(gY_),[-mdx,0]);
		xc = cvec(xcorr(gX_,X_,'unbiased'));
		xc(:,2) = cvec(xcorr(gY_,Y_,'unbiased'));
		xc(1:500,:) = 0;
		[mv mdx] = max(xc);
		mdx
		length(xc)/2
		mdx - length(xc)/2
		tc(mdx)
	
		figure(1)
		clf
		plot(X,Y);
		hold on
		plot(gX,gY)
	
		figure(2)
		clf
		plot(t,[X c.X]-pc(1))
	
		figure(3)
		clf
		plot(t,[Y c.Y]-pc(2))
	
		figure(4)
		clf();
		plot(xc)
		
		figure(5);
		clf
		plot(X,Y);
		hold on
		plot(gX,gY);
	end % if test


	[N T]  = xy2nt(c.X, c.Y, pc(1), pc(2), dir);

	switch (lower(cs.vmethod))
	case {'grid1','grid2'}
		gridN  = Grid1(0.5*dwidth*[-1 1], dw);

		gridN.build_index(N(:),'i1');

		% surface elevation
		[val err] = gridN.binop('i1',@(x) mean_man(rvec(x)),c.z_s);
		gridN.val.z_s  = val;
		gridN.err.z_s  = err;

		% depth
		[val err] = gridN.binop('i1',@(x) mean_man(rvec(x)),D);
	
		% fix gaps
		dn_max = 5;
		val = fixnan(gridN.cX1,val,dn_max);
		%fdx      = isnan(val);
		%cX1	 = gridN.cX1;
		%val(fdx) = interp1(cX1(~fdx),val(~fdx),cX1(fdx),'linear','extrap');
	
		gridN.val.depth = val;
		gridN.err.depth = err;
	
		% bed elevation
		gridN.val.z_b = gridN.val.z_s - gridN.val.depth;

	case {'reg1','reg2'}
		% mesh
		gridNr = RegularizedInterpolator1();
		nn     = round(cs.dwidth/cs.dw)-1;
		gridNr.remesh(0.5*cs.dwidth*[-1 1], nn);
		
		% bin
		bc.X   = 0.5*cs.dwidth*[-1 1];
		%bc.val = cs.topofbank*[1 1];
		bc.val = 0.*[1 1];
		gridNr.init(N(:),D(:),'depth',cs.lambda.n,[],bc);
		gridNr.init(N(:),c.z_s,'zs',cs.lambda.n,[],[]);

		% level
		% d = nanmedian(gridNr.vali.zb - calib.zb.median);
		gridNr.vali.zb = gridNr.vali.zs - gridNr.vali.depth;

		%out = gridNr;
		out.N     = gridNr.mesh.X;
		out.zb    = gridNr.vali.zb;
		out.depth = gridNr.vali.depth;
		out.zs    = gridNr.vali.zs;
	otherwise
		error('here');
	end
	
	out.t0 = t0;
	out.tstart = t0(1);

	end % else of if simple

		save(ofile,'out');
	end
end % sanggau_load_bed_level_2016

