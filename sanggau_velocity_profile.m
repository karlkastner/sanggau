% Thu 10 Aug 15:56:10 CEST 2017

	lw = 1.5;
	ps = 2.5;
	name  = 'sanggau';
	if (~exist('reload','var') || reload)
		opt = sanggau_metadata();
		% load water level
		load([ROOTFOLDER,'/src/water-level/mat/water-level.mat']);
		level.val  = Kr(nid.sanggau_merged).depth;
		level.time = Kr(nid.sanggau_merged).time;

		% load processed vadcp data
		calib = load_vadcp_discharge(opt);

 		level.val       = level.val - opt.hadcp.z0_initial;
		calib.zs0       = calib.zs0 - opt.hadcp.z0_initial;
		calib.zb.median = calib.zb.median - opt.hadcp.z0_initial;
		calib.zb.A      = calib.zb.A - opt.hadcp.z0_initial;

		% load hadcp data
	        [hadcp hlevel_ hoffset_] = sanggau_load_hadcp(opt);                       

		n_hadcp = [hadcp.N(1) hadcp.N(opt.hadcp.idmax)];
		reload  = 0;
	end
	if (~exist('pflag','var'))
		pflag = 0;
	end
	fflag = pflag;

	%nc_    = 4;
	Q     = calib.Q;
	neff  = sum(Q).^2/sum(Q.^2);

	Rh    = calib.perimeter;
	x     = calib.zs0; %Rh-midrange(Rh);
	xmid  = midrange(x);
	x0    = [midrange(x) min(x) max(x)];
	%zb = calib.zb.median;
	zb = calib.zb.A;
%	ylim_.acf = [-0.8 1]*1e-4;
%	ylim_.acf = [[-1e-4 3e-4]];
%	ylim_.acf = 1e-5*[-1.75 9.0];
	ylim_.acf = 1e-3*[-0.3 1.6];


	nf.pre  = [200 0];
	nf.post = 200;

	%vmode = 'direct';
	vmode  = 'profile';
	acfbias = true;
	vorder  = 1;
	torder  = 1;
	
	plot_velocity_profile();

