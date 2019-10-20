% Mi 3. Sep 09:52:05 CEST 2014
% Karl Kastner, Berlin

function meta = sanggau_metadata()

	meta = struct();

	meta.ofolder = [ROOTFOLDER,'/src/discharge/mat/'];
	
	meta.filename.water_level = [ROOTFOLDER,'/src/water-level/mat/water-level.mat'];
	meta.filename.bscalibration = [ROOTFOLDER,'src/water-samples/mat/backscatter-calibration.mat'];

	meta.filename.obase = { ...
				[ROOTFOLDER,'src/discharge/mat/2013-12-09-sanggau_'];
		                [ROOTFOLDER,'src/discharge/mat/2014-02-20-sanggau_'];
		                [ROOTFOLDER,'src/discharge/mat/2014-04-18-sanggau_'];
		                [ROOTFOLDER,'src/discharge/mat/2014-06-18-sanggau_'];
                                [ROOTFOLDER,'src/discharge/mat/2015-04-27-sanggau_'];
	};
	meta.filename.vadcp = { ...
		[ROOTFOLDER,'dat/kapuas/2013-12-09-sanggau-transect/mat/vadcp.mat'];
		[ROOTFOLDER,'dat/kapuas/2014-02-20-sanggau-transect/mat/vadcp.mat'];
		[ROOTFOLDER,'dat/kapuas/2014-04-18-sanggau-vadcp/mat/vadcp.mat'];
		[ROOTFOLDER,'dat/kapuas/2014-06-18-sanggau-vadcp/mat/vadcp.mat'];
                [ROOTFOLDER,'dat/kapuas/2015-04-27-sanggau-vadcp/mat/vadcp.mat'];
	};

	meta.filename.hadcp  = [ROOTFOLDER,'/src/hadcp/mat/hadcp-sanggau.mat'];

	meta.filename.bank   = [ROOTFOLDER,'/dat/kapuas/gis/landsat/bankline/kapuas-bankline-2013.shp'];
	meta.filename.centre = [ROOTFOLDER,'/dat/kapuas/gis/landsat/centreline/kapuas-centreline-2013.shp'];

	meta.flag.quiver          = 1;
	meta.flag.print = false;
	meta.shear_stress_field = 'abs';

	% serial number of last water level gauge installed in Sanggau
	%meta.sn = 15796;

	%
	% bank and centre line proccessing options
	%

	% cut off radii (should be unified)
	meta.Rcut   = 1e4;
	meta.R_max  = 4000;
	meta.ds_max = 6000;
	
	% radius of curvature at the HADCP
	meta.Rc = 1500;

	% degree of interpolation polynomial
	meta.centreline.p  = 3;

	% smoothing distance for centreline
	meta.centreline.Ri = 1000;

	meta.bank.Ri   = 500;
	meta.bank.p    = 4;
	meta.bank.smoothflag = 0;
	meta.bank.iflag = true;

	meta.bmethod = 'idw';

	meta.utm.zone    = '49M';

	%
	% hadcp processing options
	%

	% sanggau vertical reference of sensor with respect to initial level
	% 6.837m - 0.15m - 2.13m = 4.557m
	meta.delta   = 4.557;

	% time and depth the HADCP was lowered during redeployment
	% with respect to the beam centre
	meta.hadcp.redeploy.t = 735650;

	meta.hadcp.redeploy.dz = -2.075;

	% position of the HADCP in utm coordinates
	% TODO this is from Anne-Kees keller location and boundaries, seems a bit too far into the section
	% TODO, check with own measurement
	meta.hadcp.x0           =   454235.618;
	meta.hadcp.y0           =    12722.186; % was incorrectly negated
	% initial installation level of the HADCP
	meta.hadcp.z0_initial   = 7.573155;
%	meta.hadcp.z0    = 4.25;
	% final installation level of the hadcp
	meta.hadcp.z0    = meta.hadcp.z0_initial+meta.hadcp.redeploy.dz;
	meta.hadcp.range = 75;
	meta.hadcp.idmax = 75;
	meta.hadcp.squeeze.dt = 1800;

	%
	% vadcp processing options
	% 

	meta.vadcp.d_transducer = 0.25;

	meta.grid_n    = Grid1();
	meta.vmethod2    = 'notime';
	meta.vdimension    = 'nz';
	% meta.vmethod   = 'regn';
%	meta.vmethod2  = 'hermite';
%	meta.roughnessmethod = 'wake';
	meta.velocity_profile = Log_profile_with_wake();
	meta.sediment_concentration_profile     = Rouse_Profile();

	meta.dt        = [];
%	meta.dz        = 0.25;
%	meta.dw        = 1;
	meta.dz        = 1;
	meta.dw        = 8; % was 16

	% tikhonov regularization option for velocity processing
	lambda         = 1e-3;
	meta.lambda.n  = 16e3*lambda;
	meta.lambda.us_scale = 10;




	% level of rim of bank at bankfull flow
	meta.cs.topofbank = 13;

	meta.cs.T_max = 600;


	% cross section for discharge computation
	switch (4)
		case {1} % visually estimated from landsat
			coord = load([ROOTFOLDER,'dat/kapuas/gis/metadata/2013_12_06_sanggau_coordinates.csv']);
			% this is the cross section based on average velocity during peak flow
			ubar = [0.0848, 1.1484]';
			h = norm(ubar);
			s = ubar(1)/h;
			c = ubar(2)/h;
			R = [c -s;s c];
			p  = coord.';
			p2 = (p(:,1) + R*(p(:,2)-p(:,1)));
			meta.cs.left  = p(:,1);
			meta.cs.right = p(:,2);
		case {2}
			% this is the cross section generated by ADCPTOOLs revision 45 for peak flow data
			% TODO dw is a quick hack to get the same discretisation for simple comparison
			meta.cs.left  = [454655.1 12261.7];
			meta.cs.right = [454204.5 12717.6];
			meta.dw = norm(meta.left-meta.right)/641;
		case {4} % TODO merge with case 1
			% this was visually estimated from bathymetry data
			% this are the cross section left and right bank points at bankful wl
			% visually estimated
			meta.cs.left     =   [454704; 12277];
			meta.cs.right    =   [454224; 12737];

	end
	meta.cs.centre   = 0.5*(meta.cs.left + meta.cs.right);
	meta.cs.dir      = normalise(meta.cs.right - meta.cs.left);

	%
	% post-processing and plotting options
	%

	% 100m
	meta.dwf   = 100;

	% inner range of cross section that is used for data processing
	meta.nmax  = 250;

	%meta.nf    = 100;
	meta.nf_   = 21;

	% bins to smooth bed level before computing wetted perimeter
	% TODO make multiple of dn
	meta.ns = 20;

% 
% times where irregularities in the data occur
% TODO merge with water level table
%

	meta.hadcp_tinvalid0 = [
	   7.355786666666666e+05  7.355787430555555e+05
	   7.356278472222223e+05  7.356279166666666e+05
           7.356391651273147e+05, 7.356494081828702e+05;    % Sanggau dry 1
	   7.356494583333333e+05  7.356494652777777e+05
           7.356535817939815e+05, 7.356651581828702e+05;    % Sanggau dry 2
	   7.357670000000000e+05  7.357676250000000e+05
           7.356707901273147e+05, 7.356761998495370e+05;    % Sanggau dry 3
	];

meta.hadcp_tinvalid = [
   7.355881250000000   7.355883402777778
   7.356201597222222   7.356253333333334
   7.356250000000000   7.356253333333334
   7.356341250000000   7.356343402777777
   7.356502500000000   7.356504375000000
   7.356511388888888   7.356516458333334
   7.356652916666667   7.356768958333333
   7.357682361111111   7.357683611111111
   7.357881041666666   7.357883333333334
   7.358305694444445   7.358307291666666
   7.357911111111111   7.357912986111112
   7.359847291666666   7.359847986111111
   7.358187361111111   7.358272361111111
   7.358374444444445   7.358374861111111
   7.356291388888889   7.356293055555555
   7.357914444444445   7.358271666666666
   7.359625208333334   7.359632361111111
   7.359845625000000   7.359846041666667
   7.359638263888889   7.359639027777778
   7.360652361111112   7.360653194444445
   7.356352986111111   7.356769027777777
   7.359650138888889   7.359667222222223
   7.360257430555555   7.360259097222222
   7.356085694444445   7.356086111111111
   7.360812777777777   inf
   ]*1e5;

	meta.hadcp_tbreak = [0,
			   7.357683611111111 % redeployment
			   7.357948680555555 % misalignment
	                   7.356994444444445 % realignment
				,inf]*1e5;

end % sanggau_metadata
