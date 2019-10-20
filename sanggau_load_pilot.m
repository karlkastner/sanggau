% Mon Dec 15 13:00:56 CET 2014
% Karl Kastner, Berlin
% TODO use SLG class (does it interpolate coordinates on hold?)

function [out sonardata] = preload_pilot(calib)
	cs = calib.cs_(1);
	% TODO no magic numbers and names
	%zone = '49M';
	load([ROOTFOLDER,'/dat/kapuas/pilot/bathymetry/15_01_2012_Sanggau/Chart_15_1_12[15].mat']);
	H         = -sonardata.depth;
	% TODO quick fix
	H(H > 15) = NaN;
	[X Y]     = latlon2utm(sonardata.lat, sonardata.long); %, zone);

	pc     = cs.centre;
	[N T]  = xy2nt(X, Y, pc(1), pc(2), cs.dir);

	switch (lower(calib.cs_(1).vmethod))
	case {'grid1','grid2'}
		% mesh
		gridN  = Grid1(0.5*cs.dwidth*[-1 1], cs.dw);
		gridN.build_index(N(:),'i1');

		% bin
		[val err] = gridN.binop('i1',@(x) mean_man(rvec(x)),H(:));
		gridN.val.zb = -val;
		gridN.err.zb = err;

		% fix missing
		grid.val.zb = fixnan(gridN.cX1,gridN.val.zb);	
		% level
		d = nanmedian(gridN.val.zb - calib.zb.median);
		gridN.val.zb = gridN.val.zb - d;

		out.N  = gridN.cX1;
		out.zb = grid.val.zb;
	case {'reg1','reg2'}
		% mesh
		gridNr = RegularizedInterpolator1();
		nn = round(cs.dwidth/cs.dw)-1;
		gridNr.remesh(0.5*cs.dwidth*[-1 1], nn);
		
		% bin
		bc.X   = 0.5*cs.dwidth*[-1 1];
		%bc.val = cs.topofbank*[1 1];
		bc.val = 0.*[1 1];
		gridNr.init(N(:),-H(:),'zb',cs.lambda.n,[],bc);

		% level
		d = nanmedian(gridNr.vali.zb - calib.zb.median);
		gridNr.vali.zb = gridNr.vali.zb - d;

		%out = gridNr;
		out.N  = gridNr.mesh.X;
		out.zb = gridNr.vali.zb;
	otherwise
		error('here');
	end
end % preload_pilot

