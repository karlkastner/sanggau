% Mon Jul 28 09:43:05 WIB 2014
% Karl Kastner, Berlin

% question: jackknife validation and linearity
% question: how to compute z_0? 
% question: distance from the bank
% question: ubar vs. umag (direction)

% perform calibration, estimation, and validation
% this assumes that the VADCP data was processed with respect to the same cross section

% worder : order of stage dependent specific discharge
% (0 no dependence, 1 linear)
% method of discharge weight determination
% jflag : application of jackknife error estimation
function [q qh cw bin ln_z_0_A_] = sanggau_hadcp_calibration(method,worder,jflag)

	addpath('../hadcp');

	% cross sectional velocity dimension
	hudim  = 1;
	% filter order for roughness length
	forder  = 2;
	% filter width (m)
	fwidth  = 120;

	% load metadata
	meta = sanggau_metadata();
	
	% load hadcp data
	[hadcp hlevel offset] = sanggau_load_hadcp(meta);

	% filter idepth
	% TODO make this part of the HADCP class
	fdx                 = find(hadcp.idepth_m < 0.075);
	hadcp.idepth_m(fdx) = NaN;

	load([ROOTFOLDER,'/src/water-level/mat/water-level.mat']);
	level.val  = K(nid.sanggau_merged).slot.depth;
	level.time = K(nid.sanggau_merged).slot.time;

	calib = load_vadcp_discharge(meta.filename_C,level);
	%calib = sanggau_load_vadcp(meta,level);	

	% this assumes the vadcp data for each campaigns is processed on the same grid
	Nz = calib.dis(1).z0grid.cX1();
	Nb = calib.dis(1).bgrid.cX1();
	Nq = calib.dis(1).vgrid.cX1();

	% make q0 positive according to discharge direction at first calibration campaign
	calib.q0 = calib.q0*sign(calib.q0(1));	

	hadcpX = hadcp.binX;
	hadcpY = hadcp.binY;
	
	% project the bin coordinates onto the cross section
	% TODO make this part of HDischarge
%	[Nh hadcpT] = xy2nt(hadcpX, hadcpY, calib.dis(1).xpc, calib.dis(1).ypc, calib.dis(1).cs.dir);
	[hadcp.N hadcp.T] = xy2nt(hadcp.binX, hadcp.binY, calib.dis(1).xpc, calib.dis(1).ypc, calib.dis(1).cs.dir);
%	dir    =  calib.dis(1).cs.dir;
%	xrel   =  hadcpX - calib.dis(1).xpc;
%	yrel   =  hadcpY - calib.dis(1).ypc;
%	Nh     =  dir(1)*xrel + dir(2)*yrel;
%	hadcpT = -dir(2)*xrel + dir(1)*yrel;

	% rotate the HADCP velocities onto the cross section
	% TODO use a function
	hvel   = hadcp.velocity.earth(:,:,1:2);
	nbin = size(hvel,1);
	nens = size(hvel,2);
	Rcs = [ dir(2) dir(1);
               -dir(1) dir(2)];
	s = sin(mheading_rad);
	c = cos(mheading_rad);
	% TODO, check sign convention (test this on rotating the coordinates)
	Rh  = [ c,s; -s,c];
	hvel = reshape(hvel(:,:,1:2),[],2);
	% rotate to cross section
	hvel = (inv(Rcs)*hvel')';
	hvel = reshape(hvel,nbin,nens,2);

	% TODO hvel is somehow rotated in sign
	hvel(:,:,1) = -hvel(:,:,1);


	% filter z_0
	q         = quantile(calib.ln_z_0_A(:),[0.1,0.9]);
	fdx       = find(calib.ln_z_0_A < q(1) | calib.ln_z_0_A > q(2));
	ln_z_0_A_ = calib.ln_z_0_A;
	ln_z_0_A_(fdx) = NaN;
	for idx=1:length(calib.dis)
		dn     = calib.dis(idx).z0grid.dx1;
		fn     = ceil(fwidth/dn);
		ln_z_0_A_(:,idx) = winfilt(ln_z_0_A_(:,idx),fn,forder);
	end

	% TODO quick fix
	% extrapolate roughness
	for idx=1:size(ln_z_0_A_,2)
		fdx = find(isfinite(ln_z_0_A_),1);
		ln_z_0_A_(1:fdx-1,idx) = ln_z_0_A_(fdx,idx);
	end

	% filter velocity for obstructed beams
	rmask = hadcp.rmask;
	for idx=1:2
		hvel(:,:,idx) = rmask.*hvel(:,:,idx);
	end

	[qh qherr cw bin qh_ qh_err] = hadcp_calibration( ...
		calib.t0, calib.q0, Nz, ln_z_0_A_, Nb, calib.bottom_A, Nq, calib.qs, ...
		hadcp.time, hadcp.N, ...
		hlevel, offset, hvel(:,:,hudim), hadcp.D, hadcp.pitch_rad, method, worder, jflag);

	save(['mat/sanggau-hadcp-calibration-' method '-' num2str(worder) '.mat'],'qh','qherr','cw','bin','ln_z_0_A_','hadcp','hvel','hlevel','calib','Nh');	
end % sanggau_hadcp_calibration

