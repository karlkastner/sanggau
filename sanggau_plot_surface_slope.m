
	nf = 100;
	if (~exist('reload','var')||reload)
		meta = sanggau_metadata();	
		load([ROOTFOLDER,'/src/water-level/mat/water-level.mat']);
		level.val  = K(nid.sanggau_merged).slot.depth;
		level.time = K(nid.sanggau_merged).slot.time;
		calib = load_vadcp_discharge(meta.filename.discharge,level);
		reload = 0;
	end
	g = 9.81;

	D = [];	
	U = [];
	q = [];
	u_s = [];
	N = calib.cs_(1).gridN.cX1;
	for idx=1:length(calib.cs_)
		D(:,idx) = calib.cs_(idx).gridN.val.bottom;
		U(:,idx) = calib.cs_(idx).gridN.val.U;
		u_s(:,idx) = calib.cs_(idx).gridN.val.u_s;
		q(:,idx) = calib.cs_(idx).gridN.val.q.total;
	end
	iv = (D<0);
	D(iv) = NaN;
	% compute slope
	% S = us/sqrt(gH)
	S = u_s.^2./(g.*D);

	% compute chezy coefficient
	% C = sqrt(g) u/us
	C = sqrt(g)*U./u_s;

	% weighed mean
	S_ = [];
	C_ = [];
	for idx=1:length(calib.cs_)
		fdx = abs(N)<300;
		S_(idx,1) = wmean(q(fdx,idx),S(fdx,idx));
		C_(idx,1) = wmean(q(fdx,idx),C(fdx,idx));
	end
	
	figure(2);
	clf
	subplot(2,2,1)
	plot(N,S);

	subplot(2,2,2)
	plot(N,C);



	win = hanwin(1:nf)';
	
	q   = wmedfilt(win,q,1);
	S   = wmedfilt(win,S,1);

	% filter before deviding by us!
	u_s = wmedfilt(win,u_s,1);
	U   = wmedfilt(win,U,1);
	C = sqrt(g)*U./u_s;
	
%	C = wmedfilt(win,C,1);
		
	subplot(2,2,3)
	plot(N,S);
	xlim(limits(N));
	ylim(quantile(S(:),[0.1 0.9]))

	subplot(2,2,4)
	plot(N,C);
	xlim(limits(N));
	ylim(quantile(C(:),[0.1 0.9]))

	% weighed mean
	for idx=1:length(calib.cs_)
		fdx = abs(N)<300;
		S_(idx,2) = wmean(q(fdx,idx),S(fdx,idx));
		C_(idx,2) = wmean(q(fdx,idx),C(fdx,idx));
	end
	figure(1)
	clf
	subplot(2,2,1)
	plot(calib.h0,S_,'.');
	subplot(2,2,2)
	plot(calib.h0,C_,'.');

	'relative increase of slope between low and high flow'
	% median, midrange, mean, hodges_lehman, no big difference)
	PolyOLS.slope(calib.l0,S_(:,2))./midrange(S_(:,2))*range(calib.l0)

	'relative increase of chezy coefficient between low and high flow'
	PolyOLS.slope(calib.l0,C_(:,2))./midrange(C_(:,2))*range(calib.l0)
	

