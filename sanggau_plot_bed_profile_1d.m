% Fr 12. Sep 10:02:39 CEST 2014
% Karl Kastner, Berlin

% the result is not unique
% L-1 norm, so that the volume is accurately approximated?

	addpath([ROOTFOLDER,'/src/discharge/mat']);

	if (~exist('reload','var') || reload)
		
		% load hadcp data
		meta = sanggau_metadata();
		[hadcp hlevel offset] = sanggau_load_hadcp(meta);

		% load water level
		load([ROOTFOLDER,'/src/water-level/mat/water-level.mat']);
		level.val  = Kr(nid.sanggau_merged).depth;
		level.time = Kr(nid.sanggau_merged).time;
		
		htime  = hadcp.time;
		hdepth = hadcp.idepth_m;
		fdx    = find(htime > meta.hadcp.redeploy.t);
		hlevel = hdepth;
		hlevel(fdx) = hlevel(fdx) - meta.hadcp.redeploy.d;

		% load vadcp calibration data
		calib = load_vadcp_discharge(meta.filename.discharge,level);	
		
		reload = 0;
	end
	bottom_A = calib.zb.A;

	N = calib.cs_(1).gridN.cX1;
	N = N(:);
	b = -mean(bottom_A,2);
	fdx = find(b < b(end),1);
	N = N(fdx:end);
	b = b(fdx:end);

	% regress slope over all
	%p = [0 1/6 1/5 1/4 1/3];
	%p = [0 1/16 1/8 1/4];
	p = [0 0.1 0.2 0.3 0.4];
	A = [N.^0 N];
	clf();
	cc = colormap('lines');
	for odx=1:6
	subplot(2,3,odx)
	w = calib.cs_(1).dwidth;
	plot(w*N,b,'color',cc(1,:));
	hold on
	for idx=1:length(p)
		fdx = (N > p(idx) & N < 1-p(idx));
		%A2 = [N(fdx).^0 N(fdx)];
		%c = A2 \ b(fdx);
		%A2 = [N(fdx).^0 N(fdx) (N(fdx)-mean(N(fdx))).^2 (N(fdx)-mean(N(fdx))).^3 (N(fdx)-mean(N(fdx))).^4 (N(fdx)-mean(N(fdx))).^5];
		%A2 = [N(fdx).^0 N(fdx)];
		A2 = [N(fdx).^0 N(fdx)];
		for jdx=2:odx
			A2(:,end+1) = (N(fdx)-mean(N)).^(jdx+1);
		end
		c = A2 \ b(fdx);
		plot(w*N,A*c(1:2),'color',cc(idx+1,:));
	end % for idx
	end % for odx

