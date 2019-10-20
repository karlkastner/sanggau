% Thu  9 Feb 14:45:22 CET 2017
% Karl Kastner, Berlin

%
% TODO, put regression of boat vel into CrossSection
%
	meta = sanggau_metadata();
	iname = meta.filename.vadcp;
	obase = meta.filename.obase;

	if (~exist('reload','var') || reload)	

		load([ROOTFOLDER,'/src/water-level/mat/water-level.mat']);
		level.val  = K(nid.sanggau_merged).slot.depth;
		level.time = K(nid.sanggau_merged).slot.time;
		calib = load_vadcp_discharge(meta.filename.discharge,level);

	end

	figure(1)
	clf
	for idx=1:length(iname)
	load(iname{idx})
	adcp=VADCP(vadcp)

	obj=calib.cs_(1);
	centre = obj.centre;
	centre(2) = centre(2)-1e7;
	adcp.xy2nts(1,centre,obj.dir,obj.dwidth,obj.T_max);
	plot(adcp.ens.N,adcp.ens.T)
	hold on

	obj.gridN.build_index(adcp.ens.N,'i1'); 
	f=@(time,t,n,s,v) deal(nanmean(v),NaN);
	[v(:,idx) s] = obj.gridN.binop('i1',f, adcp.time, adcp.ens.T, adcp.ens.N, adcp.ens.N, ...
			hypot(adcp.btvel.ship(:,1),adcp.btvel.ship(:,2)));
	end
	figure(2)
	plot(v)
	nanmedian(v)
	nanmedian(v(:))

	
	

