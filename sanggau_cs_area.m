% 2015-03-17 13:53:00.320229886 +0100
%

np = 6;
if (~exist('reload','var') || reload )
	meta = sanggau_metadata();
	% preload water level
	load([ROOTFOLDER,'/src/water-level/mat/water-level.mat']);
	level.val  = K(nid.sanggau_merged).slot.depth;
	level.time = K(nid.sanggau_merged).slot.time;

	% preload vadcp calibrarion discharge
	calib = load_vadcp_discharge(meta.filename.discharge,level);
	
	% convert depth with respect to instrument level
	% to cs average depth
	% TODO this should be part of load calib and water level
	calib.l0  = calib.l0  + calib.lH;
	level.val = level.val + calib.lH;
	
	reload = 0;
end

% level is whith respect to mean bottom elevation
bottom = smooth(calib.bottom,20)-calib.lH;
% area function
afunc_ = @(h) afunc(h,bottom,calib.dw);
% wetted perimeter function
pfunc_ = @(h) pfunc(h,bottom,calib.dw);

l = double(cvec(linspace(min(level.val),max(level.val))));
a = afunc_(l);
p = pfunc_(l);
ap = polyfit(l,a,np);
pp = polyfit(l,p,np);
afunc_ = @(l) polyval(ap, l);
pfunc_ = @(l) polyval(pp, l);

clf()
subplot(2,2,1)
plot(l,[a afunc_(l)]);
subplot(2,2,2)
%plot(l,[p A*pp]);
plot(l,[p pfunc_(l)]);

subplot(2,2,3)
plot(l,[a./p l]);

%hold on
%plot(l,[a p afunc_(l) pfunc_(p)]);

