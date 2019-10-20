% Mi 29. Jul 20:27:05 CEST 2015
% Karl Kastner, Berlin

meta = sanggau_metadata();

for idx=1:5
	%idx = 5;
	load(meta.filename.discharge{idx});
	X = discharge.adcp.X();
	Y = discharge.adcp.Y();
	ln_z0 = discharge.adcp.ens.ln_z0.val;
	%interp2.ln_z0 = TriScatteredInterp(X,Y,ln_z0);
	fdx = isfinite(X) & isfinite(Y) & isfinite(ln_z0);
	X = X(fdx);
	Y = Y(fdx);
	ln_z0 = ln_z0(fdx);
	
	ln_z0 = meanfilt1(ln_z0,50);
	skip = 1;
	X = X(1:skip:end);
	Y = Y(1:skip:end);
	ln_z0 = ln_z0(1:skip:end);
	
	T=delaunay(X,Y);
	figure(1);
	subplot(2,3,idx)
	patch(X(T)',Y(T)',ln_z0(T)','edgecolor','none');
	q = quantile(ln_z0,[0.16 0.84]);
	q = quantile(ln_z0,[0.05 0.95]);
	q = quantile(ln_z0,[0.02 0.98]);
	q = quantile(ln_z0,[0.001 0.999]);
	caxis(q)
	view(-45,90)
	colorbar
end
