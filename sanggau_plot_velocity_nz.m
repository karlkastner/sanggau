% Fri 14 Sep 16:11:13 CEST 2018

function sanggau_plot_velocity_nz(id)

	meta = sanggau_metadata();

	iname_C = meta.filename.vadcp;
	base = iname_C{id}(1:end-4);
	name = [base, CrossSection.optstr(meta), '.mat'];
	load(name);
	% /home/pia/phd/dat/kapuas/2013-12-09-sanggau-transect/mat/vadcp-grid-notime-dn-8-dz-1.mat

	arrowcol = 0.0*[1,0,0];
	vscale = 2*4;
	wscale = 2*0.5;
	scale = [20,40];
	view_ = [0,90];
	cmap = parula(17);
	
	idx = 1;
	sN  = 1;

	N     = sN*cs(idx).N;
	U_tnz = sN*cs(idx).U_tnz();
	V_tnz = sN*cs(idx).V_tnz();
	W_tnz = cs(idx).W_tnz();
	ulim(idx,:) = [-1,1]*quantile(abs(flat(U_tnz)),[0.99]);
	zb = cs(idx).var_n('zb');
		ylim_(idx,:) = [min(zb),0];

	nlim(idx,:) = limits(cs(idx).N);

	namedfigure(1+100*idx,'Streamwise velocity');
	clf();
	%cs(idx).grid_nz.mesh.plot(mean(U_tnz,2));
	cs(idx).plot_nz(N,U_tnz,1);
	hold on
	plot(N,zb,'-k');
	switch(meta.flag.quiver)
	case {1}
		cs(idx).plot_nz_quiver(V_tnz,W_tnz,1,scale,sN);
	case {2}
		% TODO, make this a function in cross section
		% get coordinates
		dn = 30;
		n0 = innerspace(-cs(idx).transect.dwidth/2,cs(idx).transect.dwidth/2,round(cs(idx).transect.dwidth/dn));
		% resample h
		h0 = -interp1(N,zb,n0,'linear');
		nn = [];
		ss = [];
		zz = [];
		dz = 2;
		for ndx=1:length(n0)
			%nz = round(h0/dz);	
			%s0 = innerspace(0,1,nz);
			z0 = -(dz/2:dz:h0(ndx)-dz/4)';
			nz = length(z0);
			s0 = -z0/h0(ndx);
			nn = [nn;ones(nz,1)*n0(ndx)];
			ss = [ss;s0];
			zz = [zz;z0];
		end
		ss = double(ss);
		% resample v and w (same as daspect(2)/ daspect (1))
		%ws = 2*range(ylim_(idx,:))/range(nlim(idx,:));
		aspect = range(nlim(idx,:))/range(ylim_(idx,:));
		vs = vscale*range(nlim(idx,:))/dn;
		ws = wscale*range(ylim_(idx,:))/dz;
		v  = interp2(N,cs(idx).grid_nz.S0,V_tnz,nn,ss);
		w  = interp2(N,cs(idx).grid_nz.S0,W_tnz,nn,ss);
		dt = 1;
		if (1)
		plot([nn nn+vs*v*dt,NaN*nn]',[zz,zz+ws*w*dt,NaN*zz]','-','linewidth',1.5,'color',arrowcol);
		% origin marker
		plot([nn],[zz],'.','markersize',8,'color',arrowcol);
		%plot([nn],[zz],'*k');
		%plot([nn],[zz],'*k','markersize',16);
		else
			quiver(nn,zz,dt*v,dt*w,0)
		end
	end
	caxis([ulim])
	view(view_);
	colormap(cmap)
	xlabel('N (m)');
	ylabel('z (m)');
	xlim(nlim(idx,:));
	ylim(ylim_(idx,:));
	if (~meta.flag.print)
		title('Streamwise velocity');
		h = colorbar();
	else
		axis off
		figure(10000);
		clf();
		colormap(cmap);
		h = colorbar('southoutside');
		axis off
		caxis(ulim);
		h.Limits = ulim;
		%box off
	end
		h.Label.String='U (m/s)';
		h.Limits
end % sanggau_plot_velocity_nz
