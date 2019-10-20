% Sun Aug 24 19:52:16 CEST 2014
% Karl Kastner, Berlin

	
%	bankfile = [ROOTFOLDER,filesep(),'dat/kapuas/gis/riverbank/kapuas-riverbank-2013line.shp'];
%	centrefile = [ROOTFOLDER,filesep(),'dat/kapuas/gis/centreline/kapuas-centreline-2013.shp'];

	if (~exist('pflag','var'))
		pflag = 0;
	end

	nc  = 8;
	
	meta = sanggau_metadata();
	dsi       = 10;
	
	% load bank line data
	bank   = Shp.read(meta.filename.bank);
	bank   = Shp.flat(bank);
	% cut to Sanggau region
	bank   = Shp.cut(bank, meta.utm.x, meta.utm.y, meta.ds_max);
	
	% load centre line data
	centreshp    = Shp.read( meta.filename.centre );

	% consider river before and behind the adcp within distance ds_max
	centreshp = Shp.cut(centreshp, meta.utm.x, meta.utm.y, meta.ds_max);
	centre    = Centreline(centreshp);

	% create S coordinate from XY coordinate
	centre.init();

	% resample streamwise resolution
	centre.resample(dsi, meta.centreline.Ri, meta.centreline.p );

	% compute channel curvature
	centre.curvature(); %[], meta.bank.Ri, meta.bank.p, meta.bank.iflag);
	%centre.curvature([], meta.bank.Ri, meta.bank.p, meta.bank.iflag);

	% convert bank line to S-N coordinates
	[bank] = centre.channel_planemetry(bank, meta.bank.Ri, meta.bank.p, meta.bank.smoothflag);

	% find centreline point closest to the HADCP
	s0 = centre.xy2sn(meta.utm.x, meta.utm.y);
	
	% channel midpoint (coordinate origin)
	D2 = (centre.X - meta.utm.x).^2 + (centre.Y-meta.utm.y).^2;
	[mv mdx] = min(D2);
	x0_ = centre.X(mdx);
	y0_ = centre.Y(mdx);
	
	% direction along centreline (TODO use poly-derivatives)
	dir = [ centre.X(mdx+1)-centre.X(mdx-1);
		centre.Y(mdx+1)-centre.Y(mdx-1)];
	% TODO sign depends how centreline was parsed from map
	dir = -dir/norm(dir);

	% direction orthogonal to the centreline
	odir = [-dir(2); dir(1)];
	
	Rc = 1./centre.curvature;
	xr = centre.X(mdx)-odir(1)*Rc(mdx); %R(mdx);
	yr = centre.Y(mdx)-odir(2)*Rc(mdx);
	
	% plot curvature
	namedfigure(1,'Radius of curvature and width');
	clf();

	subplot(2,1,1);
	% axis limits chosen between next two points where R has a pole
	ldx = mdx - find(abs(Rc(mdx:-1:1))>meta.R_max,1)
	rdx = mdx + find(abs(Rc(mdx:end)) > meta.R_max,1)
	pdx = (ldx:rdx)'; %$find(R<0);
	S_ = -(centre.seg_S(pdx)-centre.seg_S(mdx))
	
	plot(S_,sign(Rc(mdx))*Rc(pdx),'-k');
	grid on
	%q = quantile(R(pdx),0.1);
	%axis([S_(end) S_(1) 0 4000])
	xlabel('S (m)');
	ylabel('R_c (m)');
	
	subplot(2,1,2);
	width_ = centre.width_;
	plot(S_,width_(pdx),'-k');
	grid on
	%xlim([-2000 2000]);
	xlabel('S (m)');
	ylabel('w (m)');
	
	% plot topography
	namedfigure(2,'Bathymetry 2d');
	scale = 1e-3;
	% load bathymery
	load('../bathymetry/mat/sanggau-NaN-IPoly-10-1-1-2.mat');
	clf();
	patch_man(target.T,scale*(target.X-x0_),scale*(target.Y-y0_),-target.V,'edgecolor','none');
%	for idx=1:size(target.T,1);
%		patch(target.X(target.T(idx,:))-x0_,target.Y(target.T(idx,:))-y0_,target.V(target.T(idx,:)),'edgecolor','none');
%		hold on;
%	end
	hold on
%	colormap(colormap_man('BGR',[0.1 0.9]));
	handle = colorbar();
	ylabel(handle,'m');
%	dcolorbar(8,[-10 6],nc,colormap_man('BGR',[0.1 0.9]));
	%cmap = dcolormap(8,[-10 6],colormap_man('BGR',[0.1 0.9]));	
	dcolormap(8,[-10 6],colormap_man('BGR',[0.1 0.9]));	

	caxis([-10 6]);
	c=colormap;
	c=interp1(1:size(c,1),c,linspace(1,size(c,1),8));
	colormap(c);

	% centreline
	plot(scale*(centre.X-x0_),scale*(centre.Y-y0_),'k-');
	hold on
	% bankline
	plot(scale*(bank.X(:)-x0_), scale*(bank.Y(:)-y0_),'-k');
	plot(scale*(meta.utm.x-x0_),scale*(meta.utm.y-y0_),'k*');
	% tangent circle
	plot(scale*(xr-x0_),scale*(yr-y0_),'ko','markerfacecolor',[0 0 0]);
	curvature = centre.curvature;
	Rc = 1./curvature;
	circle(scale*(xr-x0_),scale*(yr-y0_),scale*abs(Rc(mdx)),'EdgeColor',[0 0 0],'LineWidth',2,'linestyle','--');
	
	% extended cross section
	plot(scale*([meta.utm.x xr]-x0_),scale*([meta.utm.y yr]-y0_),'k--');
	
	% flow arrow
	% TODO nearest function and dir function
	d = (centre.X - x0_ - 2750).^2 + (centre.Y - y0_ - 250).^2;
	[void fdx_] = min(d);
	%fdx_ = round((1-5/10)*length(S_))';
	dx   = -diff(centre.X);
	dy   = -diff(centre.Y);
	hyp  = hypot(dx,dy);
	hq    = quiver(scale*(centre.X(fdx_)-x0_),scale*(centre.Y(fdx_)-y0_), ...
			scale*-dx(fdx_)/hyp(fdx_), scale*-dy(fdx_)/hyp(fdx_), ...
			400,'k','linewidth',5,'MaxHeadSize',1,'AutoScaleFactor',1);
%	adjust_quiver_arrowhead_size(h,10);
	
	axis equal
	axis(scale*[-1000 +3000 -2500 +1500])
	grid on
	% setting xticklabel without xtick is dangerous,
	% because ticks, but not any more the labels, change during resizing
	
	xtick = -5:5; %get(gca,'xtick');
	ytick = -5:5; %get(gca,'ytick');
	set(gca,'xtick',xtick,'xticklabel',xtick);
	set(gca,'ytick',ytick','yticklabel',ytick);
	xlabel('Easting (km)');
	ylabel('Northing (km)');
	%annotation('arrow', 'Color', atan2(dy(fdx_),dx(fdx_)),'headStyle','cback1','HeadLength',100,'HeadWidth',10);
	%=get(h, 'Children');
	%Data = get(l,'XData');
	%Data = get(l,'YData');
	%et(h, 'XData', XData{2}, 'YData', 10*YData{2});
	
	if (pflag)
		pdfprint(1,'img/sanggau-curvature');
	
		for idx=1:3
			pdfprint(2,'img/sanggau-topography',idx);
%			system('convert -density 150 img/sanggau-topography.pdf img/sanggau-topography.png');
			if (1 == idx)
				system('convert -density 150 img/sanggau-topography-crop.pdf img/sanggau-topography-crop.png');
			else
				system(['convert -density 150 img/sanggau-topography-',num2str(idx),'-crop.pdf img/sanggau-topography-',num2str(idx),'-crop.png']);
			end
		%	pdfprint(2,'img/sanggau-topography',3);
		%	system('convert -density 150 img/sanggau-topography-3.pdf img/sanggau-topography-3.png');
		end
	end
	
	printf('R_cross R_min\n');
	%R = abs(centre.Rc);
	R = abs(Rc);
	[R(mdx) min(R(pdx))]
	printf('W_min W_cross W_max\n');
	%W = centre.width;
	W = centre.width_;
	[min(W(pdx)) W(mdx) max(W(pdx))]

	% left and right bank coordinates at the HADCP
if (0)
	xl = interp1(left.S,left.X,s0,'linear');	
	yl = interp1(left.S,left.Y,s0,'linear');
	xr = interp1(right.S,right.X,s0,'linear');	
	yr = interp1(right.S,right.Y,s0,'linear');
	xc = 0.5*(xl+xr);
	yc = 0.5*(yl+yr);
	h = sqrt((xl-xr)^2+(yl-yr)^2);
	sin_ = (yr-yl)/h;
	cos_ = (xr-xl)/h;

	plot([xl xr]-x0_,[yl yr]-y0_,'-o')

	load('../hadcp/mat/hadcp-sanggau-squeezed.mat');
	hadcp  = HADCP(hadcp);

	% correct the depth for redeployment
	fdx = find(hadcp.time > meta.hadcp.redeploy.t);
	offset = zeros(size(hadcp.idepth_m));
	offset(fdx) = meta.hadcp.redeploy.d;
	level = hadcp.idepth_m - offset;

	for idx=1:length(meta.filename.discharge)
		dis_ = load(['mat' filesep filename_C{idx}],'discharge');
		dis(idx) = dis_.discharge;
		t0(idx,1)       = mean(datenum(dis(idx).adcp.time));
		l0(idx)         = interp1(hadcp.time,level,t0(idx));
		dis(idx).ens.H  = dis(idx).ens.H - l0(idx);
	end

	% almost bankful
	lmax = l0(1);
	[-l0 + lmax]

	figure(3);
	clf();
	cc = colormap('lines');
	Nl = [];
	Hl = [];
	Nr = [];
	Hr = [];
	for idx=1:length(meta.filename.discharge)
		N  = (dis(idx).ens.X - xc)*cos_ + (dis(idx).ens.Y - yc)*sin_;
%		S  = (dis(idx).ens.X - xc)*cos_ - (dis(idx).ens.Y - yc)*sin_;
%		figure()
%		plot(S,N)
%pause
		D2 = (dis(idx).ens.X - xl).^2 + (dis(idx).ens.Y - yl).^2;
		fdx = find(D2 < 75*75);
		subplot(1,2,1);
		plot(N(fdx),-(dis(idx).ens.H(fdx)+lmax),'color',cc(idx,:))
		%plot(dis(idx).ens.X(fdx),dis(idx).ens.H(fdx),'color',cc(idx,:))
		hold on
		Nl = [Nl; N(fdx)];
		Hl = [Hl; dis(idx).ens.H(fdx)];

		D2 = (dis(idx).ens.X - xr).^2 + (dis(idx).ens.Y - yr).^2;
		fdx = find(D2 < 75*75);
		subplot(1,2,2);
		plot(N(fdx),-(dis(idx).ens.H(fdx)+lmax),'color',cc(idx,:))
		%plot(dis(idx).ens.X(fdx),dis(idx).ens.H(fdx),'color',cc(idx,:))
		hold on
		Nr = [Nr; N(fdx)];
		Hr = [Hr; dis(idx).ens.H(fdx)];
	end
	A = [Nl.^3 Nl.^2 Nl Nl.^0];
	c = A \ -(Hl+lmax);
	r = roots(c)
        subplot(1,2,1);                                                 
	N = linspace(min(Nl),max(Nl),100)';
	A =[N.^3 N.^2 N N.^0];
	plot(N,A*c,'k');

	A = [Nr.^3 Nr.^2 Nr Nr.^0];
	% TODO, has to be normalised
	c = A \ -(Hr+lmax);
	r = roots(c)
        subplot(1,2,2);                                                 
	N = linspace(min(Nr),max(Nr),100)';
	A =[N.^3 N.^2 N N.^0];
	plot(N,A*c,'k');

	% the quadratics do not work well, so visually estimated
	nl =  -350
%	xl = -sin_*nl + xc
%	yl =  cos_*nl + yc
	xl =  cos_*nl + xc
	yl =  sin_*nl + yc
	nr =  315
%	xr = -s*nr + xc
%	yr =  cos_*nr    + yc
	xr =  cos_*nr + xc
	yr =  sin_*nr + yc
	figure(2)
	plot([xl xr]-x0_,[yl yr]-y0_,'-og');

end

