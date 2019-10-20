% Thu Mar 26 20:41:22 CET 2015
% Karl Kastner, Berlin

	if (~exist('reload','var') || reload)
		meta  = sanggau_metadata();
		[hadcp hlevel offset] = sanggau_load_hadcp(meta);
%		hadcp = sanggau_load_hadcp(meta,[],1800);
		load('../bathymetry/mat/sanggau-NaN-IPoly-10-1-1-2.mat');
		reload = 0;
	end
	if (exist('pflag','var'))
		pflag = 0;
	end
	ylim_ = 1e3*[-2.5 1.5];
	xlim_ = 1e3*[-1 3];
	d = 500;
	% xlim_ = [-2.5 2.5];	
	% ylim_ [-4000 1000])
	%xlim(1e3*[-2.75 2.25]);

	X = target.X - meta.utm.x;
	Y = target.Y - meta.utm.y;
if (1)
        namedfigure(1,'Bathymetry Sanggau');
        clf();                                                                  

	T = target.T;
	V = target.V;
	patch(X([T(:,1) T(:,2) T(:,3) T(:,4)]'), ...
		Y([T(:,1) T(:,2) T(:,3) T(:,4)]'), ...
		V([T(:,1) T(:,2) T(:,3) T(:,4)]'), 'edgecolor','none');

        axis equal                                                              
        colormap('jet');
	axis equal
	xlim(xlim_)
	ylim(ylim_);
%	caxis([0 35]);
	clim_ = [-5 10];
	caxis(clim_);
	hold on;
	plot(0,0,'k*');
	heading = median(hadcp.heading_rad);
	dir = [sin(heading); cos(heading)];
	hold on;
	s=450;
	plot([0 s*dir(1)],[0 s*dir(2)],'k');
	xlabel('Northing (m)');
	ylabel('Easting (m)');
	colorbar();
	grid on
	set(gca,'xtick',xlim_(1):d:xlim_(2));
	set(gca,'ytick',ylim_(1):d:ylim_(2));
	if (pflag)
		pdfprint('img/sanggau-deployment-map.pdf')
		system('convert -density 150 img/sanggau-deployment-map.pdf img/sanggau-deployment-map.png');
		pdfprint(1,'img/sanggau-deployment-map.pdf',2)
		system('convert -density 150 img/sanggau-deployment-map-2.pdf img/sanggau-deployment-map-2.png');
		pflag = 0;
	end
%	namedfigure(2,'');
%	clf();
%	X = target.X - meta.utm.x;
%	Y = target.Y - meta.utm.y;
%       for idx=1:size(target.T,1);                                     
%		patch(X(target.T(idx,:)),Y(target.T(idx,:)),target.V(target.T(idx,:)),'edgecolor','none');
%                hold on;                                                
%        end 
%	colormap gray
%	colorbar
%	caxis([-3 0])	
end

