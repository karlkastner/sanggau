% 2014-08-28 17:03:55.877756197 +0200
% Karl Kastner, Berlin

% estimate z_0 with decreasing mesh width in n-direction
% TODO, let winfilt ignore nans instead of patching it

function sanggau_z_0_campaigns()
	addpath ../bottom/
	addpath	../discharge/
	bottom = sanggau_load_gsd();
%	bottom_coordinates();

	% plot range
	prange = [-3 0];

	% quick fix
	method = 3;
	% choose methods to use
	B = 1; %[1 2 3];
	% spatial z_0 resolution
	dn = 1;
	% filter options
	fn     = 120/dn;
	forder = 2; 

%	pcolor =   [[1 0.5 0.5];
  %                  [0.5 0.5 1];
 %                   [0.5 1 0.5];
	%lcolor = {'r','b','g'};
	lcolor = colormap('lines')
%	hsv_ = rgb2hsv(lcolor);
%	hsv_(:,3) = min(1,hsv_(:,3) + 0.5);
%	scolor = hsv2rgb(hsv_)\
	scolor = min(1,lcolor+0.25);
%pause

	folder = '../discharge/mat';
	karl.name_C = { ...
		'2013-12-09-sanggau-jflag-1-R_z0-8', ...
		'2014-02-20-sanggau-jflag-0-R_z0-8', ...
		'2014-04-18-sanggau-jflag-1-R_z0-8', ...
		'2014-06-18-sanggau-jflag-1-R_z0-8' };

%		2014-04-18-sanggau-jflag-1-R_z0-8.mat
%		'2013-12-09-sanggau', ...
%		'2014-02-20-sanggau', ...
%		'2014-04-18-sanggau', ...
%		'2014-06-18-sanggau' };
	bart.name_C = {
		'2013-12-09-sanggau-bart', ...
		'2014-02-20-sanggau-bart', ...
		'2014-04-18-sanggau-bart', ...
		'2014-06-18-sanggau-bart' };

	% prepare figures
	for idx=1:10
		figure(idx);
		clf();
	end

	% for each campaign
	nc = length(karl.name_C);
	for cdx=1:nc
		leg{cdx} = karl.name_C{cdx}(1:10);

	% for each data source (bart or mine)
	for bdx=B

	if (bdx > 1)
		% bart's method
		tdx  = 1
	
		filename = [folder filesep() bart.name_C{cdx} '.mat'];
		load(filename);

		if (2 == bdx)
			vel  = msh(tdx).sec.velocity;
		else
			vel = msh(tdx).sec.vele;
		end
		nbed = msh(tdx).p.nbed;
		zbed = msh(tdx).p.zbed;
		nbin = size(vel,1);
		nens = size(vel,2); 
	
		% pseudo cross section
		width = (nbed(end)-nbed(1));
		cs.width = width;
	
		% mid-grid cell
		nbed = nbed(2:2:end-1);
		% zbed is also errornously defined, it is -depth
		zbed = -zbed(2:2:end-1);
	
		% pseudo ensembles
		N   = msh(tdx).N(1,:)/cs.width;
	
		% pseudo bins
		vel = vel(:,:,1);
		zzbed = repmat(zbed(:)',nbin,1);
		% "Z" is misnamed, it is actually -d = z - H
		Z  = zzbed + msh(tdx).Z;
		S   = Z ./ zzbed;
		for idx=1:nens
			ldx_ = find(S(:,idx)>0,1,'last');;
			if (isempty(ldx_))
				ldx_=0;
			end
			ldx(idx,1:2) = ldx_*[1 1]; 
		end
	else	% karl's method
		load([folder filesep karl.name_C{cdx} '.mat']);

		cs = discharge.cs;
		vgrid = discharge.vgrid;
		ldx = NaN(vgrid.n1-1,2);
		vel = flipud(vgrid.val(:,:,1)');
		nbin = size(vel,1);
		for jdx=1:vgrid.n1-1
			fdx = find(isfinite(vel(:,jdx)),1,'last');
			if (isempty(fdx))
				ldx(jdx,1:2) = [0 0];
			else
				ldx(jdx,1:2) = fdx*[1 1];
			end
		end
		Z = repmat(discharge.bgrid.val(:,1)',vgrid.n2-1,1) ...
			+ flipud(vgrid.cXX2');
		N = vgrid.cX1;
	end % else bdx

	z0grid = discharge.z0grid;

	if (bdx>1)
		% was flipud, why?
		u_s_   = z0grid.val.u_s(:,1);
		ln_z0_ = z0grid.val.ln_z0(:,1)/log(10);
	else
		u_s_   = z0grid.val.u_s(:,1);
		ln_z0_ = z0grid.val.ln_z0(:,1)/log(10);
	end
	u_s(:,cdx,bdx)       = u_s_;
	ln_z0(:,cdx,bdx)     = ln_z0_;
	ln_z0_err(:,cdx,bdx) = z0grid.err.ln_z0(:,1)/log(10);

	% filter the roughness length
	l = ln_z0(:,cdx,bdx);
	[f e] = winfilt(l,fn,forder);
	f_ln_z0.val(:,cdx,bdx) = f;
	f_ln_z0.err(:,cdx,bdx) = e;

	% compute the friction factor
	H = discharge.bgrid.val(:,1);
	fr(:,cdx,bdx) = (log(H) - 1 - ln_z0(:,cdx,bdx)).^(-2);

	end % for bdx (each data source)
	end % for cdx (each campaign)


	X = z0grid.cX1();
	X = X(:);
	n = length(X);
	% for each data source, average over campaigns
	for bdx=B

		% shear velocity plot
		namedfigure(1,'Shear velocity spatial distribution');
		%subplot(2,1,bdx);
		%q1 = quantile(squeeze(u_s(:,:,bdx)'),[0.16 0.5 0.84])';
		m1 = mean(u_s(:,:,bdx),2);
		e1 = serr(u_s(:,:,bdx)')';
		q1 = [m1-e1 m1 m1+e1];

		cdata = double(zeros(1,n-1,3));
		cdata(:,:,1) = scolor(bdx,1);
		cdata(:,:,2) = scolor(bdx,2);
		cdata(:,:,3) = scolor(bdx,3);
		patch(  [X(1:end-1) X(2:end) X(2:end) X(1:end-1)]', ...
			[q1(1:end-1,1) q1(2:end,1) q1(2:end,3) q1(1:end-1,3)]', ...
			cdata, 'edgecolor','none', 'facealpha', 0.25,'facecolor','flat');
		hold on
		plot(X,q1(:,2),'-','color',lcolor(bdx,:),'linewidth',2);
		ylim(quantile(q1(:,2),[0.01 0.99]))
		xlim([X(1) X(end)])

		namedfigure(2,'Average roughness length spatial distribution (unfiltered)');
		m1 = mean(ln_z0(:,:,bdx),2);
		e1 = serr(ln_z0(:,:,bdx)')';
		q1 = [m1-e1 m1 m1+e1];
		cdata = double(zeros(1,n-1,3));
		cdata(:,:,1) = scolor(bdx,1);
		cdata(:,:,2) = scolor(bdx,2);
		cdata(:,:,3) = scolor(bdx,3);
		patch(  [X(1:end-1) X(2:end) X(2:end) X(1:end-1)]', ...
			[q1(1:end-1,1) q1(2:end,1) q1(2:end,3) q1(1:end-1,3)]', ...
			cdata, 'edgecolor','none', 'facealpha', 0.25,'facecolor','flat');
		hold on
		plot(X,q1(:,2),'-','color',lcolor(bdx,:),'linewidth',2);
		ylim([-4 -1]);
		xlim([X(1) X(end)])
%		ylim(quantile(q1(:,2),[0.01 0.99]))
		xlim([X(1) X(end)])

		namedfigure(3,'Weighed average of roughness length');
		%m1 = mean(ln_z0(:,:,bdx),2);
		%e1 = serr(ln_z0(:,:,bdx)')';
		[m1 s2] = minavg(ln_z0(:,:,bdx),ln_z0_err(:,:,bdx).^2);
		e1 = sqrt(s2);
		q1 = [m1-e1 m1 m1+e1];
		cdata = double(zeros(1,n-1,3));
		cdata(:,:,1) = scolor(bdx,1);
		cdata(:,:,2) = scolor(bdx,2);
		cdata(:,:,3) = scolor(bdx,3);
		patch(  [X(1:end-1) X(2:end) X(2:end) X(1:end-1)]', ...
			[q1(1:end-1,1) q1(2:end,1) q1(2:end,3) q1(1:end-1,3)]', ...
			cdata, 'edgecolor','none', 'facealpha', 0.25,'facecolor','flat');
		hold on
		plot(X,q1(:,2),'-','color',lcolor(bdx,:),'linewidth',2);
%		ylim(quantile(q1(:,2),[0.01 0.99]))
		xlim([X(1) X(end)])
		ylim([-4 -1]);
		xlim([X(1) X(end)])

	
		namedfigure(4,'Roughness length spatial distribution, averaged before filtering');
		q = nanmedian(ln_z0(:,:,bdx),2);
		% filter
		fdx = find(~isfinite(q));
		q(fdx) = nanmedian(q);
		qf         = q;
		[qf sqf]   = winfilt(qf,fn,forder);
		qf_(:,bdx) = qf;
		qf = [qf-sqf,qf,qf+sqf];
		patch(   [X(1:end-1) X(2:end) X(2:end) X(1:end-1)]',...
			 [qf(1:end-1,1) qf(2:end,1) qf(2:end,3) qf(1:end-1,3)]', ...
			 cdata, 'edgecolor','none', 'facealpha', 0.25,'facecolor','flat');
		hold on
		plot(X,qf(:,2),'-','linewidth',2,'color',lcolor(bdx,:));
		plot(X,q,'.','color',lcolor(bdx,:));
		ylim([-4 -1]);
		xlim([X(1) X(end)])
		%ylim([-5 0]);
		xlim([X(1) X(end)])
	
		namedfigure(5,'Roughness length spatial distribution, averaged after filtering');
		%subplot(2,1,1);
	%	subplot(length(B),1,bdx);
		[m1 s2] = minavg(f_ln_z0.val(:,:,bdx),f_ln_z0.err(:,:,bdx).^2);
		e1 = sqrt(s2);
		qf = [m1-e1 m1 m1+e1];
		patch(  [X(1:end-1) X(2:end) X(2:end) X(1:end-1)]', ...
			[qf(1:end-1,1) qf(2:end,1) qf(2:end,3) qf(1:end-1,3)]', ...
			cdata, 'edgecolor','none', 'facealpha', 0.25,'facecolor','flat');
		hold on
		plot(X,qf(:,2),'-','linewidth',2,'color',lcolor(bdx,:));
		title('Spatial distribution of filtered roughness length');
		ylim([-4 -1]);
		xlim([X(1) X(end)])

		
		namedfigure(6,'Roughness length spatial distribution, weighed average');
		% average weighed with local error
		[mopt s2opt] = minavg(f_ln_z0.val(:,:,bdx),f_ln_z0.err(:,:,bdx).^2);
		sopt = sqrt(s2opt);
		qf = [mopt-sopt mopt mopt+sopt];
		patch(  [X(1:end-1) X(2:end) X(2:end) X(1:end-1)]', ...
			[qf(1:end-1,1) qf(2:end,1) qf(2:end,3) qf(1:end-1,3)]', ...
			cdata, 'edgecolor','none', 'facealpha', 0.25,'facecolor','flat');
		hold on
		plot(X,mopt,'-','linewidth',2,'color',lcolor(bdx,:));
		%ylim([-5 0]);
		xlim([X(1) X(end)])
		ylim([-4 -1]);
		xlim([X(1) X(end)])
		title('Spatial distribution of optimally filtered roughness length');
	
		figure(7);
%		subplot(length(B),2,2*(bdx-1)+1);
		plot(X,ln_z0(:,:,bdx));
		ylim(quantile(reshape(ln_z0(:,:,bdx),[],1),[0.01,0.99]));
		xlim([X(1) X(end)])
		title('Spatial distribution of unfiltered roughness length');
		ylim([-4 -1]);
		xlim([X(1) X(end)])
%		subplot(length(B),2,2*(bdx-1)+2);
%		plot(X,ln_z0.err(:,:,bdx));
%		ylim(quantile(reshape(ln_z0.err(:,:,bdx),[],1),0.01,0.99));
		xlim([X(1) X(end)])
%		title('Spatial distribution of error');

		figure(8);
		subplot(length(B),2,2*(bdx-1)+1);
		plot(X,f_ln_z0.val(:,:,bdx));
		%ylim(quantile(reshape(f_ln_z0.val(:,:,bdx),[],1),[0.01,0.99]));
		xlim([X(1) X(end)])
		title('Spatial distribution of filtered roughness length');
		ylim([-4 -1]);
		xlim([X(1) X(end)])
		subplot(length(B),2,2*(bdx-1)+2);
		plot(X,f_ln_z0.err(:,:,bdx));
		%ylim(quantile(reshape(f_ln_z0.err(:,:,bdx),[],1),[0.01,0.99]));
		xlim([X(1) X(end)])
		title('Spatial distribution of filter error');
		ylim([-4 -1]);
		xlim([X(1) X(end)])

		namedfigure(9,'CDF of shear velocity');
		clf();
		q4 = u_s(:,:,bdx);
		fdx = find(isfinite(sum(q4,2)));
		n4 = length(fdx);
		q4 = sort(q4(fdx,:));
		plot(q4,(1/n4)*(1:n4)'*ones(1,nc));
		xlim(quantile(q4(:),normcdf([-2 2])));
		xlabel('u_z')
		legend('location','northwest',leg{:});
	
		namedfigure(10,'CDF of roughness length');
		clf();
		q4 = ln_z0(:,:,bdx);
		fdx = find(isfinite(sum(q4,2)));
		n4 = length(fdx);
		q4 = sort(q4(fdx,:));
		plot(q4,(1/n4)*(1:n4)'*ones(1,nc));
		xlim(quantile(q4(:),normcdf([-2 2])));
		xlabel('log_{10} z_0')
		legend('location','northwest',leg{:});

		namedfigure(11,'CDF of friction parameter');
		clf();
		q4 = fr(:,:,bdx);
		fdx = find(isfinite(sum(q4,2)));
		n4 = length(fdx);
		q4 = sort(q4(fdx,:));
		plot(q4,(1/n4)*(1:n4)'*ones(1,nc));
		xlim(quantile(q4(:),normcdf([-2 2])));
		xlabel('f')
		legend('location','northwest',leg{:});
		
		namedfigure(12,'CDF of roughness length standard error');
		clf();
		q4 = ln_z0_err(:,:,bdx);
		fdx = find(isfinite(sum(q4,2)));
		n4 = length(fdx);
		q4 = sort(q4(fdx,:));
		plot(q4,(1/n4)*(1:n4)'*ones(1,nc));
		xlim(quantile(q4(:),normcdf([-2 2])));
		xlabel('log_{10} z_0')
		legend('location','northwest',leg{:});

	
	end % for bdx
	
	% plot sediment grain size into functions
	figure(3);
	r = 182:185;
	bottom.ln_z_0(r) = log10(0.1*2.^bottom.d84(r));
	% TODO, this should be done by a function !
	c = discharge.cs.dir(1);
	s = discharge.cs.dir(2);
	R = [c s; -s c];
	NT = (R*[(bottom.X - discharge.xpc) (bottom.Y - discharge.ypc)]')';
	plot(NT(r,1),bottom.ln_z_0(r),'ko','markerfacecolor','k');
%	ylim(prange);
	xlim([min(discharge.bgrid.cX1()) max(discharge.bgrid.cX1())]);
	ylabel({'$log_{10}(z_0)$',' '},'interpreter','latex');
	xlabel({' ','$m$'},'interpreter','latex');
	%print('img/sanggau-roughness-length.svg','
	%pdfprint('img/sanggau-roughness-length.eps');
	plot2svg('img/sanggau-roughness-length.svg',2);
	system('LD_LIBRARY_PATH= inkscape img/sanggau-roughness-length.svg --export-pdf=img/sanggau-roughness-length.pdf');

	figure(9);
	pdfprint('img/sanggau-cdf-shear-velocity.eps');

	figure(11);
	pdfprint('img/sanggau-cdf-friction-parameter.eps');

end % functions
