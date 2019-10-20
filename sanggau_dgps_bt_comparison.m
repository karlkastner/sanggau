% 2014-08-31 22:50:28 +0200
% Karl Kastner, Berlin

% TODO investigate bias dependend on spatioal distribution

function dis = sanggau_dgps_bt_comparison()

	addpath ./mat

	extended = 0;	
	
	% outlier cut off
	lim = [-1 1];
	% 2d-histogram range +-
	scale = 0.5;

	% filter options
	fn	= 0;%7;
	forder	= 0;%2;

	% magnetic deviation
	% http://www.geomag.nrcan.gc.ca, http://www.ngdc.noaa.gov/geomag-web/
	alpha_m = deg2rad(0.66);
	Rm = [ cos(alpha_m) sin(alpha_m);
	      -sin(alpha_m) cos(alpha_m)];

	filename_C = { ...
	'2013-12-09-sanggau.mat', ...
	'2014-02-20-sanggau.mat', ...
	'2014-04-18-sanggau.mat', ...
	'2014-06-18-sanggau.mat' };
	clear serr s

	edges = linspace(lim(1),lim(2),100);
	for idx=1:length(filename_C)
		% load data
		load(filename_C{idx});
		dis(idx) = discharge;
		ubar(idx,:) = discharge.cs.ubar;
		ens      = discharge.ens;
		btvel = ens.btvel(:,1:2);
		btvela = ens.btvela(:,1:2);
		X = ens.X;
		Y = ens.Y;

		% compensate for magnetic deviation
		btvel = (Rm*btvel')';

		% rotate bottom velocity from earth to cross section reference
		[void Ri] = rotvel(discharge.cs.dir,[]);
		btvel = (Ri*btvel')';
		btvela = (Ri*btvela')';
		
		% filter
%		btvela(:,1) = winfilt(btvela(:,1),fn,forder);
%		btvela(:,2) = winfilt(btvela(:,2),fn,forder);

		% dgps-vadcp velocity std and bias
		du = btvel(:,1:2)-btvela(:,1:2);
		du__ = du;

		% remove outliers
		fdx = find(min(du,[],2) > lim(1) & max(du,[],2) < lim(2));
		btvel  = btvel(fdx,:);
		btvela = btvela(fdx,:);
		du = du(fdx,:);
		X = X(fdx);
		Y = Y(fdx);
		du_ = du;
	
		if (extended)
		figure(1);
		% 2d scatter
		subplot(2,2,idx);
		plot(du(:,1),du(:,2),'.','markersize',1);
		axis(scale*[-1 1 -1 1]);
		axis square
		end

		figure(2);
		% 2d greyscale histogram
		subplot(2,2,idx);
		fdx = find(isfinite(du__(:,1).*du__(:,2)));
		e = linspace(-scale,scale,101);
		h = bin2d(du__(fdx,1),du__(fdx,2),ones(length(fdx),1),e,e,@sum);
		h = h/length(fdx);
		imagesc(e,e,1-h);
		axis square
		colormap gray
		axis(scale*[-1 1 -1 1]);
		set(gca,'xtick',[-0.5 -0.25 0 0.25 0.5]);
		set(gca,'ytick',[-0.5 -0.25 0 0.25 0.5]);
		grid on
		xlabel('\Delta u_{stream} (m/s)');
		ylabel('\Delta u_{span} (m/s)');

		% time series	
		figure(10+idx)
		for jdx=1:2;
			subplot(2,1,jdx);
			plot([btvel(:,jdx),btvela(:,jdx)],'.','markersize',1);
			ylim([-4 4]);
		end
	
		% velocities
		% TODO norm
		figure(3)
		subplot(2,2,idx)
		hist([btvel(:,1),btvela(:,1)],1000);
		title('Velocity');
	
		% magnitude quotient
		figure(4);
		nv = sqrt(btvel(:,1).^2 + btvel(:,2).^2);
		nva = sqrt(btvela(:,1).^2 + btvela(:,2).^2);
		qn = nv./nva;
		fdx_ = find(isfinite(qn));
		subplot(2,2,idx);
		hist(qn(fdx_),linspace(0,2,100))
		gdx = find(qn > 2.^-0.25 & qn < 2.^0.25);
		gdx = fdx_;
		hist(qn(gdx),linspace(0,2,100));
		xlim([0 2]);
		du     = du(gdx,:);
		btvel_ = btvel;
		btvel  = btvel(gdx,:);
		btvela = btvela(gdx,:);
		X = X(gdx);
		Y = Y(gdx);
		title('Velocity magnitude quotient');
	
		% bias and standard error
		q = quantile(du, [normcdf(-1) 0.5 normcdf(+1)])';
		bias(idx,:)   = mean(du);
		s(idx,:)      = std(du);
		bias_q(idx,1:2) = q(:,2);
		s_q(idx,:)     = 0.5*(q(:,3)-q(:,1));
		% effective number of samples
		C1 = corrcoef(du(1:end-1,1),du(2:end,1));
		C2 = corrcoef(du(1:end-2,1),du(3:end,1));
		rho(idx,1) = C1(2,1);
		rho(idx,2) = C2(2,1);
		f(idx,1) = sqrt((1+rho(idx,1))/(1-rho(idx,1)));
		n(idx,1) = (size(du,1)-1);
		serr(idx,:) = f(idx,1)/sqrt(n(idx,1))*s(idx,:);
	
		figure(5);
		subplot(2,2,idx);
		h = hist(du,edges);
		h = h/sum(h(:,1));
		bar(edges,h);
		xlim(lim);
		title('Velocity difference');
	
		% angle difference
		figure(6);
		subplot(2,2,idx);
		dphi = atan2(btvel(:,2),btvel(:,1)) - atan2(btvela(:,2),btvela(:,1));
		fdx = find(dphi>pi);
		dphi(fdx)=dphi(fdx)-2*pi;
		fdx = find(dphi<-pi);
		dphi(fdx)=dphi(fdx)+2*pi;
		dphi = rad2deg(dphi);
		hist(dphi,linspace(-0.125*180,0.125*180,100));
		bias_phi(idx,1) = mean(dphi);
		s_phi(idx,1) = std(dphi);
		bias_phi_q(idx,1) = nanmedian(dphi);
		title('Angle difference');
		xlim([-0.125*180 0.125*180]);

		% effective number of samples
		C_p = corrcoef(dphi(1:end-1,1),dphi(2:end,1));
		%C2 = corrcoef(du(1:end-2,1),du(3:end,1));
		rho_p(idx,1) = C_p(2,1);
		%rho(idx,2) = C2(2,1);
		f_phi(idx,1) = sqrt((1+rho_p(idx,1))/(1-rho_p(idx,1)));
		n_phi(idx,1) = (size(dphi,1)-1);
		serr_phi(idx,:) = f_phi(idx,1)/sqrt(n_phi(idx,1))*s_phi(idx,:);
	
		if (extended)

		% dgps velocity std
	%	du = [0.5*btvela(3:end,1:2) - btvela(2:end-1,1:2) + 0.5*btvela(1:end-2,1:2) ];
	%	du = 0.5*(-1*btvel(5:end,:) + 2*btvel(4:end-1,:) + 0 -2*btvel(2:end-3,:) + 1*btvel(1:end-4,:));
		du = (1*btvel(5:end,:)-4*btvel(4:end-1,:)+6*btvel(3:end-2,:)-4*btvel(2:end-3,:)+1*btvel(1:end-4,:));
		q = quantile(du,[normcdf(-1) normcdf(1)])';
		s_d(idx,:) = 0.5*(q(:,2)-q(:,1));
		s_d_(idx,:) = std(du);
		%rho_d(idx,:) = median(du.^2)./median(du(1:end-3,:).*du(4:end,:));
		rho_d(idx,:) = (nanmedian(du(1:end-3,:).*du(4:end,:))-(nanmedian(du(1:end-3,:)).*nanmedian(du(4:end,:)))) ...
				./ s_d(idx,:).^2;
		figure(7);
		subplot(2,2,idx);
		h = hist(du,edges);
		h = h/sum(h(:,1));
		bar(edges,h);
		xlim(lim);
		title('ADCP vel noise');
	
		% vadcp velocity std
		%du  = [0.5*btvel(3:end,1:2) - btvel(2:end-1,1:2) + 0.5*btvel(1:end-2,1:2) ];
	%	du = (1/12)*(-1*btvela(5:end,:) + 16*btvela(4:end-1,:) -30*btvela(3:end-2,:) + 16*btvela(2:end-3,:) - 1*btvela(1:end-4,:));
	%	du = 0.5*(-1*btvela(5:end,:) + 2*btvela(4:end-1,:) + 0 -2*btvela(2:end-3,:) + 1*btvela(1:end-4,:));
	%	du = (1*btvela(5:end,:)-4*btvela(4:end-1,:)+6*btvela(3:end-2,:)-4*btvela(2:end-3,:)+1*btvela(1:end-4,:));
		du = (0.5*btvela(5:end,:)-0.5*btvela(4:end-1,:)+0.5*btvela(3:end-2,:)-1*btvela(2:end-3,:)+1*btvela(1:end-4,:));
		q   = quantile(du,[normcdf(-1) normcdf(1)])';
		s_a(idx,:) = 0.5*(q(:,2)-q(:,1));
		%rho_a(idx,:) = nanmedian(du.^2)./nanmedian(du(1:end-3,:).*du(4:end,:));
		%rho_a(idx,:) = nanmedian(du(1:end-3,:).*du(4:end,:))./(nanmedian(du(1:end-3,:)).*nanmedian(du(4:end,:)));
		rho_a(idx,:) = (nanmedian(du(1:end-3,:).*du(4:end,:))-(nanmedian(du(1:end-3,:)).*nanmedian(du(4:end,:)))) ...
				./ s_a(idx,:).^2;
	
		figure(8);
		subplot(2,2,idx);
		h = hist(du,edges);
		h = h/sum(h(:,1));
		bar(edges,h);
		xlim(lim);
		title('DGPS vel noise');
		end % extended
	
	end
	%serr = sqrt(sum(serr.^2,2))
	'u_span u_stream b_span b_stream serr_span serr_stream n/1e4 f'
	[ubar bias serr n/1e4 f]
	'b_phi(deg) serr_phi n/1e4 f_phi'
	[bias_phi serr_phi n_phi/1e4 f_phi]
	%s_d
	%s_d_
	%rho_d
	%s_a
	%rho_a
	%C1
	% TODO, for C2 the dependence on C1 has to be removed first
	%C2	

	figure(2);
	pdfprint(['img/bt-vel-bias-' num2str(fn) '-' num2str(forder) '.eps']);
	figure(6);
	pdfprint(['img/bt-angle-bias-' num2str(fn) '-' num2str(forder) '.eps']);
	
end % function


