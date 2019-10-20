% 2014-10-29 00:17:50 +0100
% Karl Kastner, Berlin

% compute energy slope during calibration

	addpath mat/
	sanggau_metadata();
	R = 20;

	H = [];
	us = [];
	for idx=1:length(filename_C)
		dis = load(filename_C{idx});
		discharge(idx) = dis.discharge;
		%us = [us discharge(idx).z0grid.val(:,1)];
		us_ = conv(discharge(idx).z0grid_A(1).z0grid.val(:,1),ones(R,1)/R,'same');
		us = [us us_];
		H = [H discharge(idx).bgrid.val(:,1)];
	end
%	H = sqrt(H);
%	us = sqrt(us);
	figure(1);
	clf();
	plot(H,us.^2,'.');
	hold on
	us = real(us);
	cc = colormap('lines');
	for idx=1:length(filename_C)
		fdx = find(isfinite(us(:,idx)));
		% constant
		A = H(fdx,idx);
		c = A\us(fdx,idx).^2;
		plot([0 15], [0 15]*c(1),'color',cc(idx,:));
%		% afine (not better)
%		A = [H(fdx,idx).^0 H(fdx,idx)];
%		c = A\us(fdx,idx).^2;
%		plot([0 15], c(1) + [0 15]*c(2),'--','color',cc(idx,:));
	end
	cX1 = discharge(1).bgrid.cX1();

	% joint regression : this is non-sense, as the energy slope is not the same at each campaign
	us_ = us(:);
	fdx = find(isfinite(us_));
	A = [H(:)];
	c = A(fdx,:) \ us_(fdx).^2;
	err = (A(fdx,:)*c-us_(fdx).^2);
	plot([0 15],[0 15]*c(1),'k');
	R2 = 1 - var(err)/var(us_(fdx).^2)
	% robust
	verr = (quantile(err,0.84) - quantile(err,0.16)).^2;
	vus2 = (quantile(us_(fdx).^2,0.84) - quantile(us_(fdx).^2,0.16)).^2;
	R2   = 1 - verr/vus2
	%plot([0 15],c(1)+[0 15]*c(2),'k');
	% quantile regression
	c = median(us_(fdx).^2 ./ H(fdx));
	plot([0 15],[0 15]*c(1),'k--');
	set(gcf,'PaperUnits','points','PaperPosition',256*[0 0 2*1.096 2])
	ylabel('u_*^2 (m/s)');
	xlabel('h (m)');
	pdfprint('img/sanggau-us2-h.eps')
	figure(11);
	clf();
	subplot(3,1,2);
	plot([sort(us_(fdx))],(1:length(err))/length(err));
	subplot(3,1,3);
	plot([sort(err)],(1:length(err))/length(err));
	figure(2);
	clf
	plot(discharge(1).bgrid.cX1,us.^2./(Constant.g*H));
	[my intercept slope] = robustlinreg(discharge(1).bgrid.cX1,us.^2./(Constant.g*H));
	X = discharge(1).bgrid.cX1;
	hold on
	plot(X,repmat(intercept(:)',length(X),1)+X(:)*slope(:)');
	xlim([min(cX1) max(cX1)]);
	xlabel('N (m)');
	ylabel('S_0')
	legend(datestr( ...
		[mean(discharge(1).ens.time) ...
		 mean(discharge(2).ens.time) ...
		 mean(discharge(3).ens.time) ...
		 mean(discharge(4).ens.time)], ...
		'dd/mm/yyyy'));
	title('energy slope');
set(gcf,'PaperUnits','points','PaperPosition',256*[0 0 3*1.096 2])
	pdfprint('img/energy-slope.eps');
	figure(3)
	plot(discharge(1).bgrid.cX1,H)
	(diff(quantile(us.^2./(Constant.g*H),[0.16 0.84]))/2)
	nanmedian(us.^2./H) 
	X = discharge(1).bgrid.cX1()';	
	[my intercept slope] = robustlinreg(discharge(1).bgrid.cX1,us(:).^2./H(:));
	(diff(quantile(us.^2./H-slope*[X X X X],[0.16 0.84]))/2)
	figure(4)
	plot(us.^2./H-slope*[X X X X])

