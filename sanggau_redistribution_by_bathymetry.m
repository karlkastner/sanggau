% Fr 19. Jun 11:41:56 CEST 2015
% Karl Kastner, Berlin

	namedfigure(2,'Velocity distribution predicted by bathymetry');
	clf();
	hold on
	plot(N,filtfunc(double(urel.bathy),nf),'-');
	ylim([0.6 1.2])
	% predicted velocity distribution by bathymetry
	% this section is deprecated
	nn = size(bottom_A,1);
	l0 = calib.l0;
	h_A = bottom_A + repmat(rvec(l0),nn,1);
	u_A = (log(h_A)-1-ln_z0).*(h_A.^0.5) .* (h_A > 0);
	q_A = u_A.*h_A;
	U = (nansum(q_A)./nansum(h_A));
	urel.bathy = u_A./repmat(rvec(U),nn,1);

	namedfigure(4,'Depth during calibration');
	plot(N,h_A);

	hold on
	urel.bathy     = filtfunc(urel.bathy,nf);
	m = mean(urel.bathy, 2);
	s = std(urel.bathy, [], 2);
	plot(N,m,'g'); hold on
	plot(N,m*[1 1] + s*[-1 1],'g--');

