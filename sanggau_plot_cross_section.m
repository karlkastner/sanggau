	namedfigure(8,'HADCP Beams');
	clf();
        hold on;
	% plot bottom profile
        plot(Nb,-(calib.bottom_A-calib.lH));
	
	% TODO, this starts at first cell not at 0
	s0 = Nh(1);
	se = Nh(end);
	r  = (se-s0)*tan(hadcp.SPREADANGLE_RAD);

	hold on
	cc = colormap('lines');
	for idx=1:length(calib.t0)
		theta = interp1(hadcp.time, hadcp.pitch_rad, calib.t0(idx));
		theta = sign(Nh(end)-Nh(1))*theta;
		c = cos(theta);
		s = sin(theta);
		R = [c s;
	             -s c];
		C = ( R*[0 (se-s0) (se-s0);
	                 0 +r      -r] )';
		%fdx = find(isfinite(lh0(:,idx)),1,'last');
		%patch([s0 se se],   [-lh0(1,idx) -lh0(fdx,idx) -lh0(fdx,idx)]+[0 -r r],[1 1 1],'facecolor','none');
		patch(s0+C(:,1)',lH-lh0(1,idx)+C(:,2)',[1 1 1],'facecolor','none','edgecolor',cc(idx,:));
	end
%	patch([s0 se se],dd+[0 -r r],[1 1 1],'facecolor','none');
	% plot water level at individual campaigns
	plot([Nb(1)*[1 1 1 1]; Nb(end)*[1 1 1 1]],[1 1]'*(calib.l0(:) + lH)');
	xlim([Nb(1) Nb(end)]);
	ylim([-3 13])

