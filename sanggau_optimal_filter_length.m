% Fr 19. Jun 11:13:54 CEST 2015
% Karl Kastner, Berlin

% TODO different filter type
% TODO better handle boundaries

	if (~exist('reload','var') || reload)
		meta = sanggau_metadata();	
		% load hadcp data
		[hadcp hlevel offset] = sanggau_load_hadcp(meta);
		% load water level
		load([ROOTFOLDER,'/src/water-level/mat/water-level.mat']);
		level.val  = K(nid.sanggau_merged).slot.depth;
		level.time = K(nid.sanggau_merged).slot.time;

		% load vadcp data
		calib = load_vadcp_discharge(meta.filename.discharge,level);
%		level.val = level.val + calib.lH;
		reload = 0;
	end

	cc = colormap('lines');
	for idx=1:calib.nc
		meas.U(idx,1) = calib.dis(idx).cs.U(1);
		meas.u(:,idx) = calib.dis(idx).Un;
		% interpolate missing values with no slip bc
		N = calib.dis(1).gridN.cX1;
		meas.u(1,:)   = 0;
		meas.u(end,:) = 0;
		for jdx=1:size(meas.u,2)
			fdx = isnan(meas.u(:,jdx));
			meas.u(fdx,jdx) = interp1(cX1(~fdx),meas.u(~fdx,jdx),cX1(fdx),'linear');
		end
		%meas.urel(:,idx) = meas.u(:,idx)./calib.dis(idx).cs.U(1);
		%meas.urel = bsxfun(@times,meas.u,1./mean(meas.u));
	end % for idx

	% influence of filter length
	nf = 2*round(2.^(0.5:0.25:7))';
	fdx = abs(N) < max(N) - 0.25*nf(end);
	clear s2tot s2dis
	for idx=1:length(nf)
		%[mu stot_ sdis_] = meanfilt1(meas.urel,nf(idx));
		%[mu stot_ sdis_] = polyfilt1(meas.urel,nf(idx),2);
		%[mu stot_ sdis_] = meanfilt1(meanfilt1(meanfilt1(meanfilt1(meas.urel,nf(idx)/4),nf(idx)/4),nf(idx)/4),nf(idx)/4);
		%[mu stot_ sdis_] = filter1(meas.urel,nf(idx),0);
		[mu stot_ sdis_] = filter1(meas.u,nf(idx),1);
		s2tot(idx,:) = mean(stot_(fdx,:).^2);
		s2dis(idx,:) = mean(sdis_(fdx,:).^2);
	end
	s2noise = s2tot - s2dis;
	namedfigure(1,'Standard error of streamwise velocity vs. filter length');
	clf();
	subplot(2,1,1);
	loglog(nf,sqrt(s2tot),'.-');
	legend(datestr(cvec(calib.t0),'dd/mm/yy'));
	subplot(2,1,2);
	loglog(nf,sqrt(s2dis),'*-');
	hold on
	set(gca, 'ColorOrderIndex', 1)
	loglog(nf,sqrt(s2noise),'o-');

