% Tue Nov 11 10:59:15 CET 2014
% Karl Kastner, Berlin

%v=1e-3*sqrt(sum((squeeze(double(tools.adcp.VEL(12,:,1:3)))'-double(tools.adcp.btvel(:,1:3)')).^2,1))'; fdx = find(v < 10 & mean(tools.adcp.btrange,2) > 400); median(v(fdx)), mean(v(fdx)),plot(sort(v(fdx))), serr(v)
function tool = sanggau_process_discharge_bart(id,pv,ShipReference)

	meta = sanggau_metadata();

	if (nargin() < 1 || isempty(id))
		id = 1:length(meta.filename.vadcp);
	end
	if (nargin() < 2 || isempty(pv))
		% procTrans version
		pv = 1;
	end

	ofolder = '../discharge/mat';

	% TODO get from ifilename
	ofilename_C = {
	'2013-12-09-sanggau', ...
	'2014-02-20-sanggau', ...
	'2014-04-18-sanggau', ...
	'2014-06-18-sanggau' };


	% TODO delta z 0.25 is too small, error with Bart's sw
	deltaZ          = 0.25;
	deltaN          = 1;
%	deltaZ          = 1;
%	deltaN          = 4;
	DepthTransducer = 0.25;
	if (nargin < 3 || isempty(ShipReference))
		%ShipReference = 'gps';
		ShipReference = 'bt';
	end
	%ShipReference   = 'bt';

	Pusr = 1e5*[   4.546551162508343   4.542045272175349
                       0.122617224468175   0.127176287124800];
	UseExtHeading = false;
%	varg = {};
	varg = {'Pusr', Pusr};
	
	% load data (preloaded with ADCPtools)
	for idx=id
		fprintf('%d\n',idx);
		load(meta.filename.vadcp{idx});

		if (2 == idx)
			% quick fix for corrupted nFiles data
			vadcp.nFiles.GGA.lat = nanmean(vadcp.NMEAGGA.Lat,2);
			vadcp.nFiles.GGA.long = nanmean(vadcp.NMEAGGA.Long,2);
		end

		tool = ADCPTools(vadcp);
		tool.ConventionalProcessing = true;

		% coordinates of beam intersection
		% TODO necessary ?
		if (0)
			[x,y]	   = utmADCP(adcp);
			[dx,dy,dz] = depthADCP(adcp);
			dx	   = bsxfun(@plus,dx,x');
			dy	   = bsxfun(@plus,dy,y');
			[mx,my,mz] = mapADCP(adcp);
		end
	
		% split data - for now one transect
		tid = ones(1,size(vadcp.VEL,2));
	
		% generate a mesh and process the velocity data
		if (1 == pv)
		%msh = toolprocTrans(adcp, tid, ...
		    tool.procTrans(tid, ...
			'DeltaZ', deltaZ, ...
			'DeltaN', deltaN, ...
			'DepthTransducer', DepthTransducer, ...
			'ShipReference', ShipReference, ...
			'UseExtHeading', UseExtHeading, ...
			varg{:} );
		else
		    tool.procTrans2(adcp,tid, ...
			'DeltaZ', deltaZ, ...
			'DeltaN', deltaN, ...
			'DepthTransducer', DepthTransducer, ...
			'ShipReference', ShipReference, ...
			varg{:} );
		end

		% compute discharge
		tool.compQ();

		%ofilename = [ofolder filesep() ofilename_C{idx} '-bart-',num2str(pv),'-',ShipReference,'.mat'];
		ofilename = [ofolder filesep() ofilename_C{idx} '-bart-dz-',num2str(deltaZ),'-dw-',num2str(deltaN),'-pt-',num2str(pv),'-',ShipReference,'.mat'];
		save(ofilename,'tool');
	end % for idx

end % function sanggau_bart_process

