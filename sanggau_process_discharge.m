% 2014-11-14 16:35:49.838417817 +0100
% Karl Kastner, Berlin

% TODO preload also VADCP
% TODO merge with rasau_process_discharge and bifurcation_process discharge
function discharge = sanggau_process_discharge(opt,dataset_id)

	if (nargin() < 1 || isempty(opt))
		opt = sanggau_metadata;
	end

	iname_C = opt.filename.vadcp;

	if (nargin() < 2 || isempty(dataset_id))
		dataset_id = 1:length(iname_C);
	end


	dw = opt.dw;

	% load water level for surface level correction
	load(opt.filename.water_level);
	level      = struct;	
	level.time = Kr(nid.sanggau_merged).time;
	level.val  = Kr(nid.sanggau_merged).depth;

	load(opt.filename.bscalibration);

	% process adcp data for each measurement day
	for id=rvec(dataset_id)
		fprintf(1,'Processing %s\n',iname_C{id});
		base = iname_C{id}(1:end-4);
		iname = iname_C{id};
		% load vadcp data	
		load(iname);


		% TODO preload ?
		vadcp = VADCP(vadcp, ... % 'utm.zone',opt.utm.zone, ...
				'd_transducer', opt.vadcp.d_transducer);
		vadcp.velocity_profile = opt.velocity_profile;
		vadcp.sediment_concentration_profile = opt.sediment_concentration_profile;
		%vadcp.sediment_concentration_profile.nlflag = false;
		vadcp.backscatter_obj  = backscatter_calibration.bso(end);


		%level(idx).val = level(idx).level;
		vadcp.fill_coordinate_gaps();
		vadcp.assign_water_level(level.time,level.val);

		% define transect
		for jdx=1
		transect(jdx) = ADCP_Transect(...
			 'tdx',  jdx ...
			,'mode', 'returning' ... %mtransect.mode ...
			,'xlim', [opt.cs.left(1), opt.cs.right(1)] ...
			,'ylim', [opt.cs.left(2), opt.cs.right(2)] ...
			);
		end % for jdx

		% side-effect free assigning of transects
		% has to precede detect_crossings
		vadcp.assign_transect(transect);
	
		for jdx=1:length(transect)
			transect(jdx).detect_crossings(vadcp);
			% export mmt-file
			base_ = basename(dirname(dirname(iname)));
	
			mkdir(['mat/',base_,'/']);
			mmtname = ['mat/',base_,'/',base_,'.mmt'];
			% TODO no magic names
			pd0name = ['C:\Users\pia\Desktop\',base_,'\',base_,'.PD0'];
			transect(jdx).export_mmt(mmtname,pd0name,vadcp);
		end % for jdx

		load(opt.filename.bscalibration);
		vadcp.velocity_profile = opt.velocity_profile;
		vadcp.sediment_concentration_profile = opt.sediment_concentration_profile;
		vadcp.backscatter_obj  = backscatter_calibration.bso(end);
		% no fines, only sensed coarse sediment
		vadcp.backscatter_obj.param(3) = 0;
		vadcp.shear_stress_field = opt.shear_stress_field;
		vadcp.process(transect);

		% transect.integrate_discharge(vadcp);

		% for different grid widths
		% TODO better lambda for vmethod = v1
		for jdx=1:length(dw)
			cs   = CrossSection( ...
	                          'dw', dw(jdx) ...
			        , 'dz', opt.dz ...
				, 'bmethod', opt.bmethod ...
				, 'cdx', 1 ...
				, 'errmode', 'none' ...
				, 'lambda', opt.lambda ...
				, 'level', level ...
			        , 'topofbank', opt.cs.topofbank ...
				, 'transect', transect ...
			        , 'T_max', opt.cs.T_max ...
		                , 'vdimension',  opt.vdimension ...
				, 'grid_n',      opt.grid_n ...
		                , 'vmethod2',  opt.vmethod2 ...
			);
			% TODO return vadcp
			cs.process_discharge( vadcp );
			cs.process_backscatter( vadcp );
			discharge = cs.discharge(cs.t0(1));
			area      = cs.area(cs.t0(1));
			radius    = cs.radius(cs.t0(1));
			oname = [base, cs.optstr(cs), '.mat'];
			fprintf('Writing %s\n',oname);
			save(oname, 'cs','discharge','area','radius');
			clear cs
		end % for jdx
		save([base,'-processed.mat'], 'vadcp');
		clear vadcp
	end % for id
	dt = [];
	discharge_summary(opt,'sanggau',dt);
end % sanggau_process_discharge

