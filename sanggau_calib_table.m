% 2014-08-19 17:56:45.433877678 +0200
% Karl Kastner, Berlin

% TODO, level is not the right name, it should read depth
% summarise calibration data into a table
	addpath ../discharge/mat

	folder   = '../discharge/mat/';
	base     = 'sanggau-calib-table';
	date     = datestr(now(),'yyyy-mm-dd-HH-MM');
	txtname  = [base,'-',date,'.txt'];
	texname  = [base,'-',date,'.tex'];
	txtname_ = [base,'.txt'];
	texname_ = [base,'.tex'];

	unitL = 'm';
	unitA = 'm^2';
	unitU = 'm/s';
	unitQ = 'm^3/s'

%	diary(txtname);
	fid = fopen([folder filesep txtname],'w');

	% load metadata
	meta = sanggau_metadata();                                                      
                                                                                
	% preload water level                                                           
	load([ROOTFOLDER,'/src/water-level/mat/water-level.mat']);
	level.val  = K(nid.sanggau_merged).slot.depth;
	level.time = K(nid.sanggau_merged).slot.time;
                                                                                
	% load vadcp calibrarion discharge                                           
	calib = load_vadcp_discharge(meta.filename.discharge,level);

	% load rc discharge
	rc      = load('../discharge/mat/sanggau-stage-discharge.mat');
	gravel  = load('../bottom/mat/sanggau-gravel-data.mat');
	bedform = load('mat/sanggau-bedform-data.mat');
	load('mat/sanggau-calibration-error.mat');
	load('mat/sanggau-wududn.mat');

	% correct for average depth
	depth     = calib.l0  + calib.lH;
	level.val = level.val + calib.lH;

	% TODO put into discharge class
	% turbulent (effective) Dean number
	calib.De_t = calib.Re_t.*sqrt(calib.l0/meta.Rc);

	ubarv = calib.cs.u0;

	% average 
	lg_z0 = calib.cs.ln_z0/log(10);

	for idx=1:length(calib.t0)
		% vadcp : q, h,A, ubar, $u_bar hadcp, E u_vadcp - u_hadcp, E (u_vadcp - u_hadcp)^2, n_ens, v_boat$	\\
		date_C{idx} = datestr(calib.t0(idx),'dd/mm/yy');
		ubar_h(idx) = NaN;
		q_h(idx)    = NaN;

	    Q(idx)      = calib.dis(idx).cs.discharge.total(1);
%	    Qmid(idx)   = calib.dis(idx).cs.exval.ex.mid();
%	    Qtop(idx)   = calib.dis(idx).cs.discharge.top;
%	    Qbot(idx)   = calib.dis(idx).cs.discharge.bottom;
%	    Qleft(idx)  = calib.dis(idx).cs.discharge.left;
%	    Qright(idx) = calib.dis(idx).cs.discharge.right;
	end
	% discharge table
%	Qmtb = Q - Qleft - Qright;
%	Qmid = Q - Qleft - Qright - Qbot - Qtop;

	fprintf(fid,'%16s','Date');		fprintf(fid,'%10s',date_C{:}); fprintf(fid,'\n');
	fprintf(fid,'%16s','Depth (m)');	fprintf(fid,'%10.1f',depth); fprintf(fid,'\n');
	fprintf(fid,'%16s','Area (m^2)');	fprintf(fid,'%10.0f',calib.cs.area); fprintf(fid,'\n');
	fprintf(fid,'%16s','ubar vadcp (m/s)');	fprintf(fid,'%10.2f',calib.cs.u0); fprintf(fid,'\n');
	fprintf(fid,'%16s','ubar hadcp (m/s)');	fprintf(fid,'%10.2f',ubar_h); fprintf(fid,'\n');
%	fprintf(fid,'%12s','q vadcp (m^3)');	fprintf(fid,'%10.0f',calib.cs.q0); fprintf(fid,'\n');
	fprintf(fid,'%16s','q hadcp (m^3)');	fprintf(fid,'%10.0f',q_h); fprintf(fid,'\n');
	fprintf(fid,'%16s','Q_total');  fprintf(fid,'%10.0f', Q);      fprintf(fid,'\n');
%	fprintf(fid,'%16s','Q_mtb');    fprintf(fid,'%10.0f', Qmtb);   fprintf(fid,'\n');
%	fprintf(fid,'%16s','Q_mid');    fprintf(fid,'%10.0f', Qmid);   fprintf(fid,'\n');
%	fprintf(fid,'%16s','Q_top');    fprintf(fid,'%10.0f', Qtop);   fprintf(fid,'\n');
%	fprintf(fid,'%16s','Q_bottom'); fprintf(fid,'%10.0f', Qbot);   fprintf(fid,'\n');
%	fprintf(fid,'%16s','Q_left');   fprintf(fid,'%10.0f', Qleft);  fprintf(fid,'\n');
%	fprintf(fid,'%16s','Q_right');  fprintf(fid,'%10.0f', Qright); fprintf(fid,'\n');

	fprintf(fid,'Settings:\n');
%	fprintf(fid,'jflag '); fprintf(fid,'%d, ',calib.dis(1).jflag);  	fprintf(fid,'\n')
	fprintf(fid,'dw ');    fprintf(fid,'%f, ',calib.dw);         	        fprintf(fid,'\n')
	fprintf(fid,'dz ');    fprintf(fid,'%f, ',calib.dz);          	fprintf(fid,'\n')
	fprintf(fid,'vfunc: '); fprintf(fid,'%s ',func2str(calib.dis(1).cs.vfunc));	fprintf(fid,'\n')
	fprintf(fid,'cross section coordinates: %d %d %d %d', ...	fprintf(fid,'\n')
		int32([calib.dis(1).cs.xlim(1), calib.dis(1).cs.xlim(2), calib.dis(1).cs.ylim(1), calib.dis(1).cs.ylim(2)]));	fprintf(fid,'\n')

	
	fclose(fid);

% machine readable table (latex)

fid = fopen([folder filesep texname],'w');


% for each campaign
l = 'ABCDE';
for idx=1:length(calib.t0)
	defprint(fid,['date',l(idx)],'%s', datestr(calib.t0(idx),'dd/mm/yyyy'));
	defprint(fid,['cvDepthAvg' l(idx)], '{%2.1f}',depth(idx),unitL);
	defprint(fid,['cArea' l(idx)], '{%1.2f \\, 10^3}',calib.cs.area(idx)/1e3,unitA);
	defprint(fid,['cvQ' l(idx)], '{%1.2f \\, 10^3}',calib.cs.q0(idx)/1e3,unitQ);
	defprint(fid,['cvUAvg' l(idx)], '{%1.2f}',calib.cs.u0(idx),unitU);
	defprint(fid,['cLgZ' l(idx)], '{%1.1f}',lg_z0(idx));
	defprint(fid,['cLnZ' l(idx)], '{%1.1f}',calib.cs.ln_z0(idx));
	defprint(fid,['cZ' l(idx)], '{%0.2f}',100*exp(calib.cs.ln_z0(idx)),'cm');
	defprint(fid,['cChezy' l(idx)], '{%2.0f}',calib.cs.Chezy(idx));
	defprint(fid,['cCfi' l(idx)], '{%2.0f}',1./sqrt(calib.cs.C_f(idx)));
	defprint(fid,['cDeant' l(idx)], '{%3.1f}',calib.De_t(idx));
	defprint(fid,['cRet' l(idx)], '{%2.0f}',calib.Re_t(idx));
	defprint(fid,['cWidth' l(idx)], '{%3.0f}',calib.cs.width(idx),unitL);
	defprint(fid,['cP' l(idx)], '{%3.0f}',calib.cs.perimeter(idx),unitL);
	defprint(fid,['wududn' l(idx)], '{%1.2f}',wududn.val(idx));
end
% Max calibration
	defprint(fid,'cvDepthAvgMax',  '%2.1f',max(depth),unitL);
	defprint(fid,'depthAvgBankful','%2.1f',max(depth)+meta.bankfull_offset,unitL);
	defprint(fid,'cAreaMax',       '%1.2f',max(calib.cs.area)/1e3,unitA);
	defprint(fid,'cvQMax',         '%1.2f \\, 10^3',max(calib.cs.q0)/1e3,unitQ);
	defprint(fid,'cvUAvgMax',      '%1.2f',max(calib.cs.u0),unitU);
% Min
	defprint(fid,'cvDepthAvgMin',  '%2.1f',min(depth),unitL);
	defprint(fid,'cAreaMin',       '%1.2f \\, 10^3',min(calib.cs.area)/1e3,unitA);
	defprint(fid,'cvQMin',         '%1.2f \\, 10^3',min(calib.cs.q0)/1e3,unitQ);
	defprint(fid,'cvUAvgMin',      '%1.2f',min(calib.cs.u0),unitU);
% Range
	defprint(fid,'cvDepthAvgRange','%2.1f',max(depth)-min(depth),unitL);

	defprint(fid,'zHadcpA',        '%2.1f',calib.lH,unitL);
	defprint(fid,'zHadcpB',        '%2.1f',calib.lH-meta.hadcp.redeploy.d,unitL);
	defprint(fid,'meshWidth',      '%1.1f',calib.dw,unitL);
	defprint(fid,'meshHeight',     '%1.2f',calib.dz,unitL);
% discharge and level for all time
	defprint(fid,'manningrcParam', '%0.2f',rc.manningrc.lin.param.val);
	defprint(fid,'manningrcRR',    '%0.3f',rc.manningrc.lin.R2);
	fdx = ~isfinite(rc.level);
	defprint(fid,'hMin',           '%2.1f',min(rc.level),unitL);
	defprint(fid,'hMax',           '%2.1f',max(rc.level),unitL);
	rc.level(fdx) = 0;
	q = rc.discharge;
	defprint(fid,'hMed',           '%2.1f',median(rc.level),unitL);
	defprint(fid,'qMin',           '%2.1f \\, 10^3',min(q)/1e3,unitQ);
	defprint(fid,'qMax',           '%2.1f \\, 10^3',max(q)/1e3,unitQ);
	q(fdx) = 0;
	defprint(fid,'qMed',           '%2.1f \\, 10^3',median(q)/1e3,unitQ);
	defprint(fid,'rcParam',        '%2.1f, %1.2f, %1.2f',rc.powerrc.lin.param.val);
% roughness data
	defprint(fid,'DuneLnZzero',  '%2.1f', bedform.ln_z0);
	defprint(fid,'DuneL',        '%2.1f', bedform.length);
	defprint(fid,'DuneH',        '%2.1f', bedform.height);
	defprint(fid,'GravelLnZzero','%2.1f', gravel.ln_z0);

% parameters of the error distribution of the predicted vertical profile
	defprint(fid, 'flnzSerrZero', '%0.3f',       sqrt(vertical.roughness(1).model.s2));
	defprint(fid, 'flnzSerrOne', '%0.3f',        sqrt(vertical.roughness(2).model.s2));
	defprint(fid, 'flnzSerrZeroDirect', '%0.3f', sqrt(vertical.direct(1).model.s2));
	defprint(fid, 'flnzSerrOneDirect', '%0.3f',  sqrt(vertical.direct(2).model.s2));
	defprint(fid, 'flnzRhoZero', '%0.3f',        sqrt(vertical.roughness(1).model.rho));
	defprint(fid, 'flnzRhoOne', '%0.3f',         sqrt(vertical.roughness(2).model.rho));
	defprint(fid, 'flnzRhoZeroDirect', '%0.3f',  sqrt(vertical.direct(1).model.rho));
	defprint(fid, 'flnzRhoOneDirect', '%0.3f',   sqrt(vertical.direct(2).model.rho));

% parameters of the error distribution of the predicted transversal profile
	defprint(fid, 'fnSerrZero', '%0.3f', sqrt(transversal(1).model.s2));
	defprint(fid, 'fnSerrOne', '%0.3f', sqrt(transversal(2).model.s2));
	defprint(fid, 'fnRhoZero', '%0.3f', sqrt(transversal(1).model.rho));
	defprint(fid, 'fnRhoOne', '%0.3f', sqrt(transversal(2).model.rho));

% cross correlation
	defprint(fid, 'rhoFzFnZero', '%0.2f', rho_nz(1));
	defprint(fid, 'rhoFzFnOne', '%0.2f', rho_nz(2));

% processing parameter
	defprint(fid,'nfRoughness','%2d',meta.nf_,unitL);
% deployment specific
	defprint(fid,'dredeploy','%2.1f',meta.hadcp.redeploy.d,unitL);

% roughness averages
	defprint(fid,'lnZzeroMean',      '%1.1f', roughness.mean);
	defprint(fid,'lnZzeroStd',       '%1.1f', roughness.std);
	defprint(fid,'lnZzeroInnerMean', '%1.1f', roughness.inner.mean);
	defprint(fid,'lnZzeroOuterMean', '%1.1f', roughness.outer.mean);
	
	

% close file
	fclose(fid);
	system(['ln -f -s ' txtname ' ' folder filesep txtname_]);
	system(['ln -f -s ' texname ' ' folder filesep texname_]);
	type([folder filesep txtname_]);

