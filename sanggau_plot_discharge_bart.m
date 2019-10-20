% Tue Aug 26 17:15:49 CEST 2014
% Karl Kastner, Berlin

% process sanggau velocity with bart's software
function [tool discharge] = sanggau_plot_discharge_bart(id)
% TODO quick fix
	rev = 100;
	ofolder = '../discharge/mat';

	dz = '0.25';
	dw = '1';	

	opt = '1-gps';
%	opt = ['dz-' dz '-dw-' dw '-pt-1-gps';
	opt = ['dz-' dz '-dw-' dw '-pt-1-bt'];
	opt_ = ['-dw-' dw '-dz-' dz '-discharge'];

	filename_C = {
	'2013-12-09-sanggau', ...
	'2014-02-20-sanggau', ...
	'2014-04-18-sanggau', ...
	'2014-06-18-sanggau' };

	filename = [ofolder filesep() filename_C{id} '-bart-' opt '.mat'];
	
	load(filename);
	tool.compQ();

	meta = sanggau_metadata();	
	iname = [meta.filename.discharge{id}(1:end-4) opt_ '.mat']
	load(iname);
        discharge.integrate_discharge();

        field = 'sec';
	Q1 = tool.msh(1).(field).Q;
	Q2 = discharge.cs.discharge;
	A1 = tool.msh(1).(field).A;
	A2 = discharge.cs.area;
	Q1.mtb   = Q1.mid+Q1.bot+Q1.top;
	Q2.mtb   = Q2.total - Q2.left - Q2.right;
	Q2.mid   = Q2.mtb - Q2.bottom - Q2.top;
	A1.mtb   = A1.mid+A1.bot+A1.top;
	A2.mtb   = A2.total - A2.left - A2.right;
	A2.mid   = A2.mtb - A2.bottom - A2.top;

	title('Bart Karl	')
	datestr(discharge.adcp.t0,'dd/mm/yyyy')

	fprintf('A_total  '); fprintf('%5.0f %5.0f', [A1.total A2.total(1)]); fprintf('\n');
	fprintf('A_mtb    '); fprintf('%5.0f %5.0f', [A1.mtb A2.mtb]); fprintf('\n');
	fprintf('A_mid    '); fprintf('%5.0f %5.0f', [A1.mid A2.mid]); fprintf('\n');
	fprintf('A_top    '); fprintf('%5.0f %5.0f', [A1.top A2.top]); fprintf('\n');
	fprintf('A_bottom '); fprintf('%5.0f %5.0f', [A1.bot A2.bottom]); fprintf('\n');
	fprintf('A_left   '); fprintf('%5.0f %5.0f', [A1.left  A2.left]); fprintf('\n');
	fprintf('A_right  '); fprintf('%5.0f %5.0f', [A1.right A2.right]); fprintf('\n');

	fprintf('Q_total  '); fprintf('%5.0f %5.0f', [Q1.total Q2.total(1)]); fprintf('\n');
	fprintf('Q_mtb    '); fprintf('%5.0f %5.0f', [Q1.mtb Q2.mtb(1)]); fprintf('\n');
	fprintf('Q_mid    '); fprintf('%5.0f %5.0f', [Q1.mid Q2.mid(1)]); fprintf('\n');
	fprintf('Q_top    '); fprintf('%5.0f %5.0f', [Q1.top Q2.top]); fprintf('\n');
	fprintf('Q_bottom '); fprintf('%5.0f %5.0f', [Q1.bot Q2.bottom]); fprintf('\n');
	fprintf('Q_left   '); fprintf('%5.0f %5.0f', [Q1.left  Q2.left]); fprintf('\n');
	fprintf('Q_right  '); fprintf('%5.0f %5.0f', [Q1.right Q2.right]); fprintf('\n');

	% number of transect (there is only 1)
	tdx = 1;
	
	% velocity
	veli = tool.Unz();

	% patch boundaries
	N       = tool.msh(tdx).p.N;
	Z       = tool.msh(tdx).p.Z;
	nbed    = tool.msh(tdx).p.nbed;
	zbed    = tool.msh(tdx).p.zbed;
	
	namedfigure(1,'Velocity NZ');
	clf();
	for idx=1:3
		if (rev <= 20)
			fgood   = msh(tdx).p.fgood;
			fdx	= find(fgood);
		else % rev 45
			fdx   = tool.msh(tdx).p.fgood_3(:,idx);
		end
		subplot(3,2,2*idx-1)
		% velocity
		% TODO somehow N/Z changed from vec to grid, which version
		patch(N,Z,veli(fdx)','linestyle','none');
		colorbar();
		hold on
		% bottom
		plot(nbed,zbed,'color',[0.5 0.5 0.5],'linewidth',2);
		xlim([nbed(1) nbed(end)])
		subplot(3,2,2*idx)
		hist_man(veli(fdx),100);
		%quiver(msh(ct).N,msh(ct).Z,msh(ct).cs.velocity(:,:,2)*qscale,msh(ct).velocity(:,:,3)*qscale,'autoscale','off','color','k')
	end

	colormap jet

	namedfigure(2,'Velocity N');
	Nb  = tool.msh(tdx).N;
	Un = tool.Un();
	plot(Nb,Un);

	% TODO histogram for vel and vele difference (not large)

	namedfigure(3,'Bottom profile');
	clf();
	plot(nbed-0.5*(nbed(end)-nbed(1)),zbed,'k');
	hold on
	plot(discharge.gridN.cX1,-discharge.gridN.val.bottom,'b')
	ylim([-15 0]);
	fdx = isfinite(discharge.adcp.N4.*discharge.adcp.T4.*discharge.adcp.H4);
	interp = TriScatteredInterp(discharge.adcp.N4(fdx),discharge.adcp.T4(fdx),double(discharge.adcp.H4(fdx)));
	b = interp(discharge.gridN.cX1,0.*discharge.gridN.cX1);
	plot(discharge.gridN.cX1,-b,'r');
	legend('Bart', 'Karl 1D', 'Karl 2D');

	% compare to my method
	% N
	cX1 = discharge.gridNZ.cX1();
	% Z
	cX2    = discharge.gridNZ.cX2();
	NN     = repmat(rvec(cX1),length(cX2),1)';
	ZZ     = repmat(cvec(cX2),1,length(cX1))';
	fdx  = tool.msh(tdx).p.fgood_3(:,1);
	N      = [tool.msh(tdx).N(:) ];
%                  tool.msh(tdx).Ntop
%                  tool.msh(tdx).Ntop
%                  tool.msh(tdx).Ntop
	N      = N - 0.5*(max(tool.msh.N(:)) - min(tool.msh.N(:)));
	Z      = tool.msh(tdx).Z;

	veli_ = -tool.VELnz();
	f = {'U','V','W','Mag'};
	N_ = discharge.gridNZ.cX1();
	Z_ = discharge.gridNZ.cX2();
	for idx=1:4
	%figure(100+idx);
	namedfigure(100+idx,['Velocity differences ' f{idx}]);
	clf();
	if (idx < 4)
	U_   = discharge.gridNZ.val.(f{idx})';
	veli = veli_(:,:,idx);
	else
	U_   = sqrt(discharge.gridNZ.val.(f{1}).^2 + discharge.gridNZ.val.(f{2}).^2)';
	veli = sqrt(veli_(:,:,1).^2 + veli_(:,:,2).^2);
	end
	interp = TriScatteredInterp(N(fdx),Z(fdx),veli(fdx));
	U  =  interp(NN(:),ZZ(:));
	U  = reshape(U,length(cX1),length(cX2))';
	subplot(3,1,1);
	imagesc(N_,Z_,U);
	axis xy;
	title('Bart');
	colorbar();
	caxis([0 1.5]);
	subplot(3,1,2);
	imagesc(N_,Z_,U_);
	axis xy;
	title('Karl');
	colorbar();
	caxis([0 1.5]);
	subplot(3,1.5,4);
	imagesc(N_,Z_,(U-U_));
	axis xy;
	colorbar();
	caxis([-1 1]);
	title('Bart - Karl');
	subplot(3,3,9);
	hist_man(U(:)-U_(:));
	colormap jet
	end

	% plot mesh
	figure(1000);
	clf();
	N = tool.msh.p.N';
	Z = tool.msh.p.Z';
	T = [1,2,6; 2,5,6; 2,3,5; 3,4,5];
	A = zeros(size(N,1),6);
	for idx=1:4
		A(:,idx) = triarea([N(:,T(idx,1)), N(:,T(idx,2)), N(:,T(idx,3))], ...
		                   [Z(:,T(idx,1)), Z(:,T(idx,2)), Z(:,T(idx,3))] );
	end
	max(A(:))
	min(A(:))
	sum(A(:))
	A = sum(A,2);
	Q = tool.msh.sec.vele(tool.msh.p.fgood)'*A

%	for idx=1:size(N,1)
%		for jdx=1:7
%			plot(N(idx,jdx),Z(idx,jdx),'.');
%			text(N(idx,jdx),Z(idx,jdx),num2str(jdx));
%			hold on
%		end
%		%patch([N(idx,1) N(idx,2) N(idx,3)
%	end
%	namedfigure(5,'Velocity difference histogram');
%	fdx = find(isfinite(U-U_));

end % sanggau_bart

