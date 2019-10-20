% 2014-08-28 17:03:55.877756197 +0200
% Karl Kastner, Berlin

% estimate z_0 with decreasing mesh width in n-direction

	figure(3);
	clf();

	% quick fix
	method = 3; %[2 3 4 5];
	dn = 2.^(-(2:0.25:7));

	for bdx=1:2

	if (1 == bdx)
		id   = 1;
		tdx  = 1
	
		ofolder = 'mat';
		ofilename_C = {
		'2013-12-09-sanggau', ...
		'2014-02-20-sanggau', ...
		'2014-04-18-sanggau', ...
		'2014-06-18-sanggau' };
	
		ofilename = [ofolder filesep() ofilename_C{id} '-bart.mat'];
		load(ofilename);
	
		vel  = msh(tdx).sec.velocity;
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
		%bin.S   = bin.Z_ ./ zzbed;
		for idx=1:nens
			ldx_ = find(bin.S(:,idx)>0,1,'last');;
			if (isempty(ldx_))
				ldx_=0;
			end
			ldx(idx,1:2) = ldx_*[1 1]; 
		end
	
		% pseudo bottom grid
		bgrid.cX1     = nbed;
		bgrid.cX1.val = zbed;
	
		% pseudo velocity grid
		vgrid.val = msh.sec.velocity;
	
		dn_ = abs(nbed(2)-nbed(1));
		dn = dn_*(1:40)';
	else
		file_C = { ...
		'mat/2013-12-09-sanggau.mat', ...
		'mat/2014-02-20-sanggau.mat', ...
		'mat/2014-04-18-sanggau.mat', ...
		'mat/2014-06-18-sanggau.mat' };
		load(file_C{1});

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
	
%		dn_ = abs(N(2)-N(1))*cs.width;
		dn = dn_*(1:40)';
	end % else bdx

	for idx=1:length(dn)
	for mdx=1:length(method)
	idx

	if (4 ~= mdx)
	[cs z0grid] = regress_log_profile( ...
			vel, ...
			Z, ...
			N, ...
			ldx, ...
			cs, ...
			dn(idx), ...
			method(mdx) ...
			);
	z0grid_A(idx,mdx,bdx).cs = cs;
	z0grid_A(idx,mdx,bdx).z0grid = z0grid;
	else
		z0grid.val = 0.5*( z0grid_A(idx,1,bdx).z0grid.val + z0grid_A(idx,2,bdx).z0grid.val );
	end

	% interpolate u_s and z_0 back to the main grid
	u_s = interp1(z0grid.cX1,z0grid.val(:,1),N,'constant')
	ln_z_0 = interp1(z0grid.cX1,z0grid.val(:,2),N,'constant')
	u_s = real(u_s);
	ln_z_0 = real(ln_z_0);
%	vz(idx,mdx) = var(ln_z_0);
	vz(idx,mdx) = diff(quantile(ln_z_0(:),[0.75 0.25]))*0.5/norminv(0.25);
%	vus(idx,mdx) = var(u_s);
	vus(idx,mdx) = diff(quantile(u_s(:),[0.75 0.25]))*0.5/norminv(0.25);
		
	uu_s    = repmat(u_s,nbin,1);
	lln_z_0 = repmat(ln_z_0,nbin,1);
	U_reg = (1/Constant.KAPPA)*uu_s.*(log(Z) - lln_z_0);
	U_reg = real(U_reg);
	Ur{idx,mdx,bdx} = U_reg;
	Zr{idx,mdx,bdx} = Z;
	Nr{idx,mdx,bdx} = N;
	% difference
	U_diff = U_reg - vel;
	bias(idx,mdx) = nanmedian(U_diff(:));
	serr(idx,mdx) = diff(quantile(U_diff(:),[0.75 0.25]))*0.5/norminv(0.25);
	l_(idx,mdx) = median(ln_z_0);
	u_(idx,mdx) = median(u_s);

%	complicated
	if (2 == bdx)
%		dl =	
%		du = 
%		N1 = repmat(Nr{idx,mdx,1}(:)',nbin,1);
%		n2 = repmat(N(:)',nbin,1);
%		U1 = interp2(N1,Zr{idx,mdx,1},Ur{idx,mdx,1},NN,ZZ,'linear');
%		U_diff = U_reg - 1;
%		cbias(idx,mdx) = nanmedian(U_diff(:));
%		cerr(idx,mdx)  = diff(quantile(U_diff(:),[0.75 0.25]))*0.5/norminv(0.25);
	end
	end % mdx
	end % idx
figure(bdx);
clf();
subplot(3,1,1);
loglog(dn,abs(bias),'.');
subplot(3,1,2);
loglog(dn,serr(:,:),'.');
subplot(3,2,5)
semilogy(dn,vz); ylim([1 1e2])
subplot(3,2,6)
semilogy(dn,vus); ylim([1 1e2])

u_
l_
% semilogy(abs(fft(z0grid_A(1,2).z0grid.val(4:end-2,2)/log(10))),'.')
figure(3);
color = {'b', 'r'};
hold on
subplot(2,1,1)
color = {'b', 'r'};
plot(dn,bias,['.' color{bdx}]);
hold on
subplot(2,1,2)
plot(dn,serr,['.' color{bdx}]);
%serr_C(bdx) = serr;
%dn_C(bdx) = dn;
	end
