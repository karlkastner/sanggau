% Thu Aug 21 08:51:21 CEST 2014
% Karl Kastner, Berlin

% TODO one could also rotate the profiles locally
% TODO regress profile parameter

meta = sanggau_metadata();

if (~exist('reload','var') || reload)

	% targe coordinate
	n = 100;
	S = (1:n)'/(n+1);
	imethod = 'pchip';
	
	load([ROOTFOLDER,'/src/water-level/mat/water-level.mat']);
	level.val  = K(nid.sanggau_merged).slot.depth;
	level.time = K(nid.sanggau_merged).slot.time;

	calib = load_vadcp_discharge(meta.filename.discharge,level);
	
	% for each data set
	filename_C = meta.filename.discharge;
	Uval = NaN(n,length(filename_C),'single');
	Ustd = NaN(n,length(filename_C),'single');
	Vval = NaN(n,length(filename_C),'single');
	Vstd = NaN(n,length(filename_C),'single');
	parfor idx=1:length(filename_C)
		% load data
		discharge = load(filename_C{idx});
		discharge = discharge.discharge;
		%nbins = max(discharge.adcp.nbins);
		U2d_ = NaN(n,discharge.gridN.n1-1,'single');
		V2d_ = NaN(n,discharge.gridN.n1-1,'single');
		S_ = discharge.adcp.bin.S;
		U_ = discharge.adcp.velocity.cs(:,:,1);
		V_ = discharge.adcp.velocity.cs(:,:,2);
		% for each slot in the cross section
		for jdx=1:discharge.gridN.n1-1
			id = discharge.gridN.id.i1.id(jdx).id;
			if (length(id) > 0)
			Uj = NaN(n,length(id),'single'); 
			Vj = NaN(n,length(id),'single'); 
			% interpolate profiles to common grid
			for kdx=1:length(id)
				if (sum(isfinite(U_(:,id(kdx)))) > 1)
					Uj(:,kdx) = interp1(S_(:,id(kdx)),U_(:,id(kdx)),S,imethod,NaN);
				end
				if (sum(isfinite(V_(:,id(kdx)))) > 1)
					Vj(:,kdx) = interp1(S_(:,id(kdx)),V_(:,id(kdx)),S,imethod,NaN);
				end
			end % for kdx
			% average profiles
				U2d_(:,jdx) = nanmean(Uj,2);
				V2d_(:,jdx) = nanmean(Vj,2);
			end
		end % for jdx
		U2d{idx} = U2d_;
		V2d{idx} = V2d_;
		% average profiles over entire cs
		% TODO area (momentum) weighed
		w = hanwin(calib.N/range(calib.N));
% TODO, make this a wmean/wstd function
		Uval(:,idx) = nansum(bsxfun(@times,U2d_,rvec(w)),2) ...
				./ nansum(bsxfun(@times,isfinite(U2d_),rvec(w)),2);
		Uvar(:,idx) = nansum(bsxfun(@times, ...
				bsxfun(@minus,U2d_,Uval(:,idx)).^2,rvec(w)),2) ...           
                                ./ (nansum(bsxfun(@times,isfinite(U2d_),rvec(w)),2)-1);
		Uverr(:,idx) = nansum(bsxfun(@times, ...
				bsxfun(@minus,U2d_,Uval(:,idx)).^2,rvec(w)),2) ...           
                                ./ (nansum(bsxfun(@times,isfinite(U2d_),rvec(w)),2)-1).^2;

		fdx = isfinite(nansum(V2d_));
		Vval(:,idx) = nansum(bsxfun(@times,V2d_,rvec(w)),2) ...
				./nansum(bsxfun(@times,isfinite(V2d_),rvec(w)),2);
		Vvar(:,idx) = nansum(bsxfun(@times, ...
				bsxfun(@minus,V2d_,Vval(:,idx)).^2,rvec(w)),2) ...           
                                ./ (nansum(bsxfun(@times,isfinite(V2d_),rvec(w)),2)-1);
		Vverr(:,idx) = nansum(bsxfun(@times, ...
				bsxfun(@minus,V2d_,Vval(:,idx)).^2,rvec(w)),2) ...           
                                ./ (nansum(bsxfun(@times,isfinite(V2d_),rvec(w)),2)-1).^2;

% todo, this has to be computed after rotation
		UVval(:,idx) = nansum(bsxfun(@times,V2d_.*U2d_,rvec(w)),2) ...
				./nansum(bsxfun(@times,isfinite(V2d_),rvec(w)),2);
	end % for idx
	reload = 0;
end % if reload
Ustd = sqrt(Uvar);
Vstd = sqrt(Vvar);
Userr = sqrt(Uverr);
Vserr = sqrt(Vverr);

% fit velocities to the surface
U     = Uval;
V     = Vval;
UV    = UVval;
%V_(fdx,:) = (S(fdx)-0.5)*c';
fdx = S>0.50 & S<=0.85;
A = vander_1d(S(fdx),1);
p_u = A\U(fdx,:);
p_v = A\V(fdx,:);
p_uv = A\UV(fdx,:);
fdx = S>0.85;
A = vander_1d(S(fdx),1);
U(fdx,:) = A*p_u;
V(fdx,:) = A*p_v;
UV(fdx,:) = A*p_uv;

% rotate velocities
if (0)
smid = 0.5; % 1/exp(1)
umid = interp1(S,Uval,smid);
vmid = interp1(S,Vval,smid);
else
umid = mean(U);
vmid = mean(V);
end
mag = sqrt(umid.^2 + vmid.^2);
si = vmid./mag;
co = umid./mag;

% no rotation
if (0)
	si(:) = 0; %mean(si);
	co(:) = 1; %mean(co);
end

% TODO actually the variance must be rotated as well cov(X*R) = R'*cov(X*R)
for idx=1:size(Uval,2)
	R = [co(idx) -si(idx); si(idx) co(idx)];
	UVr = [U(:,idx) V(:,idx)]*R;
	U(:,idx) = UVr(:,1);
	V(:,idx) = UVr(:,2);

	% TODO 2D rotation for UV
end

% normalise velocities

%plot(s,U.val);
SS = repmat(S,1,size(Uval,2));
%errorlines(S,Uval,Ustd);
U_ = bsxfun(@times,U,1./mag);
V_ = bsxfun(@times,V,1./mag);

% plot profiles
namedfigure(1,'Vertical profile of streamwise velocity');
clf();
plot(U_,SS);
legend('location','northwest',datestr(calib.t0,'dd/mm/yyyy'));

% only first three campaigns not to clutter the graphs
% these campaigns are min, max, centre
namedfigure(1000,'Vertical profile of streamwise velocity');
clf();
plot(U_(:,1:3),SS(:,1:3));
legend('location','northwest',datestr(calib.t0,'dd/mm/yyyy'));


namedfigure(2,'Vertical profile of spanwise velocity');
clf();
plot(V_,SS);

% regress slope to vertical profile of spanwise velocity
fdx = (S > 0.15 & S < 0.85);
c = ((S(fdx)-0.5) \ V_(fdx,:))';

namedfigure(3,'Strength of secondary current vs. depth');
clf
% fit a line
poly = PolyOLS(1);
poly.fit(calib.h0, 0.5*c);
[h0 sdx] = sort(calib.h0);
cc = colormap('lines');
for idx=1:length(sdx)
	plot(calib.h0(sdx(idx)),0.5*c(sdx(idx)),'ko','markerfacecolor',cc(idx,:));
	hold on
end
h = linspace(4,13,100);
plot(h,poly.predict(h'),'k');
text(5,7.5/100,sprintf('y=%1.2fx+%1.2f, R^2=%0.2f',poly.param*100,poly.R2))
legend('location','southeast',datestr(calib.t0(sdx),'dd/mm/yyyy'));
xlabel('z_s (m)');
ylabel({'$\frac{1}{2}\frac{1}{u}\frac{\partial v}{\partial \eta}|_{\eta = 0.5}$', '$\;\;(\%)$'},'interpreter','latex','rot',0);
percenttick('y');
[r p] = corr(calib.h0,c,'type','Pearson')

namedfigure(4,'Flow direction at mid depth');
clf
alpha = rad2deg(atan2(si,co));
%plot(h0,alpha,'o');
for idx=1:length(sdx)
	plot(calib.h0(sdx(idx)),alpha(sdx(idx)),'o','markerfacecolor',cc(idx,:));
	hold on
end
legend('location','northeast',datestr(calib.t0(sdx),'dd/mm/yyyy'));
ylabel('deviation $^\circ$','interpreter','latex');
xlabel('H (m)');

a0 = rad2deg(25);
d  = deg2rad(2);
R0 = [ cos(a0) sin(a0)
      -sin(a0) cos(a0) ];
Rd = [ cos(a0+d) sin(a0+d)
      -sin(a0+d) cos(a0+d) ];

R0 \ Rd
n

% Note: the plot was wrong in the previous version,
%       the depth cannot be scaled to unit depth without scaling the roughness length
namedfigure(5,'Vertical distribution of sw velocity');
clf
%plot(nanmean(U,2),S,'k');

mU = nanmean(U_,2);
sU = nanstd(U_,[],2);

Z = S(:,1)*calib.h0(:)';
mask = true(size(U));
n = round(0.84*size(mask,1));
mask(n:end,:) = false;
[us ln_z0 serr] = fit_log_profile(U,S*ones(1,size(U,2)),mask,cvec(calib.h0))
Z_  = repmat(max(Z,[],2),1,size(Z,2));
U__ = log_profile(rvec(us.val),rvec(ln_z0.val),S*ones(1,size(U,2)),rvec(max(calib.h0)))
mU_ = mean(U__);
sU_ = std(U_);

s = 1./mean(U__);
plot(bsxfun(@times,U,s),bsxfun(@times,Z,1));
%plot(bsxfun(@times,U,s),bsxfun(@times,Z,1),'k');
hold on
%plot(Z,bsxfun(@times,U,1),'k');
%plot(mU-sU,Z,'k--')
%plot(bsxfun(@times,U-Userr,s),Z,'k--')
%hold on
%plot(bsxfun(@times,U+Userr,s),Z,'k--')
hold on

%plot(bsxfun(@times,U__,s),Z_);


vline(1,'k:')
xlabel('$u/\bar u$','interpreter','latex');
%ylabel('$\eta$','interpreter','latex','rot',0);
ylabel('$z_s (m)$','interpreter','latex','rot',0);
%ylim([0 1.2])
%xlim([0.2 1.2]);
xlim([0.6 1.2]);
%ylim([0 1])
ylim([0 12.25])
legend('location','northwest',datestr(calib.t0,'dd/mm/yyyy'));


%figure(55);
%m=nanmean(U_,2);
%s = nanstd(U_,[],2);
%errorlines(S,m,m-s,m+s,'k'); 
%view([90,-90])


% cross term

%fdx = S > 0.85;
%V_ = V;
%U_ = U;
%for idx=1:calib.nc
%	gdx = isfinite(V_(:,idx));
%	fdx = ~gdx;
%	V_(fdx,idx) = interp1([S(gdx); 1],[V_(gdx,idx);0],S(fdx),'linear');
%	gdx = isfinite(U_(:,idx));
%	fdx = ~gdx;
%	U_(fdx,idx) = interp1([S(gdx); 1],[U_(gdx,idx);0],S(fdx),'linear');
%end

namedfigure(6,'UV');
plot(U.*V,SS);

figure(7);
clf
uv_avg = mean(U_.*V_)./mag.^2;
for idx=1:length(sdx)
        plot(calib.h0(sdx(idx)),uv_avg(sdx(idx)),'o','markerfacecolor',cc(idx,:));
        hold on                                                                 
end  

UVval_ = UVval;
figure(8);
clf();
plot(UVval_,S)
[v sdx]=sort(calib.h0);
m=mean(UVval_);
%plot(calib.h0(sdx),m(sdx)./mag);

if (exist('pflag','var') && pflag)
	pdfprint(3,'img/sanggau-dv-deta-vs-stage');
	pdfprint(3,'img/sanggau-dv-deta-vs-stage',2);
	pdfprint(4,'img/sanggau-flow-direction-vs-stage');
	pdfprint(5,'img/sanggau-u-div-ubar-vs-eta',2);
	pdfprint(5,'img/sanggau-u-div-ubar-vs-eta',3);
	pflag = 0;
end

