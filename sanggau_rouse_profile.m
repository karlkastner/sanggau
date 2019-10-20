% Wed 17 May 18:13:42 CEST 2017

%load /home/pia/phd/dat/kapuas/2013-12-09-sanggau-transect/mat/vadcp.mat
%vadcp = VADCP(vadcp);

% possible causes of low rouse number in fit:
% 	schmitd number not 1 (unlikely)
%	u_s underestimated (unlikely)
%	diameter of suspended sediment overestimated (maybe, but by a factor of 4?)
%	attenuation towards bed
%	outliers near the bed?
%	bend flow???

meta = sanggau_metadata();
ps = 3;

if(~exist('pflag','var'))
	pflag = 0;
	np = 1;
else
	np = 3;
end
fflag = pflag;

if (~exist('reload','var') || reload)
	cs_   = CrossSection();
	adcp_ = VADCP();
	Ro = [];
	us = [];
	for idx=1:length(meta.filename.discharge)
		clear cs	
		load(meta.filename.discharge{idx});
		cs_(idx)   = cs;
		clear vadcp
		load(meta.filename.vadcp{idx});
		adcp_(idx) = VADCP(vadcp);
	

	end
	% sort
	cs     = cs_;
	[t0 sdx] = sort(arrayfun(@(x) head(x.t0),cs));
	cs = cs(sdx);
	vadcp  = adcp_;
	[t0 sdx] = sort(arrayfun(@(x) head(x.t0),vadcp));
	vadcp = vadcp(sdx);
	H  = cvec(arrayfun(@(x) x.radius,cs));
	Hmax  = arrayfun(@(x) quantile(x.H,0.99),vadcp);
	Q  = cvec(arrayfun(@(x) x.discharge,cs));
	U  = cvec(arrayfun(@(x) x.Ut,cs));

	Ro = struct();
	us = struct();

	% fetch valus
	for idx=1:length(cs)
		% rouse number
		Ro.A(:,idx) = cs(idx).gridN.val.rouse;
		% shear velocity
		us.A(:,idx) = cs(idx).gridN.val.u_s;
	end
	% settling velocity
	ws.A     = 0.4*Ro.A.*us.A;

	N      = cs(1).gridN.cX1;
	inner = abs(N) < meta.nmax;
	%Ro     = real(Ro);
	%us     = real(us);
	reload = 0;
end

splitfigure([2 2],[1 1],fflag);
plot(N,Ro.A);
legend(datestr(arrayfun(@(x) x.t0,vadcp)));
ylabel('Rouse number')
xlabel('N (m)');

splitfigure([2 2],[1 2],fflag);
plot(N,us.A);
ylabel('Shear velocity (m/s)')
xlabel('N (m)');

splitfigure([2 2],[1 3],fflag);
plot(N,ws.A);
ylabel('Settling veloscity (m/s)')
xlabel('N (m)');

for idx=1:0 
	% length(cs)
	Ro = settling_velocity(0.6)/(0.4*nanmean(cs(idx).gridN.val.u_s))
	cs(idx).discharge.total
	subplot(2,2,2)
	plot(Q(idx),nanmedian(real(cs(idx).gridN.val.rouse)),'*')
	hold on
	set(gca,'colororderindex',get(gca,'colororderindex')-1);
	plot(Q(idx),Ro,'o');
	legend('fitted','theoretic')
end



% bin vertical profiles
nz=50;
%Hmax  =max(vadcp(1).H);
%Hmax  = arrayfun(@(x) max(x.H),vadcp);
gridS = Grid1([1/nz 1-1/nz],1/nz);
gridZ = Grid1(max(Hmax)*[1/nz 1-1/nz],max(Hmax)/nz);
valZ  = [];
valS  = [];
valZ_  = [];
valS_  = [];

%H = [];
% average profiles
% TODO inner region
for idx=1:length(vadcp)
	%cdx = ones(vadcp(idx).nt,1);
	cdx = 1;
	centre = cs(1).centre;
	centre(2) = centre(2)-1e7;
	vadcp(idx).xy2nts(cdx,centre,cs(1).dir,cs(1).dwidth,cs(1).T_max);

	S  = vadcp(idx).bin.S;
	Z  = vadcp(idx).bin.Z;

	bs_ = nanmean(vadcp(idx).beta,3);
	[bs_ C_ beta_] = backscatter_to_concentration2(Z(:,1),bs_);
%	bs_ = C_;

	msk = all(vadcp(idx).mask,3);
	msk = bsxfun(@and,msk,rvec(abs(vadcp(idx).ens.N))<meta.nmax);
	bs = bs_;
	bs(~msk | isinf(bs)) = NaN;

	gridZ.build_index(Z,'i1');
	gridS.build_index(S,'i1');
	%val(:,idx) = gridZ.binop('i1',@nanmedian,bs(:))
	valZ(:,idx) = gridZ.binop('i1',@nanmean,bs(:));
	valS(:,idx) = gridS.binop('i1',@nanmean,bs(:));
	
	% fit profiles for part of the cross section
	%m = 3;
	%N_ = linspace(-meta.nmax,meta.nmax,m+1);
	N_ = [-250,100,250];
	for jdx=1:length(N_)-1;
		%bs = nanmean(vadcp(idx).beta,3);
		bs  = bs_;
		msk = all(vadcp(idx).mask,3);
		msk = bsxfun(@and,msk,rvec(vadcp(idx).ens.N>N_(jdx) & vadcp(idx).ens.N>N_(jdx+1)));
		bs(~msk | isinf(bs)) = NaN;
		bs__ = C_;
		bs__(~msk | isinf(bs__)) = NaN;
		
		valZ_(:,idx,jdx) = gridZ.binop('i1',@nanmedian,bs(:));
		valS_(:,idx,jdx) = gridS.binop('i1',@nanmedian,bs(:));
		valZ__(:,idx,jdx) = gridZ.binop('i1',@nanmedian,bs__(:));
		valS__(:,idx,jdx) = gridS.binop('i1',@nanmedian,bs__(:));
	end
end
valS(1:5,:) = NaN;
valS(end-3:end,:) =NaN;
%H = rvec([cs.radius]);

rpS  = Rouse_Profile();
rpZ  = Rouse_Profile();
S    = cvec(gridS.cX1);
Z    = cvec(gridZ.cX1);
msk  = true(size(S));

rpS.fit(valS,S*rvec(H),H,isfinite(valS));
rpZ.fit(valZ,Z,Hmax,isfinite(valZ));
valSp = rpS.predict(S*rvec(H),H);
valZp = rpZ.predict(Z*ones(size(Hmax)),Hmax);

Ro_ =  [nanmean(Ro.A(inner,:))' rpS.c(2,:)' rpZ.c(2,:)'];
us_ =   nanmean(us.A(inner,:))';
ws_ =  [nanmean(ws.A(inner,:))', 0.4*(us_.*rpS.c(2,:)'), 0.4*(us_.*rpZ.c(2,:)')];
leg.Q = arrayfun(@(x) sprintf('Q=%g m^3/s',x),round(abs(Q)),'uniformoutput',false);
leg.U = arrayfun(@(x) sprintf('%g',x),U,'uniformoutput',false);

splitfigure([2 2],[2 1],fflag);
cla();
plot(valS,S*rvec(H));
hold on;
set(gca,'colororderindex',1);
plot(valSp,S*rvec(H),'--');
%legend(datestr(arrayfun(@(x) x.t0,vadcp)))
h=legend(leg.Q{:});
xlabel('Backscatter');
ylabel('Z (m)')
%ylabel('S (z/H)');
%ylabel('S');
%view(90,-90);

splitfigure([2 2],[2 2],fflag);
cla();
% TODO average along z
plot(valZ,Z*ones(1,5));
%plot(Z,valZ);
hold on
set(gca,'colororderindex',1);
%plot(Z,valZp,'--');
%view(90,-90);
valZp(isnan(valZ)) = NaN; % invalidate above surface
plot(valZp,Z*ones(1,5),'--');
xlabel('Basckatter');
ylabel('z (m)');

splitfigure([2 2],[2 3],fflag);
cla
for idx=1:size(valS_,3)
	plot(valS_(:,:,idx),S*rvec(H));
	hold on
	set(gca,'colororderindex',1);
	plot(valS__(:,:,idx),S*rvec(H),'--');
	set(gca,'colororderindex',1);
'hink'
pause
end
xlabel('Backscatter');
ylabel('Z (m)')

splitfigure([2 2],[2 4],fflag);
cla
for idx=1:size(valS_,3)
	set(gca,'colororderindex',1);
	plot(valZ_(:,:,idx),Z);
	hold on
	plot(valZ__(:,:,idx),Z,'--');
	set(gca,'colororderindex',1);
end
xlabel('Backscatter');
ylabel('Z (m)')
h=legend(leg.Q{:});

splitfigure([2 2],[3 3],fflag);
cla();
plot(H,Ro_(:,1:np),'o');
% fit and predict
pls = PolyOLS(1);
pls.fit(H,Ro_(:,1:np));
hold on;
set(gca,'colororderindex',1);
H100 = linspace(0,max(H))';
Ro__ = pls.predict(H100);
plot(H100,Ro__(:,1:np));
ylabel('Rouse number');
xlabel('Discharge (10^3 m^3/s)');
xlim([0 1.1*max(H)]);

splitfigure([2 2],[3 4],fflag);
cla();
plot(H,ws_(:,1:np),'*');	% use Q?
pls = PolyOLS(1);
pls.fit(H,ws_(:,1:np));
hold on;
set(gca,'colororderindex',1);
ws__ = pls.predict(H100);
plot(H100,ws__(:,1:np));
ylim([0 1.1*max(ws_(:))]);
xlabel('Hydraulic radius (m)');
ylabel('Settling velocity w_s (m/s)');
xlim([0 1.1*max(H)]);
%h=legend(leg.Q{:});

if (pflag)
	pdfprint(11,'img/sanggau-rouse-number-vs-n',ps);
	pdfprint(13,'img/sanggau-rouse-number-vs-n',ps);
	pdfprint(23,'img/sanggau-rouse-number-vs-q',ps);
	pdfprint(24,'img/sanggau-settling-velocity-number-vs-q',ps);
	pflag = 0;
end

