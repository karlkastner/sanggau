% Do 14. Mai 11:43:31 CEST 2015
% Karl Kastner, Berlin

if (~exist('reload','var') || reload)
	meta = sanggau_metadata();


	load([ROOTFOLDER,'/src/water-level/mat/water-level.mat']);
	level.val  = Kr(nid.sanggau_merged).depth;
	level.time = Kr(nid.sanggau_merged).time;

	calib = load_vadcp_discharge(meta.filename.discharge,level);
	
	n = length(calib.cs_);
	clear adcp_;
	for idx=1:n
		load(meta.filename.vadcp{idx})
		vadcp = VADCP(vadcp);
		cs = calib.cs_(idx);
		centre = cs.centre;
		centre(2) = centre(2)-1e7;
		dir = cs.dir;
		% rotate, such that [dzb/ds; dzb/dn] = [1;0]
		% the optimal rotation angle is difficult to determine from the data,
		% due to bed forms, different water level and the low number of cross sections
		% because this is a bend, the linear transformation is also not exact,
		% the coordinates should be transformed along an arc
		%s = -(0.0370+0.0163+0.0105+0.0102+0.0102-5.3343e-04+0.0017+0.0050+0.0023); c = sqrt(1-s.^2);
		%s = -20/300;
		% iterative optimal angle (not converging)
		a = -(0.0370+0.0167+0.0100+0.0089+0.0081+0.0026+0.0030+0.0026);
		%a = -0.045;
		a = -0.06;
		s = sin(a);
		c = sqrt(1-s^2);
		R = [c s; -s c];
		dir = R*dir;
		vadcp.xy2nts(cs.cdx,centre,dir,cs.dwidth,cs.T_max);
		% TODO quick fix
		calib.dis(idx).ens.level4  = calib.l0(idx);
		adcp_(idx) = vadcp;
		clear  vadcp;
	end
	% bring into same order
	f=@(x) x.time(1);
	t=arrayfun(f,adcp_)';
	for idx=1:n
		[mv mdx] = min(abs(calib.t0(idx)-t))
		calib.dis(idx).adcp = adcp_(mdx);
	end
% concatenate data
%joint = struct('X',[],'Y',[],'H',[]);
dat  = struct('X',[],'Y',[],'H',[]);
dat(n+1).name = 'combined';
for idx=1:n
	dat(idx).name = datestr(calib.t0(idx),'dd/mm/yyyy');
	H   = calib.dis(idx).adcp.H4(:) - calib.dis(idx).ens.level4(:);
	X   = calib.dis(idx).adcp.T4(:);
	Y   = calib.dis(idx).adcp.N4(:);
	fdx = isfinite(H.*X.*Y);
	dat(idx).H = H(fdx);
	dat(idx).X = X(fdx);
	dat(idx).Y = Y(fdx);
	dat(n+1).H = [dat(n+1).H; dat(idx).H];
	dat(n+1).X = [dat(n+1).X; dat(idx).X];
	dat(n+1).Y = [dat(n+1).Y; dat(idx).Y];
end

	reload = 0;
end

% qunatile of plot range
p = 0.99;

% TODO offset at individual campaigns


%namedfigure(1,'Combined bathymetry')
%clf();
%T=delaunay(joint.X,joint.Y);
%patch(joint.X(T)',joint.Y(T)',joint.H(T)','edgecolor','none');
%d = max(diff(quantile(X,[0.01,0.99])),diff(quantile(Y,[0.01,0.99])));
%r = sqrt(quantile( (X-median(X)).^2 + (Y-median(Y)).^2, 0.99));
%xlim(median(X)+r*[-1 1]);
%ylim(median(Y)+r*[-1 1]);
%q = quantile(joint.X,[p,1-p]);
%xlim(q);
%colormap jet
%colorbar
%	caxis([10 20])

% plot 2d bathymetry
for idx=1:n+1
%	H = calib.dis(idx).adcp.H4(:) - calib.dis(idx).ens.level4(:);
%	X = calib.dis(idx).adcp.T4(:);
%	Y = calib.dis(idx).adcp.N4(:);
%	fdx = isfinite(H.*X.*Y);
%	H = H(fdx);
%	X = X(fdx);
%	Y = Y(fdx);
%	joint.H = [joint.H; H];
%	joint.X = [joint.X; X];
%	joint.Y = [joint.Y; Y];

	namedfigure(idx,['Bathymetry ',dat(idx).name]);
	clf();
	%subplot(2,2,1)
	T=delaunay(dat(idx).X,dat(idx).Y);
	patch(dat(idx).X(T)',dat(idx).Y(T)',dat(idx).H(T)','edgecolor','none');
	%d = max(diff(quantile(X,[0.01,0.99])),diff(quantile(Y,[0.01,0.99])));
	%mx = median(dat(idx).X);
	%r = sqrt(quantile( (X-median(X)).^2 + (Y-median(Y)).^2, 0.99));
	r = quantile(dat(n+1).X,[1-p,p]);
	r = max(abs(r));
	xlim(r*[-1 1]);
	%xlim(median(X)+r*[-1 1]);
	%ylim(median(Y)+r*[-1 1]);
	%xlim(median(X)+r*[-1 1]);
	colormap jet
	colorbar
end

% bed slope at individual campaign
dn = 1;
N = [-315:dn:315];
P = [];
% regress local bed slope
for idx=1:n+1
	%H = calib.dis(idx).adcp.H4(:) - calib.dis(idx).ens.level4(:);
	%X = calib.dis(idx).adcp.T4(:);
	%Y = calib.dis(idx).adcp.N4(:);
	%fdx = isfinite(H.*X.*Y);
	%H = H(fdx);
	%X = X(fdx);
	%Y = Y(fdx);
	param = [];
	for jdx=1:length(N)-1
		fdx  = (dat(idx).Y>N(jdx)) & (dat(idx).Y<N(jdx+1));
		poly = PolyOLS(1);
		poly.fit(dat(idx).X(fdx),dat(idx).H(fdx));
		param(jdx,:) = poly.param;
	end
	P(:,:,idx) = param;
%	figure(10);
%	subplot(2,2,1)
%	plot(N(1:end-1),param(:,1));
%	hold on
%	subplot(2,2,2)
%	plot(N(1:end-1),param(:,2));
%	hold on
end
figure(10);
clf
subplot(2,2,1)
plot(N(1:end-1),squeeze(P(:,1,:)));
%hold on
subplot(2,2,2)
plot(N(1:end-1),squeeze(P(:,2,:)));
%hold on
subplot(2,2,3);
dh_dn = cdiff(P(:,1,end))./dn;
dh_ds = P(:,2,end);
h = hypot(dh_dn,dh_ds);
s = dh_dn./h;
c = dh_ds./h;
d = [s c];
mean(c)
mc = median(c)
ms = sqrt(1-mc.^2);
s_ = sqrt(1-mc.^2);
ma  = pi/2-atan2(ms,mc)
a = pi/2-atan2(s_,c);
[ma sma lma uma] = median_man(a)
%d = [dh_dn./h, dh_ds./h];
%s = d(:,1)./
%plot(N(1:end-1),dh_dn);
plot(N(1:end-1),d);
%calib.cs_(1).gridN.cX,

subplot(2,2,4)
m = [P(:,2,end) median(P(:,2,:),3) mean(P(:,2,:),3)];
mean(m)
median(m)
plot(N(1:end-1),m);
%P(:,2,end));
%hold on
%plot(N(1:end-1),median(P(:,2,:),3));
%plot(N(1:end-1),mean(P(:,2,:),3));

method = {'nearest','linear','natural'}
%method = {'natural'};
%method = {'linear'};

for idx=1:n+1

X = dat(idx).X;
Y = dat(idx).Y;
H = dat(idx).H;

%namedfigure(length(calib.dis)+2,'Longitudinal profiles combined data');
for jdx=1:length(method)
	interp(jdx).z = TriScatteredInterp(X,Y,double(H),method{jdx});

	interp(jdx).x = TriScatteredInterp(X,Y,double(X),method{jdx});
	interp(jdx).y = TriScatteredInterp(X,Y,double(Y),method{jdx});

end

t = cvec((median(X)-r:median(X)+r));
nn = median(Y) + 100*(-2:2);
nn = 100*(-2:2);
h = [];
% TODO no magic numbers
nf = 25;
dmax = 50;
for kdx=1:length(nn)
	n_ =  nn(kdx)*ones(size(t));
	for jdx=1:length(method)
		h_ = interp(jdx).z(t,n_);
		% invalidate values out of range
		%dx = interp(1).x(t,n_)-t;
		%dy = interp(1).y(t,n_)-n_;
		[id dis]=knnsearch([X Y],[t,n_]);
		%$ = hypot(dx,dy);
		fdx = dis > dmax;
		h_(fdx) = NaN;
		h(:,kdx,jdx)  = h_;
		% TODO, invalidate where distance too large
		hf(:,kdx,jdx) = meanfilt1(h(:,kdx,jdx),nf);
	end
end
%h = hf;

figure(1000+idx);
clf
plot(t,h(:,:,1))
hold on
set(gca, 'ColorOrderIndex', 1)
plot(t,hf(:,:,1),'--')
%plot(t,h2,'--')
%set(gca, 'ColorOrderIndex', 1)
%lot(t,hf(:,:,1),'--')
%et(gca, 'ColorOrderIndex', 1)
%lot(t,hf(:,:,3),'.')
%egend(num2str(cvec(nn)));

end % for idx

% extremum detection
%{
d = sign(diff(h(:,:,3)));
n_ = 10;
hold on
f = {}
dnd = [];
for idx=1:size(d,2)
	dnd(:,idx) = d(1:end-n_,idx)-d(n_+1:end,idx);
	f{idx} = find(dnd(:,idx));
	plot(t(f{idx}),h(f{idx},idx,3),'ok');
end
%}

nn = (-250:1:250);
h = [];
nf = 20;
dx = [];
dy = [];
for idx=1:length(nn)
	n_ =  nn(idx)*ones(size(t));
	x0 = interp(1).x(t,n_);
	x1 = interp(2).x(t,n_);
	y0 = interp(1).y(t,n_);
	y1 = interp(2).y(t,n_);
%	h(:,idx,2)  = interp1(t,n);
	h(:,idx) = interp(3).z(t,n_);
	h(:,idx) = meanfilt1(h(:,idx),nf);

	% for distance to nearest point
	dx(:,idx) = x1-x0;
	dy(:,idx) = y1-y0;
end
ds2 = dx.^2 + dy.^2;
ds = sqrt(ds2);
figure(100);
subplot(2,2,1)
imagesc(ds);
colorbar
subplot(2,2,2)
plot([nanmedian(abs(dx)); nanmedian(abs(dy)); nanmedian(ds); nanmax(ds)]')
subplot(2,2,3);
plot(sum(isfinite(ds)))
median(sum(isfinite(ds)))

inner=true(size(t)); %t>=-30&t<=70;
inner=t>=-100&t<=200;
% subtract linear part
A = vander_1d(cvec(t(inner)),1);
p = A\h(inner,:);
d = h(inner,:) - A*p;
	%d(:,idx)=bsxfun(@minus,h(inner,:,3),nanmedian(h(inner,:,3)));
figure(11)
clf();
dt = 1;
h_ = h(inner,:);
for idx=1:size(d,2)
	gdx=detect_extreme(d(:,idx),0.2);
	% discard first
	gdx = gdx(2:end);
	if (mod(idx,50) == 1)
	fdx=find(inner);
	plot(t(inner),d(:,idx));
	hold on;
	plot(t(fdx(gdx)),d(gdx,idx),'ko');
	grid on;
%	pause
	end
%	dune.h{idx} = abs(0.5*h_(gdx(1:end-2))-h_(gdx(2:end-1))+0.5*h_(gdx(3:end)));
	dune.level{idx} = d(gdx,idx);
	dune.id{idx} = gdx;
	%dune.h{idx} = abs(d(gdx(1:end-2),idx)-d(gdx(2:end-1),idx));
	% note that this is attenuated somewhat by the filter
	dune.h{idx} = abs(0.5*d(gdx(1:end-2),idx)-d(gdx(2:end-1),idx)+0.5*d(gdx(3:end),idx));
	dune.l{idx} = dt*(gdx(3:end)-gdx(1:end-2));
end

h = cellfun(@mean,dune.h);
sh = cellfun(@serr,dune.h);
l = cellfun(@mean,dune.l);
sl = cellfun(@serr,dune.l);
h = cvec(h);
l = cvec(l);


% correct height of dunes for attenuation by smoothing
a = abs(frequency_response_boxcar(nf,1./l))
h = h./a

figure(1);
subplot(2,2,1)
%plot(nn,h)
errorbar(nn,h,sh,'o')
subplot(2,2,2)
errorbar(nn,l,sl,'o')
subplot(2,2,3)
plot(a)

d90    = 5e-4;

[ks ln_z0] = roughness_length_rijn(h,l,d90)

nanmedian(ks.bedform(nn<100))

h = nanmedian(h(nn<100))
l = nanmedian(l(nn<100))

[ks ln_z0] = roughness_length_rijn(h,l,d90)

height = h;
length = l;
save('mat/sanggau-bedform-data.mat','height','length','ln_z0');
clear length

if (exist('pflag','var') && pflag)
	print(length(calib.dis)+1,'img/sanggau-2d-bathymetry.png','-dpng')
	pflag = 0;
end


