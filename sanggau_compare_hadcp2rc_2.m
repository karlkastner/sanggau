% Tue 14 Mar 11:09:36 CET 2017
%TODO calibrate each beam individually
%TODO determine variance along beam


meta   = sanggau_metadata();
if (~exist('reload','var') || reload)
	[hadcp hlevel offset] = sanggau_load_hadcp(meta,0,1800);
	rc     = load('../discharge/mat/sanggau-stage-discharge.mat');
	reload = 0;
end
if (~exist('pflag','var'))
	pflag = 0;
end
mode = 0;
printscale = 2.5;

% numbers of bins to use
% (beyond beam 60 bed interference occured)
m      = 60;
% maximum gap length for fixing data series
dt_max = 3/24;
% TODO rotate to cs
v = -hadcp.velocity.earth(:,:,1);
v0 = v;

% fix short gaps in water level and rating curve discharge
t     = rc.time;
level = fixnan(t,rc.level,dt_max,'linear');
A     = rc.area_exact_f(level);
Qr    = rc.qe_lin;
Qr    = fixnan(t,Qr,dt_max,'linear');

% fix VADCP along beam
nbin = size(v,1);
dnmax = 3;
v = fixnan(1:nbin,v,dnmax,'linear','extrap');

figure(1000);
V = nanmean(v(1:m,:));
d=0.5*(V(2:end-1)-(V(1:end-2)+V(3:end))/2);
d=d(isfinite(d));
qstd(d),
q=quantile(d,[0.01 0.99]); kernel1d(linspace(q(1),q(2)),d); vline(quantile(d,normcdf([-2:2])))

% interpolate VADCP to level time slots
v     = interp1_limited(hadcp.time,v',t,dt_max,'linear')';
p     = interp1_limited(hadcp.time,hadcp.pressure_bar,t,dt_max,'linear');

% bulk average velocity
if (0 == mode)
	V  =  nanmean(v(1:m,:))';
	V0 = V;
else
	V  = v(1:m,:).';
	V0 = V;
end

% set velocity data to NaN, when HADCP was not submerged
fdx = p < 0.01; % 10cm
V(fdx,:) = NaN;
h_tinvalid = meta.h_tinvalid;
for idx=1:size(h_tinvalid,1)
	fdx = (t>=h_tinvalid(idx,1)) & (t<=h_tinvalid(idx,2));
	V(fdx,:) = NaN;
end

if (0)
	namedfigure(1,'VADCP');
	clf();
	plot(t,[V0 V]);
	legend('All','Filtered');
end

%
% HADCP calibration
%

field = {'A','AV','AVl_1','AVl_2', 'AVl_3', 'AVl_4','O', 'Adl_dt'}; label = field;
%field = {'AV','AVl_1','AVl_2', 'AVl_3', 'AVl_4','A','O', 'Adl_dt'}; label = field;
%field = {'AV','AVl_1','AVl_2', 'A', 'O', 'Adl_dt'};
%label = {'AV^{ }','AVL^{ }','AVL^2','A','c','A^{ }dL/dt','Residual^{ }'};
%field = {'AV','AVl_1', 'A', 'Adl_dt'};
%field = {'O','AV','AVl_1','AVl_2', 'AVl_3', 'AVl_4','A', 'Adl_dt'};
%label = {'','AV^{ }','AVL^{ }','AVL^2','AVL^3','AVL^4','A^{ }','A^{ }dL/dt','Residual^{ }'};

nt      = length(Qr);

Qh_   = NaN(nt,length(field),size(V,2));
Q1    = NaN(nt,4,size(V,2));
% scales = inverse velocity profile
fi     = NaN(nt,5,size(V,2));

% select valid samples
dt = t(2)-t(1);
level_f = meanfilt1(level,round(1./dt));
dl_dt   = cdiff(level_f)./cdiff(t);
fdx    = isfinite(Qr) & all(isfinite(V),2) & isfinite(dl_dt) & isfinite(level);
% area dl/dt
Adl_dt  = A.*dl_dt;

% constant
O	= ones(nt,1);
O(fdx)  = O(fdx)/norm(O(fdx));
% cross section area
A       = A;

l = level;
l = l-midrange(l(fdx));
l = l/max(l(fdx));

% for each ADCP bin (if not lumped)
for jdx=1:size(V,2)
disp(jdx)

% area HADCP velocity
AV      = A.*V(:,jdx);
% area velocity level
AVl_1   = A.*V(:,jdx).*l;  
% area velocity level squared    
AVl_2   = A.*V(:,jdx).*l.^2;
% cubed
AVl_3   = A.*V(:,jdx).*l.^3;
% cubed
AVl_4   = A.*V(:,jdx).*l.^4;

% set up regression matrix
X     = NaN(nt,length(field));
for idx=1:length(field)
	% this precludes a QR factorisation
	X(:,idx)   = eval(field{idx});
end
% orthogonalise columns with respect to each other
% TODO should the columns better be orthogonalised per deployment period?
Xq      = NaN(nt,length(field));
[Xq(fdx,:) Xr] = qr(X(fdx,:),0);

% deployment periods
tbreak = meta.h_tbreak;
% regress coefficients for each deployment period
c = [];
for idx=1:length(tbreak)-1
	fdx_i = (t>=tbreak(idx)) & (t<tbreak(idx+1)) & fdx;
	if (sum(fdx_i)>0)
		% regression coefficients for the current deployment period
		c(:,idx) = Xq(fdx_i,:) \ Qr(fdx_i);
		% predicted discharge during current deployment period (column wise)
		Qh_(fdx_i,:,jdx)  = bsxfun(@times,Xq(fdx_i,:),c(:,idx).');
		% velocity scales
%		fi(fdx_i,:,jdx)    = bsxfun(@times,Qh_(fdx_i,2:6,jdx),1./V(fdx_i,jdx));
		for kdx=1:size(Xq,2);
			id = 1:kdx;
			Q1(fdx_i,kdx,jdx) = Xq(fdx_i,id)*c(id,idx);
		end % for kdx
	end % if sum(fdx) > 0
end % for idx

end % for jdx

% harmonic mean estimate
switch (mode)
case {0}
	% nothing to do, velocity was averaged along range beforehand
case {1}
	% arithmetic mean along beam
	Qh_ = mean(Qh_,3);
	Q1  = mean(Q1,3);  
case {2}
	% harmonic mean
	% this is more complicated, first, the scales have to be averaged harmonically and the velocity arithmetically,
	% non-velocity terms are averaged arithmetically

	% sum up parts of the velocity scale (columns)
	fi_ = fi;
	fi  = sum(fi,2);
	% harmonic mean of velocity profile along profiling range = arithmetic mean of scales
	fi = mean(fi,3);

	Qh  = fi.*mean(V,2) + sum(mean(Qh_(:,[1 7 8],:),3),2);

	mean(sum(sum(shiftdim(bsxfun(@times,shiftdim(fi_(fdx,:,:),2),V(fdx,:)'),1),3),2))
%	mean(sum(bsxfun(@times,f_(fdx,:,:),V(fdx,:)),3))
	mean(sum(fi(fdx).*mean(V(fdx,:),2)))
	Qh_ = mean(Qh_,3);
	mean(sum(Qh_(fdx,2:6),2))

	figure(100)
%	plotyy(t,fi.*mean(V,2),t,sum(mean(shiftdim(bsxfun(@times,shiftdim(fi_(:,:,:),2),V(:,:)'),1),3),2))
	plot(t,[fi.*mean(V,2),sum(mean(shiftdim(bsxfun(@times,shiftdim(fi_(:,:,:),2),V(:,:)'),1),3),2)])
	
	Qh_ = sum(Qh_,3);
	% this is not distinguished
	Q1  = mean(Q1,3);  
end

Q1(~fdx,:,:) = NaN;

% prediction of total discharge
Qh   = sum(Qh_,2);
res  = Qh - Qr;
res1 = bsxfun(@minus,Q1,Qr);
%res1 = meanfilt1(res1,round(0.25/dt));

% total variance of the rating curve discharge
s2        = struct();
s2.total  = var(Qr(fdx));

% variance of columns
s2res = struct();
% incremental
for idx=1:length(field)
	res_i            = Qr - sum(Qh_(:,1:idx),2);
	s2res_i          = rms(res_i(fdx)).^2;
	s2explained_incremental(idx) = s2.total - s2res_i;
	if (1==idx)
		s2.(field{idx}) = s2explained_incremental(idx);
	else
		s2.(field{idx}) = s2explained_incremental(idx)-s2explained_incremental(idx-1);
	end
	
%	s2.(field{idx}) = var(Qh_(fdx,idx)); % rms^2
end

% residual (error or unexplained) variance (var = rms, as mean is zero)
s2.res   = rms(res(fdx)).^2;
label{end+1} = 'Residual';

% short term variation (high frequent noise)
if (0)
dt = t(2)-t(1);
n = round(1./dt);
%wres    = meanfilt1(res,n);
wres     = medfilt1(double(res),n);
s2.noise = nanrms(res-wres).^2;
end

namedfigure(2,'HADCP vs Rating Curve Discharge');
clf
plot(t,[Qr Qh]);
%plot(t,[U U0 U1]);

namedfigure(3,'Q_HADCP - Q_RC vs Q');
plotyy(t,Qr,t,res);
%plot(t,res);
vline(tbreak(2:end),'r')

namedfigure(4,'Q_HADCP - Q_RC');
clf
subplot(2,2,1);
x = linspace(-max(abs(res)),max(abs(res)));
s = 10;
y = kernel1d(x,res(isfinite(res)),s);
plot(x,y);
xlabel('Q_{HADCP}-Q_{RC} (m^3/s)')
q = quantile(res,normcdf(-2:2));
lim = max(quantile(res,normcdf([-3 3])));
vline(q);
%xlim([-lim lim]);
%semilogy(x,y);

subplot(2,2,2);
plot(Qr,medfilt1(Qh,51))

subplot(2,2,3);
plot(dl_dt,res);

figure(5);
clf();
Qh__ = [Qh_(:,1), sum(Qh_(:,[1 4]),2) sum(Qh_(:,1:4),2)];
res_ = bsxfun(@minus,Qh__,Qr);
plot(Qr/1e3,res_/1e3);
xlabel('Q_{RC} (10^3 m^3/s)');
ylabel('Q_{HADC} - Q_{RC} (10^3 m^3/s)');
xlim([0 10]);
vline(min(Qr)/1e3,'k--');
%legend('location','southwest','AV','AV+AVL','AV+AVL+AVL^2+A','AV+AVL+AVL^2+A+AVL^3+AVL^4');
legend('location','southwest','AV','AV+A','AV+A+AVL+AVL^2');


hold on
plot([0 10],[0  15*(1e-1)],'k-')
plot([0 10],[0 -15*(1e-1)],'k-')
plot([0 10],[0  10*(1e-1)],'k-')
plot([0 10],[0 -10*(1e-1)],'k-')
plot([0 10],[0  5*(1e-1)],'k-')
plot([0 10],[0 -5*(1e-1)],'k-')
plot([0 10],[0 0],'k-')
ylim([-2 1.5]);
addy(gca)
ytick(-15:5:15); ylim([-20 15]); yticklabel(sprintf('%d%%\n',(-15:5:15)'))
xtick([]);
%box off

namedfigure(6,'Variance of discharge');
clf
s2_ = struct2array(s2);
s2_ = 100*double(s2_(2:end))/s2.total;
ax = axis;
hold on
bar(s2_);
%bar(1:6,s2_(1:6))
%bar(7:length(s2_),[s2_(7:end)],'r');
%fn = fieldnames(s2);
fn = label;
xticklabel(fn{1:end});
ylabel('VAR(Q_x)/VAR(Q_{RC}) (%)');
text((1:length(s2_))-0.25, ...
		double(s2_(1:end))+3,num2str(cvec(s2_(1:end)),'%2.2f'));
xlim([0.1 length(s2_)-0.1+1])

% note: magnitude of regression coefficients are proportinal to explained variance, as columns are normalised

% TODO operation readiness

%spearman_to_pearson(corr(dl_dt(fdx),res(fdx),'type','spearman'));
%nanr2(Qr,Qh)
%for idx=1:size(Qh,2)
%	nanr2(Qr,Qh(:,idx))
%end


%plot(t,[U1])
%ylim([0 1.3]);

% plot monthly intervals
if (0)
dt = t(2)-t(1);
n = round(30/dt);
N = length(t);
for idx=1:ceil(N/n)
	figure(100+idx);
	fdx = (idx-1)*n+1:min(N,idx*n);
	plot(t(fdx),V(fdx));
	weektick();
end
end


if (pflag)
	pdfprint(6,'img/var-q-composition.pdf',printscale);
	pdfprint(5,'img/qhadcp-qrc-vs-qrc.pdf',printscale);
	pflag = 0;
end


