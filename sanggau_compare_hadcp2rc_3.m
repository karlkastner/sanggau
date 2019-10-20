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

label = {''};
field_C = { {'A'},
            {'A', 'AV'},
            {'A', 'AV','AVL'},
            {'A', 'AV', 'AVL', 'AVL2'}
            {'A', 'AV', 'AVL', 'AVL2','Adl_dt'}
	};
%lab = {'A','A+AV','A+AV+AVL','A+AV+AVL^2','A+AV+AVL+AVL^2'
lab = {'A','+AV','+AVL','+AVL^2','+A dL/dt'}
		
nt      = length(Qr);

%Qh_   = NaN(nt,length(field),size(V,2));
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
AVL   = A.*V(:,jdx).*l;  
% area velocity level squared    
AVL2   = A.*V(:,jdx).*l.^2;
% cubed
AVL3   = A.*V(:,jdx).*l.^3;
% cubed
AVl_4   = A.*V(:,jdx).*l.^4;

% deployment periods
tbreak = meta.h_tbreak;
% regress coefficients for each deployment period
c = [];
Qh_ = NaN(nt,length(field_C));;
for idx=1:length(tbreak)-1
	fdx_i = (t>=tbreak(idx)) & (t<tbreak(idx+1)) & fdx;
	if (sum(fdx_i)>0)
		for kdx=1:length(field_C)
			field = field_C{kdx};
			% set up regression matrix
			X = [];
			for jdx=1:length(field)
				X(:,jdx) = eval(field{jdx});
			end % for jdx
			% regression coefficients for the current deployment period
			c = X(fdx_i,:) \ Qr(fdx_i);
			% predicted discharge during current deployment period (column wise)
		%	q = bsxfun(@times,X(fdx_i,:),c.');
		%	Qh_(fdx_i,size(q,2),kdx)  = q;
			Qh_(fdx_i,kdx) = X(fdx_i,:)*c;
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
Qh = Qh_;
Q1(~fdx,:,:) = NaN;

% prediction residual
res  = bsxfun(@minus,Qh,Qr);

rmse = nanrms(res);

namedfigure(2,'HADCP vs Rating Curve Discharge over time');
clf
plot(t,[Qr Qh]);
%plot(t,[U U0 U1]);

namedfigure(3,'Residual over time');
plotyy(t,Qr,t,res);
%plot(t,res);
vline(tbreak(2:end),'r')

namedfigure(4,'Residual distribution');
clf
if (0)
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
end

%subplot(2,2,2);
%plot(Qr,medfilt1(Qh,51))

subplot(2,2,3);
plot(dl_dt,res);

namedfigure(5,'Polynomial discharg estimate');
clf();
res_ = bsxfun(@minus,Qh,Qr);
w = kaiserwin(1:3*144).';
res_f = wmeanfilt(w,res_(:,1:3));
plot(Qr/1e3,res_f/1e3);
hold on

for idx=-15:5:15
	plot([0 10],[0  idx*(1e-1)],'k-')
end
xlim([0 10]);
%ylim([-15 15])
ylim([-1.5 1.5]);
vline(min(Qr)/1e3,'k--');
ax = gca;
addy(gca)
ytick(-15:5:15);
ylim([-15 15]);
yticklabel(sprintf('%d%%\n',(-15:5:15)'))
xtick([]);
xlabel(ax(1),'Q_{RC} (10^3 m^3/s)');
ylabel(ax(1),'Q_{HADC} - Q_{RC} (10^3 m^3/s)');
legend(ax(1),'location','southwest','A','A+AV','A+AV+AVL');

%box off

namedfigure(6,'RMSD of discharge');
clf
%s2_ = struct2array(s2);
%s2_ = 100*double(s2_(2:end))/s2.total;
ax = axis;
hold on
%bar(s2_);
rmse = nanrms(res_);
bar(rmse);
%bar(1:6,s2_(1:6))
%bar(7:length(s2_),[s2_(7:end)],'r');
%fn = fieldnames(s2);
fn = label;
xticklabel(fn{1:end});
%ylabel('VAR(Q_x)/VAR(Q_{RC}) (%)');
ylabel('RMS(Q_{HADCP} - Q_{RC})')
%text((1:length(s2_))-0.25, ...
%		double(s2_(1:end))+3,num2str(cvec(s2_(1:end)),'%2.2f'));
%xlim([0.1 length(s2_)-0.1+1])
xtick(1:5);
xticklabel(lab{:})
ax = gca;
yl = ylim;
ax(2) = addx(ax);
xtick([]);
Qrmax = max(Qr);
ylim(ax(2),100*yl/Qrmax);
ylabel('% of max(Q)');

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
	pdfprint(5,'img/rms-qhadcp-qrc-vs-qrc.pdf',printscale);
	pdfprint(6,'img/rms-hadcp-methods.pdf',printscale);
	pflag = 0;
end


