% Mon Nov  3 09:36:44 CET 2014
% Karl Kastner, Berlin

%
% TODO chezy's and mannigs are just special cases of the powerrc
%

% TODO: account for instationarity
% Note: neither the polynomial not the exponential rating curve is forced
%       through zero, as it is unclear what the zero level is due to the non
%       rectangular shape of the cross section
% TODO: what is the difference between jacknifing prediction and the parameter estimate?
% Note: interestingly the jackknife removes the quadratic term
%       this indicates that the quadratic variation in the samples is mostly due to uncertainty and not due to physics
%       (note that this is not due to the low number of samples, the jacknife will not remove the quadratic term in general, as can be seen by sampling four points of a quadratic without noise in which case jacknife estimate is ident to the true quadratic
%	matrix is ident and the quadratic term remains, 

function ret = sanggau_stage_discharge()

meta = sanggau_metadata();
base = [ROOTFOLDER,'/src/discharge/mat/sanggau-stage-discharge'];

% order of polynomial rating curve
porder = 2;

% degree of fit to cs area and perimeter
% for 6 the perimeter decreases at the end (not realistic)
% for 5 and 4 the linear term is negative (unrealistic)
np = 2;

% TODO constrained fit to enforce dA/dh and dP/dh > 0
%	i.e. no extremum or saddle point (dA/dh = 0) in the stage range

mcmc = load('mat/sanggau-mcmc.mat');
% preload water level
% preload water level
load([ROOTFOLDER,'/src/water-level/mat/water-level.mat']);
level.val  = Kr(nid.sanggau_merged).depth;
level.time = Kr(nid.sanggau_merged).time;
% preload vadcp calibrarion discharge
calib      = load_vadcp_discharge(meta.filename.discharge,level);
calib.dis  = [];

% level whith respect to sea level
calib.zs0  = calib.zs0  - meta.z_offset;
level.val  = level.val  - meta.z_offset;
z_b        = smooth(calib.zb.median,meta.ns)-meta.z_offset;
l0         = double(calib.zs0);

% re-evaluate cross section geometry after smoothing
% area and perimeter for each calibration
% note: for consistence it is necessary to take the average bed level
for idx=1:calib.nc
	%b = smooth(calib.zb.A(:,idx),ns);
	b = z_b;
	a0(idx,1) = csarea(l0(idx),b,calib.dw);
	p0(idx,1) = csperimeter(l0(idx),b,calib.dw);
	r0(idx,1) = csradius(l0(idx),b,calib.dw);
end

% area function
area_exact_f = @(l) csarea(l,z_b,calib.dw);

% wetted perimeter function
perimeter_exact_f = @(l) csperimeter(l,z_b,calib.dw);

% polynomial fit for fast evaluation
% TODO there is actually quite a difference, how come?
if (1)
	%l  = double(cvec(linspace(min(level.val),max(level.val),100)));
	l  = double(cvec(linspace(min(z_b)+0.1,max(level.val),5)));
	a  = cvec(area_exact_f(l));
	p  = cvec(perimeter_exact_f(l));
	apoly = polyfit(l,a,np);
	%ppoly = polyfit(l,p,np);
	ppoly  = pchip(l,p);
else
	lmin = -max(z_b); 
	apoly  = polyfit([lmin; l0],[0; a0],np);
	%ppoly = polyfit([lmin; l0],[0; p0],np);
	ppoly  = pchip([lmin; l0],[0; p0]);
end

afunc_ = @(l) polyval(apoly, l);
%pfunc_ = @(l) polyval(ppoly, l);
pfunc_ = @(l) ppval(ppoly, l);

% regression of polynomial rating curve
t0 = calib.t0;
q0 = double(calib.discharge);
q0 = q0*sign(q0(1));
l  = level.val;
time = level.time;

gdx = isfinite(l);
fdx = ~gdx;
l_fix = l;
l_fix(fdx) = interp1(time(gdx),l(gdx),time(fdx),'linear');
dh_dt =[cdiff(l_fix)./cdiff(Constant.SECONDS_PER_DAY*time)];
% filter dh_dt
help = zeros(size(dh_dt));
% two day filter
T = 2;
d = median(diff(time));
for idx=1:length(dh_dt)
	fdx = max(1,floor(idx-T/d)):min(length(time),ceil(idx+T/d));
	w = hanwin((time(fdx)-time(idx))/T);
	w = w/sum(w);
	help(idx) = w'*dh_dt(fdx);
end
dh_dt = help;
dh_dt0 = interp1(time,dh_dt,calib.t0,'linear')
fdx = isnan(dh_dt0);
gdx = isfinite(dh_dt);
dh_dt0(fdx) = interp1(time(gdx),dh_dt(gdx),t0(fdx),'linear','extrap');

% Weighing with peak discharge does not reduce uncertainty at peak flow a lot
% we also do not weigh here, as n_eff = sum(w)^2/sum(w.^2) only approx. 3
% W = q0/norm(q0);
W = ones(size(q0));

%
% polynomial rating curve
%

polyrc = PolyRatingCurve(porder);
polyrc.fit(t0,q0,l0);
[qp_lin  void biasp_lin  errp_lin]   = polyrc.predict(time,l);

%
% exponential rating curve
%

% Note: the result is very sensitive with respect to the start values,
% as one would expect when regressing 3 parameters with only 4 calibrations
% the linear model as a start value is less good, empricial start parameters
%param0 = [300 0 1.5];
%param0 = [600 0 1.5];
flag = 0;

% rc = struct();

% power rating curve with free exponent
powerrc = PowerRatingCurve('afunc',afunc_,'pfunc',pfunc_);
powerrc.fit(t0,q0,l0);
if (flag) powerrc.jk.param.hat0 = powerrc.lin.param.hat; end
powerrc.jkfit(t0,q0,l0);
[qe_lin  void biase_lin  erre_lin]   = powerrc.predict(time,l);
[qe_jk  void biase_jk erre_jk]  = powerrc.jkpredict(time,l);

% Chezy: power rating curve with zero offset and ~ R^1/2
if (0)
chezyrc = PowerRatingCurve('afunc',afunc_,'pfunc',pfunc_,'c2',0.5);
chezyrc.fit(t0,q0,l0);
if (flag) chezyrc.jk.param.hat0 = chezyrc.lin.param.hat; end
powerrc.jkfit(t0,q0,l0);
[qe_chezy  void biase_chezy erre_chezy]  = powerrc_chezy.predict(time,l);
[qe_chezy_jk  void biase_chezy_jk erre_chezy_jk]  = chezyr.jkpredict(time,l);
end

% Q = A C R^1/2 S^1/2, C determined from R and roughness length z_0
keuleganrc = KeuleganRatingCurve('afunc',afunc_,'pfunc',pfunc_);
keuleganrc.fit(t0,q0,l0);
if (flag) keuleganrc.jk.param.hat0 = keuleganrc.lin.param.hat; end
keuleganrc.jkfit(t0,q0,l0);
[qk void biask errk]          = keuleganrc.predict(time,l);
[qk_jk void biask_jk errk_jk] = keuleganrc.jkpredict(time,l);

% Q = 1/n A R^5/6 S^(1/2)
% TODO, make this a powerratingcurve with c2 = 0 and c3 = 5/6
manningrc = ManningRatingCurve('afunc',afunc_,'pfunc',pfunc_);
manningrc.fit(t0,q0,l0);
if (flag) manningrc.jk.param.hat0 = manningrc.lin.param.hat; end
manningrc.jkfit(t0,q0,l0);
[qma void biasma errma] = manningrc.predict(time,l);

dpowerrc = DynamicPowerRC('afunc',afunc_,'pfunc',pfunc_);
dpowerrc.fit(t0,q0,l0, dh_dt0);
if (flag) dpowerrc.jk.param.hat0 = dpowerrc.lin.param.hat; end
%dpowerrc.jkfit(t0,q0,l0, dh_dt0);
[dqe_lin  void dbiase_lin  derre_lin] = dpowerrc.predict(time,l,dh_dt);
%[dqe_jk   void dbiase_jk   derre_jk]  = dpowerrc.jkpredict(time,l);
dqe_jk = [];
dbiase_jk = [];
derre_jk = [];

dmanningrc = DynamicManningRC('afunc',afunc_,'pfunc',pfunc_);
dmanningrc.fit(t0,q0,l0, dh_dt0);
if (flag) dmanningrc.jk.param.hat0 = dmanningrc.lin.param.hat; end
dmanningrc.jkfit(t0,q0,l0, dh_dt0);
[qdma void biasdma errdma] = dmanningrc.predict(time,l,dh_dt);

dmanningrc = DynamicManningRC('afunc',afunc_,'pfunc',pfunc_);
dmanningrc.fit(t0,q0,l0, dh_dt0);
if (flag) dmanningrc.jk.param.hat0 = dmanningrc.lin.param.hat; end
dmanningrc.jkfit(t0,q0,l0, dh_dt0);
[qdma void biasdma errdma] = dmanningrc.predict(time,l,dh_dt);

% dynamic ratinc curves (slope S dependent on dh/dt)
dkeuleganrc = DynamicKeuleganRC('afunc',afunc_,'pfunc',pfunc_);
dkeuleganrc.fit(t0,q0,l0,dh_dt0);
if (flag) dkeuleganrc.jk.param.hat0 = dkeuleganrc.lin.param.hat; end
dkeuleganrc.jkfit(t0,q0,l0,dh_dt0);
[qdk void biasdk errdk] = dkeuleganrc.predict(time,l,dh_dt);
[qdk_jk void biasdk_jk errdk_jk] = dkeuleganrc.jkpredict(time,l,dh_dt);



% dynamic rating curve
%S0 = 1e-5;
%w  = 600;
%dynrc = DynamicPowerRC(S0,w);
%dynrc.fit(t0,q0,l0,dh_dt0);
%[qd_lin  void biasd_lin  errd_lin]  = dpowerc.predict(time,l,dh_dt);

% export rating curve including coefficients as mat
t             = datestr(now(),'yyyy-mm-dd');
time          = level.time;
level_        = level;
level         = l;
discharge     = qma;
bias          = biasma;
if (isempty(bias))
	bias = NaN(size(qma));
end
serr	      = errma; 
%discharge     = qe_fix23;
%bias	      = biase_fix23;
%serr	      = erre_fix23;

%polyrc  = struct(polyrc);
%powerrc = struct(powerrc);
%save([base '.mat'],'time','level','discharge','bias','serr', ...
%			'polyrc','powerrc');
save([base '-' t '.mat'],'time','level','discharge','bias','serr', ...
			'dh_dt', ...
			'l0','q0', 'a0', 'p0', 'dh_dt0', ...
			'polyrc', 'qp_lin', 'errp_lin', ...
			'powerrc', 'qe_lin', 'erre_lin', 'qe_jk',  'erre_jk', ...
			'dpowerrc', 'dqe_lin', 'derre_lin', 'dqe_jk',  'derre_jk', ...
			'keuleganrc', 'qk','errk', 'qk_jk', 'errk_jk', ...
			'dkeuleganrc', 'qdk','errdk', 'qdk_jk', 'errdk_jk', ...
			'manningrc', 'qma', 'errma', ...
			'dmanningrc', 'qdma', 'errdma', ...
			... % 'powerrc_fix23', 'qe_fix23', 'erre_fix23', 'qe_fix23_jk', 'erre_fix23_jk', ...
			... % 'powerrc_fix3', 'qe_fix3', 'erre_fix3', ...
			... % ... 'dynrc', 'qd_lin', 'errd_lin', ...
			'afunc_', ...
			'pfunc_', ...
			'area_exact_f', ...
			'perimeter_exact_f', ...
			'z_b' ...
			);
system(['ln -s -f ',basename(base),'-',t,'.mat',' ',base,'.mat']);

% export rating curve as csv
dv  = datevec(time);
fid = fopen([base,'-',t,'.csv'],'w');
fprintf(fid,'yyyy mm dd HH MM H Q bias std\n');
fprintf(fid,'%4d %2d %2d %2d %2d %5.2f %5d %5d %5d\n', [dv(:,1:5) level round([discharge bias serr])]');
fclose(fid);
system(['ln -s -f ',basename(base),'-',t,'.mat ',base,'.mat']);

% tabulate
fprintf(1,'     expected estimated\n');
% 1/2 comes from C * sqrt(S)
param = powerrc.lin.param.val;
% TODO this is also not true should be w C sqrt(S0)
fprintf(1,'c_0: %f %f\n',mean(calib.width)*1/2, param(1));
% TODO 0 is not true, should be level - mean zb level
fprintf(1,'c_1: %f %f\n',            0, param(2));
%fprintf(1,'c_2: %f %f\n',          1.5, param(3));
%paramej
%disp('q- -qe0')
%q0e = interp1(time,qe,calib.t0,'linear');
%[q0 q0e]
%std(q0-q0e)/sqrt(length(q0)-3)
fprintf('serr0 %f\n',powerrc.lin.serr0);
fprintf('R2    %f\n',powerrc.lin.R2);

ret.time    = time;
ret.level   = level; % - meta.delta
%ret.qe_fix23  = qe_fix23;
ret.powerrc = powerrc;
ret.polyrc  = polyrc;
ret.calib   = calib;
%ret.dynrc   = dynrc;

end % sanggau_rating_curve()

