% 2015-02-23 21:42:52.527937395 +0100
% Karl Kastner, Berlin
% uses toolbox of H. Haario, M. Laine, A. Mira and E. Saksman, 2006. DRAM: Efficient adaptive MCMC, Statistics and Computing 16, pp. 339-354.

addpath([ROOTFOLDER,filesep(),'src-external/mcmc/stat']);

% MCMC settings
options.nsimu       = 1e4;
options.updatesigma = 1;

fprintf(1,'Loading data\n');

if (~exist('reload','var') || reload)

% preload water level
meta = sanggau_metadata();
% preload water level
	load([ROOTFOLDER,'/src/water-level/mat/water-level.mat']);
	level.val  = K(nid.sanggau_merged).slot.depth;
	level.time = K(nid.sanggau_merged).slot.time;

% preload vadcp calibrarion discharge
calib = load_vadcp_discharge(meta.filename.discharge,level);

% correct for average depth
% TODO this should be part of load calib and water level
h0        = calib.l0  + meta.lH;
level.val = level.val + meta.lH;
%x = linspace(min(rc.h0),max(rc.h0),100)';
x         = linspace(min(level.val),max(level.val),100)';

% regression of polynomial rating curve
h0 = double(h0);
q0 = double(calib.cs.q0);

reload = 0
end


data.xdata = h0;
data.ydata = q0;

rc = PowerRatingCurve();

% fit initial parameter
fprintf(1,'Fitting initial parameters\n')
rc.fit(h0,q0);

% initial paramter
tmin = rc.param.val0;

% initial covariance matrix
options.qcov = rc.lin.C;

% initial error variance
model.sigma2 = rc.lin.serr0.^2;

% objective function
model.ssfun    = @(c,void) sum((rc.objective(c',[],h0,q0)).^2);

% initial parameters
% 'name', initial_value, min_value, max_value, prior_mu, prior_sigma, and
% targetflag

params = {
    {'c1', tmin(1), -Inf},
    {'c2', tmin(2), -Inf},
    {'c3', tmin(3), -Inf}
%   {'c3', tmin(3), 0}
%    {'c3', 2, 0}
    };

fprintf(1,'Generating markov chain\n')

% generate MCMC chain
[res,chain,s2chain] = mcmcrun(model,data,params,options);

fprintf(1,'Prediction uncertainty');
modelfun = @(h,c) rc.rcfunc(c,[],h);
out = mcmcpred(res,chain,[],x,modelfun);

oname = ['sanggau-mcmc-' datestr(now,'yyyy-mm-dd-HH-MM-SS') '.mat'];
save(['mat' filesep() oname],'res','chain','s2chain','model','data','options','out','rc');
system(['ln -s -f ' oname ' mat/sanggau-mcmc.mat']);

