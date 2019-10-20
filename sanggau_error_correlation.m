% 2015-07-31 13:38:58.928394535 +0200
% Karl Kastner, Berlin

meta  = sanggau_metadata();
load([ROOTFOLDER,'/src/water-level/mat/water-level.mat']);
level.val  = K(nid.sanggau_merged).slot.depth;
level.time = K(nid.sanggau_merged).slot.time;

calib = load_vadcp_discharge(meta.filename.discharge,level);
N     = calib.N;
inner = abs(N) < meta.nmax;

load('mat/sanggau-calibration-error.mat');
trans = transversal;
%%% mode 2 : direct residuals of f_z and not those given by lienarisation are used
% mode 4 : prediction of f from ln_z0
vert  = vertical.roughness;
W = vert(1).W;

% constant and first order
for pdx=1:2
	zres = (vert(pdx).model.res_f(inner,:));
	tres = (trans(pdx).model.res(inner,:));

%	rho  = corr(zres,tres);
	rho_nz(pdx) = spearman_to_pearson(corr(flat(zres),flat(tres),'type','spearman'));
	rho_nz_(pdx) = spearman_to_pearson(corr(flat(zres*W'),flat(tres*W'),'type','spearman'));

	fprintf('%d-order Spearman correlation eps_fn eps_fz correlation %f\n',pdx-1,rho_nz(pdx));
end
rho_nz
rho_nz_
save('mat/sanggau-calibration-error.mat','-append','rho_nz');

