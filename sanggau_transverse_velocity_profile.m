% Fr 31. Jul 11:04:11 CEST 2015
% Karl Kastner, Berlin

if (~exist('pflag','var'))
	pflag = 0;
end
fflag = pflag;

% TODO no magic numbers
nf_  = 2;
nf__ = 100;
% TODO use lanczos, global fiter parameter
%win = lanczoswin(1:round(1/0.3062*nf__))';
win = rectwin(1:nf__);
%win = triwin(1:round(nf__/0.75));
dflag = 1;
%win = hanwin(1:round(3/2*nf__))';
%win = rectwin(nf__,1:round(3/2*nf__))';
if (~exist('reload','var') || reload)
	meta = sanggau_metadata();
	load([ROOTFOLDER,'/src/water-level/mat/water-level.mat']);
	level.val  = Kr(nid.sanggau_merged).depth;
	level.time = Kr(nid.sanggau_merged).time;
	calib = load_vadcp_discharge(meta);
	[hadcp hlevel offset] = sanggau_load_hadcp(meta);
	reload = 0;
end
hadcp_range = meta.hadcp.range;

%cs  = calib.cs_;
% bring into same order
%t = arrayfun(@(x) x.t0(1), cs)';
%f=@(x) x.time(1);
%t=arrayfun(f,adcp_)';
%cs_ = CrossSection();
[t sdx] = sort(calib.tstart);
cs = calib.cs_(sdx);
U_tn = calib.U_tn(:,sdx);
U_t = calib.U(sdx);
%for idx=1:length(t)
%	[mv mdx] = min(abs(calib.tstart(idx)-t))
%	cs_(idx) = cs(mdx);
%	%calib.cs(idx).adcp = adcp_(mdx);
%end
%cs = cs_;
% t0 = arrayfun(@(x) x.t0(1),cs);
% [t0 sdx] = sort(t0);
% cs = cs(sdx); 
%error sort calib also

N     = calib.N;
%gridN.cX1;

% concatenate transect data off all campaigns into a single matrix

% transverse velocity profile
fn_  = [];
w_   = [];
h_   = [];
id_  = [];
resi = [];
as2  = [];
as2.individual.acov = [];
W = [];
E = [];
fn   = bsxfun(@times,U_tn,1./U_t.');
fn_  = meanfilt1(fn,nf_);
for idx=1:length(cs)
	% cs(idx).gridN.val.U./cs(idx).U;
	%Y = cs(idx).pseudo.gridN.U./cs(idx).U;
	% TODO use a low pass filter with hanwin here
	nt   = cs(idx).transect.n;
	if (0 == nt)
		nt = 1;
	end
	id  = [idx*ones(1,nt)];
	id_ = [id_, id];
	w_  = [w_, ones(1,nt)/nt];
	W(end+1,end+1:end+nt) = 1/nt;
	E(end+1,end+1:end+nt) = 1;
	fdx = isfinite(level.val);
	h  = interp1(level.time(fdx),level.val(fdx),mean(cs(idx).t0),'nearest');
	h_(idx,1) = h;
%	h_(idx,1)  = calib.h0(idx)*ones(nt,1);
end
% make fn_ positive (sign depends on cross section definition)
fn_ = sign(nanmedian(flat(fn_)))*fn_;

w_ = w_/sum(w_);

%namedfigure(1,'Transversal velocity profile and ACF');
%clf();
tit = {'constant','linear'};
for pdx=1:2
	transverse(pdx) = transverse_profile_parameter(N, fn_, h_, w_, W, E, meta.nmax, nf_, pdx-1);
	% plot profile
	splitfigure([2 2],[1 pdx],fflag);
	%subplot(2,2,pdx)
	plot(N,fn_);
	hold on
	set(gca,'colororderindex',1)
	plot(N,transverse(pdx).poly.predict(h_)','--');
	% plot acf
	%subplot(2,2,2+pdx);
	splitfigure([2 2],[1 2+pdx],fflag);
	trans = transverse(pdx);
	a   = [trans.total.acov, trans.camp.acov trans.model.acov];
	a   = [a trans.model.s2*acfar1(trans.model.rho,size(a,1),1:size(a,1))];
	plot(a,'.');
	hline(a(1,:)*exp(-1))
	title(tit{pdx});
end
N_ = N;

save('mat/sanggau-calibration-error.mat','-append','transverse');

v=transverse;
figure(20);
clf
n = 500;
N = (1:n-1)';
plot([v(1).model.acov v(2).model.acov],'.')
hold on
set(gca, 'ColorOrderIndex', 1)
plot(N,[v(1).model.s2*acfar1(v(1).model.rho,n-1,1:(n-1)) v(2).model.s2*acfar1(v(2).model.rho,n-1,1:(n-1))])
legend('measured constant','measured linear','model constant', 'model linear');
%ylabel('$E[x_i x_{i+l}]$','interpreter','latex','rot',0);
ylabel('$a_t$','interpreter','latex','rot',0);
xlabel('l (m)');

namedfigure(30,'Prediction error of the tranversal velocity profile');
clf();
serr = [];
for pdx=1:size(fn_,2)
	fn_(:,pdx)   = wmeanfilt(win,fn_(:,pdx),dflag);
	%fn_(:,pdx)   = meanfilt1(fn_(:,pdx),nf__);
end
fn__ = [];
for pdx=1:2
	transverse(pdx) = transverse_profile_parameter(...
				N, fn_, h_, w_, W, E, meta.nmax, nf__, pdx-1);
	%serr(:,pdx)     = transverse(pdx).model.serr;
	relrms(:,pdx)     = transverse(pdx).model.relrms;
	%fn__(:,pdx)     = transverse(pdx).poly.predict(h_);
end


f=transverse(2).poly.predict([min(h_),midrange(h_),max(h_)]')';
if (0)
figure(100);
plot(N_,100*relrms);
vline(hadcp.N(1),'k--');
vline(hadcp.N(1)-hadcp_range,'k--');
xlabel('N (m)');
ylim([0 20]);
xlim(330*[-1 1]);
lh = legend('constant','linear');
ylabel('Relative prediction error (%)');
%set(get(lh,'title'),'string','Prediction order');

hold on
ax(1)=gca;
ax(2)=addx();
plot(ax(2),N_,f);
linkaxes(ax,'x')
else

p = [2 1 3];
%[ax h1 h2] = plotyy(N_,f(:,p),N_/range(N_),100*relrms);
%[ax h1 h2] = plotyy(N_,f(:,p),N_/range(N_),100*relrms);
h1 = plot(N_,f(:,p));
hold on
ax = gca
% ylim has to precede vline
ax(2) = addy();
axoff=ax;
hold on
h2 = plot(ax(2),N_/range(N_),100*relrms);

%[ax h1 h2] = plotyy(N_,100*relrms,N_,f(:,p));
set(h1,'color','k');
set(h2,'color','r');
ls = {'-','--','-.'}
for idx=1:length(h1)
	set(h1(idx),'linestyle',ls{idx});
end
for idx=1:length(h2)
	set(h2(idx),'linestyle',ls{idx});
end

drawnow
if (0)
	lh = legend('location','best','f_t mid','f_t low', 'f_t high','RMS constant','RMS linear');
	set(lh,'position',[0.4 0.3 0.2 0.2])
end

set(ax(1),'box','off');
set(ax(1),'ytick',0.6:0.2:1.2,'ycolor','k');
xlabel(ax(1),'N (m)');
xlim(ax(1),330*[-1 1]);
ylabel(ax(1),'Velocity profile $f_t = \bar u / U$','interpreter','latex');
ylim(ax(1),[0.6 1.2]);
vline(ax(1),hadcp.N(1),'k--');
vline(ax(1),hadcp.N(1)-hadcp_range,'k--');


ax(2).XAxisLocation = 'top';
ax(2).YAxisLocation = 'right';
set(ax(2),'xtick',[-0.4:0.2:0.4]);
set(ax(2),'ytick',5:5:20,'ycolor','r');
xlabel(ax(2),'\xi')
xlim(ax(2),[-1 1]/2);
ylab = 'Relative prediction error $\frac{1}{\hat f_t}\sigma_{\varepsilon_{\hat f_t}}$ (\%)',
ylabel(ax(2),ylab,'interpreter','latex');
ylim(ax(2),[0 20]);


%set(get(lh,'title'),'string','Prediction order');
%linkaxes(ax,'x')
end

if (pflag)
	ps = [2 2.5 3];
	for idx=1:length(ps)
		pdfprint(20,'img/acf-eps-f_xi',ps(idx));
		pdfprint(30,'img/sanggau-u-div-U-vs-n-prediction-error',ps(idx));
		% ylab = get(axoff(2),'ylabel');
		ylabel(axoff(2),'');
		pdfprint(30,'img/sanggau-u-div-U-vs-n-prediction-error-no-right-ylabel',ps(idx));
		ylabel(axoff(2),ylab,'interpreter','latex');
	end
	pflag = 0;
end

