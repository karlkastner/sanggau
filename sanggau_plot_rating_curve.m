% Mon Nov  3 09:36:44 CET 2014
% Karl Kastner, Berlin

function sanggau_plot_rating_curve(pflag)

if(nargin() < 1)
	pflag = 0;
end

load('../discharge/mat/sanggau-stage-discharge.mat')

mcmc = load('mat/sanggau-mcmc.mat');
% preload water level
meta = sanggau_metadata();
% load water level
	load([ROOTFOLDER,'/src/water-level/mat/water-level.mat']);
	%level_.time = Kr(nid.sanggau_merged).time;
	%level_.val  = Kr(nid.sanggau_merged).depth-0*meta.hadcp.z0_;

% preload vadcp calibrarion discharge
	calib       = load_vadcp_discharge(meta); %,level_);
	level_.time = Kr(nid.sanggau_merged).time;
	level_.val  = Kr(nid.sanggau_merged).depth-meta.hadcp.z0_;

% correc for instrument elevation
%level_.val  = K(nid.sanggau_merged).slot.depth;

% some functions do not accept single input
q0 = double(q0);
% make discharge positive in downstream direction
sig = sign(q0(1));



nc = calib.nc;
abc  = char('A'+(0:nc-1)');
abc_ = [ones(calib.nc,1)*'  ', abc];
leg = {};
legQ = {};
for idx=1:calib.nc
	leg{idx} = sprintf('%s  %s',abc(idx), datestr(cvec(calib.tstart(idx)),'dd/mm/yyyy'));
	legQ{idx} = sprintf('%s  %s Q=%g', abc(idx), datestr(cvec(calib.tstart(idx)),'dd/mm/yyyy'),sig*round(calib.Q(idx)));
end

l_ = linspace(min(level),max(level),100)';
a_=area_exact_f(l_);
p_=perimeter_exact_f(l_);
time_ = [];

[qp_lin_ void biasp_lin_ errp_lin_]  = polyrc.predict(time_,l_);
[qe_lin_ void biase_lin_ erre_lin_]  = powerrc.predict(time_,l_);
%[qe_fix23_ void biase_fix23_ erre_fix23_]  = powerrc_fix23.predict(time_,l_);
%[qe_fix23_jk_ void biase_fix23_jk_ erre_fix23_jk_]  = powerrc_fix23.jkpredict(time_,l_);
%[qe_fix3_ void biase_fix3_ erre_fix3_]  = powerrc_fix3.predict(time_,l_);
[qe_jk_ void biase_jk_ erre_jk_]        = powerrc.jkpredict(time_,l_);
[qk_ void biask_ errk_]                 = keuleganrc.predict(time_,l_);
[qk_jk_ void biask_jk_ errk_jk_]        = keuleganrc.jkpredict(time_,l_);
[qdk_ void biasdk_ errdk_]              = dkeuleganrc.predict(time_,l_,0.*l_);
[qdk_jk_ void biasdk_jk_ errdk_jk_]     = dkeuleganrc.jkpredict(time_,l_,0.*l_);
%[qma void biasma errma epma]            = manningrc.predict(time,l_);
[qma_ void biasma_ errma_ epma_]        = manningrc.predict(time_,l_);
[qdma_ void biasdma_ errdma_]           = dmanningrc.predict(time_,l_,0.*l_);

if (0)
qtrue = qma;
qtrue_ = qma_;
param = [manningrc.lin.param.val 2/3];
R2 = manningrc.lin.R2;
serr0 = manningrc.lin.serr0;
else
qtrue  = qe_lin;
qtrue_ = qe_lin_;
param = powerrc.lin.param.val;
R2    = powerrc.lin.R2;
serr0 = powerrc.lin.serr0;
end

qe_fix23 = NaN(size(qtrue));
qe_fix23_ = NaN(size(l_));
%qe_fix23_jk = NaN(size(qtrue));
qe_fix3_ = NaN(size(l_));
erre_fix23_ = NaN(size(l_));

erre_lin = abs(erre_lin);

cc = colormap('lines');

namedfigure(103,'time derivative of stage');
plot(time,[dh_dt]); % help]);

namedfigure(104,'Pairwise parameter distribution');
clf
% plot ellipse for each parameter pair
%N = nchoosek(1:3,2);
N = nchoosek(1:2,2);
for idx=1:1
	subplot(1,3,idx);
	id = N(idx,1);
	jd = N(idx,2);
	% plot estimated parameter by linearisation
	p = powerrc.lin.param.val;
	C = powerrc.lin.param.C;
	c = [C(id,id), C(id,jd), C(jd,jd), p(id), p(jd)];
	hold on
	plot(p(id), p(jd),'.b');
%	plot_ellipse(c,[],'b');
	% jk
	p = powerrc.jk.param.val;
	C = powerrc.jk.param.C;
	c = [C(id,id), C(id,jd), C(jd,jd), p(id), p(jd)];
	hold on
	plot(p(id), p(jd),'.g');
%	plot_ellipse(c,[],'g');
	% mcmc
	% exluding "outliers" does not significantly change the result
	% s = (sort(sum((mcmc.chain-repmat(mean(mcmc.chain),size(mcmc.chain,1),1)).^2,2))); n=length(s);plot(s,(1:n)/n); q=quantile(s,0.85); id=find(s<q);
	p = mean(mcmc.chain); % theta, median?
	C = mcmc.res.cov; % == cov(chain)
	c = [C(id,id), C(id,jd), C(jd,jd), p(id), p(jd)];
	hold on
	plot(p(id), p(jd),'.r');
%	plot_ellipse(c,[],'r');

	xlabel(id)
	ylabel(jd)
	%axis equal
end

namedfigure(1,'Discharge from rating curve');
clf();
%subplot(2,1,1);
[ax h1 h2] = plotyy(time,level,time,qtrue/1e3); %[qp qe]);
%[ax h1 h2] = plotyy(level.time,level.val,level.time,[qp qe]);
set(h2,'linewidth',2)
set(ax(2),'linewidth',2) 
hold(ax(1),'on')
hold(ax(2),'on')
 %set(h1,'color','k')
 %set(h2,'color','k');
 %set(ax(1),'ycolor','k') 
 %set(ax(2),'ycolor','k');
%plot(ax(1),calib.tstart,l0,'bo','markerfacecolor','b');
%plot(ax(2),calib.tstart,q0/1e3,'o','color',cc(2,:),'markerfacecolor',cc(2,:));
% plot dynamic rating curve
%plot(ax(2),time,qd_lin,'color',cc(3,:));
%set(ax(1),'ylim',[0 15]);
set(ax(1),'ylim',[1 14]);
set(ax(1),'ytick',0:2:13);
set(ax(2),'ylim',[0 13]);
set(ax(2),'ytick',(0:2:15));
linkaxes(ax,'x');
vline(calib.tstart,'--k');
% the horizontal lines are confusing due to non-linear relation between stag and discharge
%hline([calib.lH, calib.lH - meta.hadcp.redeploy.d], 'k--')
%datetick(ax(1),'x');
set(ax(2),'xtick',[]);
tx  = monthspace(time(1),time(end));
dstr = datestr(tx,'mm/yy');
dstr(2:2:end,:) = ' ';
set(ax(1),'xtick',tx,'xticklabel',dstr);
xlim([min(time) max(time)]);
ylabel(ax(1),'z_s (m WGS84)');
ylabel(ax(2),'Q (10^3 m^3/s)');

namedfigure(10,'Time series of different rating curve models');
clf();
plot(time,[qe_fix23 qk qdk qma qdma]/1e3);
%,[erre_fix23_,errk_,errdk_,errma_, errdma_]/1e3);
legend('location','southeast','constant C','Keluaga','Dyn. Keluaga','Manning','Dyn Manning');

namedfigure(2,'Distribution and stage dependence of error and bias');
clf();

% some precomputation
relerre = abs(erre_lin./qe_lin);
relerrp = abs(errp_lin./qp_lin);
relerre_ = abs(erre_lin_./qe_lin_);
relerrp_ = abs(errp_lin_./qp_lin_);
l_m   = mcmc.out.data{1};
qm_   = mcmc.out.predlims{1}{1}(3,:)';
low   = mcmc.out.predlims{1}{1}(2,:)';
up    = mcmc.out.predlims{1}{1}(4,:)';
errm_ = 0.5*(up-low);
errm = interp1(l_m,errm_,level,'linear');

% CDF of discharge error
% Note: for the distributions the entire time series has to be considered, not the resampled
errp_lin = double(errp_lin);
erre_lin = double(erre_lin);
%l    = double(l);
subplot(2,2,1);
[x void y] = hist_man(erre_lin,100);
plot(x,y,'-','color',cc(1,:));
hold on
[x void y] = hist_man(erre_jk,100);
plot(x,y,'-','color',cc(2,:));
[x void y] = hist_man(errm,100);
plot(x,y,'-','color',cc(3,:));
[x void y] = hist_man(errp_lin,100);
plot(x,y,'-','color',cc(4,:));

legend('location','southeast','exp lin','exp jn','exp mcmc','poly jn');
xlabel('std_q (m^3/s)');
title('CDF of discharge standard error');

% Stage dependence of discharge error
subplot(2,2,2)
plot(l_,[erre_lin_ erre_jk_ errm_ errp_lin_],'-');
hold on
xlim([0 max(l_)]);
ylim([0 1000]);
vline(l0);
xlabel('z_s (m WGS84)');
ylabel('std_q (m^3/s)')
title('Stage dependence of discharge error');

subplot(2,2,3);
%[x void y] = hist_man(biasp,100);
%plot(x,y,'-','color',cc(3,:));
hold on
%[x void y] = hist_man(biase,100);
%plot(x,y,'-','color',cc(3,:));
title('CDF of bias of uncorrected discharge');
xlabel('bias_q (m^3/s)');
% stage vs bias
subplot(2,2,4)
%plot(l_,[biasp_ biase_ biasmc_],'-','color',cc(3,:));
xlim([0 max(l_)]);
ylim(500*[-1 1])
vline(l0);
title('Stage dependence of bias of uncorrected discharge');
xlabel('z_s (m WGS84)');
ylabel('bias_q (m^3/s)')

namedfigure(3,'Rating curve');
clf();
subplot(2,2,1);
%errorlines(l_,[qe_lin_,qp_lin_ qe_fix23_ qe_fix3_]/1e3,[errp_lin_,erre_lin_,erre_fix23_,erre_fix3_]/1e3);
%errorlines(l_,[qe_lin_,qp_lin_ qe_fix23_ qe_fix3_ qk_ qma_ qk_dyn_]/1e3,[errp_lin_,erre_lin_,erre_fix23_,erre_fix3_,errk_,errma_, errk_dyn_]/1e3);
errorlines(l_,[qe_fix23_ qk_ qdk_ qma_ qdma_]/1e3,[erre_fix23_,errk_,errdk_,errma_, errdma_]/1e3);
hold on
plot(l0,q0/1e3,'ok','markerfacecolor','k');
xlabel('z_s (m WGS84)');
ylabel('Q (1e3 m^3/s)');
xlim([0 l_(end)]);
ylim([0 qe_lin_(end)/1e3*1.05]);
vline(l0);
%vline(calib.lH);
%vline(calib.lH - meta.hadcp.redeploy.d)
title('Rating curve');
%legend('location','southeast','exponential','polynomial','offset and exp fixed','exp fixed','keluegan','manning','dyn keluegan');
legend('location','southeast','constant C','Keluaga','Dyn. Keluaga','Manning','Dyn Manning');

namedfigure(30,'Rating curve');
clf();
%dots(l0,q0/1e3,'ok')
plot(l0,q0/1e3,'ok','markerfacecolor','k');
hold on
% line plot must follow dots for legend
plot(l_,qtrue_/1e3,'k');
%vline([calib.lH, calib.lH - meta.hadcp.redeploy.d],'k--');
%xlim([0 l_(end)]);
xlim(limits(l_));
ylim([0 qe_lin_(end)/1e3*1.05]);
%vline([calib.lH - meta.hadcp.redeploy.d],'k--');
xlabel('z_s (m WGS84)');
ylabel('Q (10^3 m^3/s)');
text(6,9,sprintf('Q = %1.2f A{R^{%1.2f}}\nR^2= %1.3f\nRMSE=%3d m^3/s', param, R2, round(serr0)));
vline(meta.hadcp.z0,'--')
text(meta.hadcp.z0+0.25,4,'HADCP level','rot',90)
%text(l0,q0/1e3,abc_);
text(l0+0.15,q0/1e3-0.15,abc);

[void sdx] = sort(calib.tstart);
C = {};
C{1,1} = ' ';
C{1,2} = ' Date';
C{1,3} = '   Discharge';
C(2:nc+1,1)=num2cell(abc(sdx));
C(2:nc+1,2)=cellfun(@(x) datestr(x,'dd/mm/yyyy'),num2cell(cvec(calib.tstart(sdx))),'uniformoutput',false);
C(2:nc+1,3)=cellfun(@num2str,num2cell(int32(q0(sdx))),'uniformoutput',false);
C = C';
%text(9,2,sprintf('%c %+7s %+15s\n',C{:}))
text(9,2,sprintf('%-3c %-12s %-s\n',C{:}))


namedfigure(4,'Rating curve');
clf();
errorlines(l_,[qe_lin_ qe_jk_ qm_]/1e3,[erre_lin_,erre_jk_,errm_]/1e3);
vline(l0);
%hline(calib.lH);
%hline(calib.lH - meta.hadcp.redeploy.d)
hold on;
plot(l0,q0/1e3,'ok','markerfacecolor','k');

%plot(l_,qe_lin_/1e3,'-','color',cc(1,:));
%hold on
%plot(l_,qe_jk_/1e3,'-','color',cc(2,:));
%plot(l_,qe_lin_/1e3*[1 1]+2*erre_lin_/1e3*[-1 1],'--','color',cc(1,:));
%plot(l_m,q_m/1e3,'-','color',cc(3,:));
%plot(l_,qe_jk_/1e3*[1 1]+2*erre_jk_/1e3*[-1 1],'--','color',cc(2,:));
%plot(l_m,1e-3*([low up]),'--','color',cc(3,:));
legend('location','southeast','Linear','JN','MCMC');
xlabel('z_s (m WGS84)')
ylabel('q (m^3/s)')
%{
subplot(2,2,2);
plot(l_,relerre_,'-','color',cc(2,:));
hold on;
plot(l_,relerrp_,'-','color',cc(3,:));
title('Relative error')
ylabel('s_q/\mu_q');
xlabel('z_s (m WGS84)');
xlim([0 l_(end)]);
ylim([0 max([relerre_; relerrp_])*1.05]);
vline(l0);
pdfprint('img/sanggau-rating-curve-2.eps')
%}

namedfigure(100,'Individual parameter density');
clf();
rc = powerrc;
tit = {'$Scale (C \sqrt{S} w)$', 'Offset', 'Exponent'}
for idx=1:2
	subplot(2,2,idx);
	gauss = 1;
	y = [];
	[y(:,3),x] = density(mcmc.chain(:,idx),[],gauss);
	hold on;
	y(:,1) = normpdf(x,rc.lin.param.val(idx),sqrt(rc.lin.param.C(idx,idx)));
	y(:,2) = normpdf(x,rc.jk.param.val(idx),sqrt(rc.jk.param.C(idx,idx)));
	plot(x,y);
	title(tit{idx},'interpreter','latex');
end

namedfigure(102,'Cross section area and perimeter');
clf();
%[ax h1 h2] = plotyy(l_,a_/1e3,l_,p_);
plot(l_,a_/1e3);
ax = gca();
hold(ax(1),'on')
ax(2)=addy();
hold(ax(2),'on')
plot(ax(2),l_,p_,'r');
set(ax(2),'xaxislocation','top','yaxislocation','right','ycolor','r');
set(ax(2),'xticklabel','');

plot(ax(1),l0,a0/1e3,'ko','markerfacecolor','k');
plot(ax(2),l0,p0,'ro','markerfacecolor','r');
ylabel(ax(1),'Area (10^3 m^2)');
ylabel(ax(2),'Perimeter (m)');
xlabel(ax(1),'z_s (m WGS84)');
set(ax(2),'ytick',400:20:800)
ylim(ax(2),[550 700]);
xlim(ax(1),limits(l_))
xlim(ax(2),limits(l_))
vline(meta.hadcp.z0,'k--')
text(meta.hadcp.z0+0.25,4,'HADCP level','rot',90,'parent',ax(1))
text(l0+0.1,a0/1e3-0.1,abc,'parent',ax(1));
text(l0+0.1,p0-5,abc,'parent',ax(2));

% date P A
[t0_ sdx] = sort(cvec(calib.tstart));

C = {};
C{1,1} = ' ';
C{1,2} = ' Date';
C{1,3} = '    Area';
C{1,4} = ' Perimeter';
C(2:nc+1,1)=num2cell(abc(sdx));
C(2:nc+1,2)=cellfun(@(x) datestr(x,'dd/mm/yyyy'),num2cell(cvec(calib.tstart(sdx))),'uniformoutput',false);
C(2:nc+1,3)=cellfun(@num2str,num2cell(int32(a0(sdx))),'uniformoutput',false);
C(2:nc+1,4)=cellfun(@num2str,num2cell(int32(p0(sdx))),'uniformoutput',false);
C = C';
text(8.5,3,sprintf('%-3c %-12s %-5s %-5s\n',C{:}),'parent',ax(1))
set(ax(1),'box','off')
set(ax(2),'box','off')
set(ax(2),'XaxisLocation','top')

namedfigure(101,'Stationary vs instationary rating curve');
subplot(2,2,1)
%hist_man((qd_lin - qe_lin)/1e3, 25); %3,100,'.');
title('Absolute');
subplot(2,2,2)
%hist_man(2*(qd_lin - qe_lin)./(qd_lin+qe_lin),25,'-');
title('Relative');
%fprintf(1,'median(q_dyn - q_stat): %f 10^3 m^3/s\n', nanmedian((qd_lin - qe_lin)/1e3));
%fprintf(1,'1/2 q84 - q16 | (q_dyn - q_stat): %f 10^3 m^3/s\n', 0.5*diff(quantile(qd_lin - qe_lin,normcdf([-1 1])))/1e3);
%fprintf(1,'mean(q_dyn - q_stat): %f 10^3 m^3/s\n', nanmean((qd_lin - qe_lin)/1e3));
%fprintf(1,'std(q_dyn - q_stat): %f 10^3 m^3/s\n', nanstd((qd_lin - qe_lin)/1e3));
%fprintf(1,'median 2*(q_dyn - q_stat)/(q_d + q_s): %f/s\n', nanmedian(2*(qd_lin - qe_lin)./(qd_lin+qe_lin)));
%fprintf(1,'mean 2*(q_dyn - q_stat)/(q_d + q_s): %f \n', nanmean(2*(qd_lin - qe_lin)./(qd_lin+qe_lin)));
%fprintf(1,'std 2*(q_dyn - q_stat)/(q_d + q_s): %f \n', nanstd(2*(qd_lin - qe_lin)./(qd_lin+qe_lin)));
%fprintf(1,'1/2 q84 - q16 | 2 (q_dyn - q_stat): %f \n', diff(quantile((qd_lin - qe_lin)./(qd_lin+qe_lin),normcdf([-1 1]))));
fprintf('\n\n\n');
%p23 = powerrc_fix23;
k   = keuleganrc;
dk  = dkeuleganrc;
m   = manningrc;
dm  = dmanningrc;
field = {'lin','jk'}
for idx=1:2
f = field{idx};
fprintf('            C S_0^1/2, z_0 (mm), S_0 (10^-5), s_err (m^3/s), R^2  AIC\nc')
%fprintf('constant C    %0.3f %1.2f %1.2f %3.0f %0.3f %2.1f %2.1f\n', p23.(f).param.val(1),                 NaN,                      NaN,   p23.(f).serr0, p23.(f).R2, p23.(f).aic, p23.(f).bic);
fprintf('Keulegan      %0.3f %1.2f %1.2f %3.0f %0.3f %2.1f %2.1f\n', NaN,               k.(f).param.val(1)*1e3,   k.(f).param.val(2)*1e5, k.(f).serr0, k.(f).R2, k.(f).aic, k.(f).aic);
fprintf('dyn. Keulegan %0.3f %1.2f %1.2f %3.0f %0.3f %2.1f %2.1f\n', NaN,              dk.(f).param.val(1)*1e3,  dk.(f).param.val(2)*1e5, dk.(f).serr0,dk.(f).R2,dk.(f).aic,dk.(f).bic );
fprintf('Manning       %0.3f %1.2f %1.2f %3.0f %0.3f %2.1f %2.1f\n', m.(f).param.val(1),   NaN,            NaN, m.(f).serr0,m.(f).R2,m.(f).aic,m.(f).bic );
fprintf('dyn Manning   %0.3f %1.2f %1.2f %3.0f %0.3f %2.1f %2.1f\n', NaN,dm.(f).param.val(1)*1e3,  dm.(f).param.val(2)*1e5, dm.(f).serr0,dm.(f).R2,dm.(f).aic,dm.(f).bic );

end

corrcov(powerrc.lin.param.C)
figure(20);
hist_man(dh_dt,100);
[[max(l0) min(l0)]-min(level)]/(max(level)-min(level))
'area polynomial'
fh = functions(afunc_);
fh.workspace{1}.apoly
'perimeter polynomial'
fh = functions(pfunc_);
fh.workspace{1}.ppoly

namedfigure(21,'Relatice predicition error of Manning rating curve');
clf();
dots(manningrc.h0,abs(manningrc.lin.res0./manningrc.q0));
hold on
%col = colormap('lines');
%for idx=1:length(calib.tstart)
%	plot(manningrc.h0(idx),abs(manningrc.lin.res0(idx)./manningrc.q0(idx)),'ko','markerfacecolor',col(idx,:))
%	hold on
%	%text(l0(idx),double(abs(manningrc.lin.res0(idx)./manningrc.q0(idx))), ['  ',abc(idx)]);
%end
text(l0,q0/1e3,abc_);

plot(l_,epma_./qtrue_,'k');
xlabel('z_s (m WGS84)');
percenttick('y');
ylabel('$\frac{s_{Q_R}}{Q_{R}}$ (\%)','interpreter','latex','rot',0);
xlim(limits(l_));
%hline(,'--');
%vline([calib.lH - meta.hadcp.redeploy.d],'k--');
ylim([0 max(epma_./qtrue_)]);


fprintf('manning r2 %f serr %f c %f\n',manningrc.lin.R2,manningrc.lin.serr0,manningrc.lin.param.val);

if (pflag)
	ps = [2 2.5 3];
	for idx=1:length(ps)
		scale = ps(idx)
		pdfprint(  1,'img/sanggau-rating-curve-2-time-series',scale-1,1/2);
		pdfprint( 30,'img/sanggau-rating-curve',scale);
		pdfprint(102,'img/sanggau-area-perimeter-vs-stage',scale);
		pdfprint( 21,'img/sanngau-rc-prediction-error',scale);
	end
	pdfprint(  2,'img/sanggau-rating-curve-2-errors');
	pdfprint(100,'img/sanggau-rc-parameter-distribution');
	pdfprint(101,'img/sanggau-rc-stationary-vs-instationary');
	pdfprint(104,'img/sanggau-rc-pairwise-parameter-distribution');
end

end % function

