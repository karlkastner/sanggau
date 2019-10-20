% Wed Mar 18 10:01:07 CET 2015
% Karl Kastner, Berlin

if (~exist('reload','var') || reload)
meta = sanggau_metadata();
% get HADCP velocity
[hadcp hlevel]= sanggau_load_hadcp(meta,[],1800);
hl.time = hadcp.time;
hl.val  = hlevel;
calib = load_vadcp_discharge(meta.filename.discharge,hl);
% rotate hadcp velocity do the cross section
hadcp.velocity.cs = rotvel(calib.dir, hadcp.velocity.earth);

% load rating curve discharge
rc = load('../discharge/mat/sanggau-stage-discharge.mat');
reload = 0;
end
% TODO quick fix
calib.bottom.median(end) = calib.bottom.median(end-1);

figure(1);
clf();
subplot(2,2,1)
% select time spans of the HADCP to use
plot(hadcp.time,hlevel);
%plot(hadcp.time,hadcp.idepth_m);
hline(0);
hline(-meta.hadcp.redeploy.d);

T = 1e5*[	7.358279430292825,   7.358563805292825;
		7.359471439499999,   7.360159772833333];

% select valid cells
fdx = (hadcp.time > T(1,1) & hadcp.time < T(1,2)) | (hadcp.time > T(2,1) & hadcp.time < T(2,2));

time  = hadcp.time(fdx);
%level = hlevel - calib.lH;
level = hlevel(fdx);
U = hadcp.velocity.cs(:,fdx,1);
E = hadcp.E(:,fdx,:);
um1  = nanmedian(U,1);
um2  = mean(U,2);
U_ = bsxfun(@minus,U,um1);
us2  = std(bsxfun(@minus,U,um1),[],2);
subplot(2,2,2)
errorlines(1:length(um2),um2,us2);

% calibrate the IVM

% interpolate discharge and area to HADCP samples
%qrc  = interp1(rc.time,rc.qe_fix23,time);
qrc  = interp1(rc.time,rc.qma,time);
%area = interp1();

load('mat/fgood.mat');

ivm_ma = IVM();
ivm_ma.fit(level(fgood),U(:,fgood)',qrc(fgood),[],calib.bottom.median, [], [], []);
qh  = ivm_ma.predict(level,U');
fprintf('qh-qivm: %f\n',mean(qh(fgood)-qrc(fgood)))

ivm_ = IVM('h');
level_ = cvec(level);
a      = afunc(level_,calib.bottom.median,calib.dw);
ma     = max(level_(fgood));
mi     = min(level_(fgood));
%level_ = level_/(ma-mi);
%u      = mean(U(:,fgood))';
%level_ = (level_ - mean(u.*a(fgood).*level_(fgood))/(mean(u.*a(fgood))));
level_ = (level_ - 0.5*(ma+mi))/(ma-mi);
ivm_.fit(level(fgood),U(:,fgood)',qrc(fgood),[],calib.bottom.median, [], [], [], level_(fgood));
qh_  = ivm_.predict(level,U',level_);

qdk  = interp1(rc.time,rc.qdk,time);
ivm = IVM();
ivm.fit(level(fgood),U(:,fgood)',qdk(fgood),[],calib.bottom.median, [], [], []);
qh_dk  = ivm.predict(level,U');

% compare discharges
subplot(2,2,3)
plot(time(fgood),[qrc(fgood) qdk(fgood) qh(fgood) qh_dk(fgood)],'.');

subplot(2,2,4)
plot(qh(fgood)/1e3, qrc(fgood)/1e3,'.b','markersize',5);
hold on;
plot(qh_dk(fgood)/1e3, qdk(fgood)/1e3,'.g','markersize',5);
plot([0 1e4], [0 1e4],'k');
hold off

figure(3)
%subplot(2,2,4)
%plot(qh(fgood)/1e3, qrc(fgood)/1e3,'.k','markersize',5);
plot(qrc(fgood)/1e3, [qh(fgood) qh_(fgood)]/1e3,'.','markersize',5);
hold on;
plot([0 10], [0 10],'k');
hold off
ylabel('Q_{hadcp} (10^3 m^3/s)');
xlabel('Q_{manning} (10^3 m^3/s)');
%plot(qrc, qh);
pdfprint('img/sanggau-discharge-rc-vs-ivm');

figure(2);
clf();
hist_man(qh-qrc);
hold on;
hist_man(qh_dk-qdk,[],'r');

if (0)
figure(3)
clf
imagesc(U)

figure(4)
clf()
for idx=1:3; subplot(3,1,idx); em1 = nanmedian(E(:,:,idx),1); plot(bsxfun(@minus,E(:,:,idx),em1)); end

figure(5)
clf()
for idx=1:3; subplot(3,1,idx); em1 = nanmedian(E(:,:,idx),1); plot(bsxfun(@minus,E(:,:,idx),em1)); end

figure(6)
clf()
C = hadcp.dat.CORR;
for idx=1:3; subplot(3,1,idx); em1 = nanmedian(E(:,:,idx),1); plot(bsxfun(@minus,C(:,:,idx),0)); end
end

d   = double(qh_dk - qdk);
fgood = abs(d-medfilt1(d,250)) < 550;
save('mat/fgood.mat','fgood');
%d_=d;
%d_(fdx)=NaN; plot([d d_]);

1-var(d(fgood))/var(qdk)
[level(fgood)'.^0 level(fgood)'] \ d(fgood)
[qh_dk(fgood) level(fgood)']\qdk(fgood)

d   = double(qdk - qrc);
1-var(d(fgood))/var(qrc)

d      = double(qh - qrc);
rmse   = sqrt(mean(d.^2));
R2     = 1-var(d(fgood))/var(qrc);
level_ = cvec(level) + calib.lH;
[level_(fgood).^0 level_(fgood)] \ d(fgood)
fprintf('rmse %f\n',rmse);
fprintf('R2   %f\n',R2);
fprintf('ivm param: %f\n',ivm_ma.param.val0);
reg    = Theil();
reg.regress(level_(fgood),d(fgood));
reg.param

