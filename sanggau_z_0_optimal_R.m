% 2014-11-12 18:15:51.070651560 +0100
% Karl Kastner, Berlin

base = '../discharge/mat/2013-12-09-sanggau';
base = '../discharge/mat/2014-02-20-sanggau';
%2014-02-20

R = 2.^(-1:5)'; % 8 for publication
%zerr = [];
%zval = [];
for idx=1:length(R)
	fprintf('progress: %d%%\n',round(100*(idx-1)/length(R)));
	fname = [base '-Rz-' num2str(R(idx),'%05.1f') '.mat'];
	if (exist(fname,'file'))
		load(fname);
	else
		load([base '.mat']);
		discharge.jflag = 0;
		discharge.R_z0 = R(idx);
		discharge.regress_log_profile();
		save(fname,'discharge');
	end
	dis(idx) = discharge;
%	discharge.R_b = R(idx);
%	discharge.regress_bottom_profile();
%	zval(:,idx) = discharge.bgrid.val(:,1);
%	zerr(:,idx) = discharge.bgrid.err(:,1);

	zval(:,idx) = discharge.z0grid.val.ln_z0(:,1);
	zerr(:,idx) = discharge.z0grid.err.ln_z0(:,1);
%	zval_(:,idx) = discharge.z0grid.val.ln_z0(:,2);
%	zerr_(:,idx) = discharge.z0grid.err.ln_z0(:,2);
	mse0(:,idx) = discharge.z0grid.err.mse0(:,1);
	mse1(:,idx) = discharge.z0grid.err.mse1(:,1);
	rho(:,idx,:) = [discharge.z0grid.err.rho(:,1:3);];
	f(:,idx,:) = [discharge.z0grid.err.f(:,1:2)];
	U(:,idx) = discharge.z0grid.val.u_s/Constant.KAPPA.*(log(discharge.bgrid.val(:,1)) - 1 - discharge.z0grid.val.ln_z0);
end

namedfigure(1,'ln_z0');
subplot(3,1,1);
plot(R,nanmedian(zval),'o');
title('Median value')
subplot(3,1,2);
fdx_ = isfinite(sum(zerr,2));
plot(R,quantile(zerr(fdx_,:),normcdf(-2:2)),'o');
title('Error quantiles')

subplot(3,1,3);
%quantile(mse0,normcdf(-2:2))
fdx = min(isfinite(mse0-mse1),[],2);
nm = [mean(mse0(fdx,:)); mean(mse0(fdx,:)-mse1(fdx,:)); mean(mse1(fdx,:))];
plot(R,nm,'o');
title('Error decomposition');
legend('total','linear','higher order');

figure(2);
q = median(zerr(fdx_,:));
%q = mean(zerr(fdx_,:));
plot(R(2:end),q(2:end)./q(1:end-1))

for idx=1:length(R)
	namedfigure(3,'Roughness length over the cross section');
	subplot(length(R),1,idx);
	plot([zval(:,idx) zval(:,idx)*[1 1]+3*zerr(:,idx)*[-1 1]]);
	title(R(idx));
	ylim([-8 -2]);
	namedfigure(4,'Error decomposition in linear and higher order terms');
	subplot(length(R),1,idx);
	plot([mse0(:,idx), mse0(:,idx)-mse1(:,idx), mse1(:,idx)])
	title(R(idx));
	namedfigure(5,'Error autocorrlelation');
	subplot(length(R),1,idx);
	%plot(rho(:,idx))
	plot(squeeze(rho(:,idx,:)))
	title(R(idx));	
	namedfigure(6,'Standard error amplification due to reduced effective sample size');
	subplot(length(R),1,idx);
	plot(squeeze(f(:,idx,:)))
	title(R(idx));	
	figure(8);
	subplot(length(R),2,2*idx-1);
	n_ = sum(sum(dis(idx).bin.valid));
	plot((1:n_)/n_,sort(reshape(dis(idx).z0grid.err.Res(dis(idx).bin.valid),1,[])))
	title(R(idx));
	rmse(idx) = sqrt(mean(dis(idx).z0grid.err.Res(dis(idx).bin.valid).^2))
	subplot(length(R),2,2*idx);
	n_ = length(dis(1).z0grid.val.bias_ln_z0);
	plot((1:n_)/n_,sort(dis(idx).z0grid.val.bias_ln_z0));
	ylim(quantile(dis(idx).z0grid.val.bias_ln_z0,[0.05 0.95]));
	bias(idx) = nanmedian(dis(idx).z0grid.val.bias_ln_z0)
end
%figure(3);
%plot(sort(val(:,1) - val(:,end)))

figure(7);
q = sqrt(nanmedian(mse0));
plot(q);
set(gca,'xticklabel',R);
ylabel('Median of RMSE (m/s)');
xlabel('Radius of local estimate (m)');

figure(9);
imagesc(dis(1).z0grid.err.Res)
colorbar();



figure(1);
pdfprint('img/roughness-length-radius.eps');


