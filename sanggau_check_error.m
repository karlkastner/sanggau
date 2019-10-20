% 2015-07-30 19:00:01.183374865 +0200

% error of ln_z0
	for idx=1:size(a,1);
		for jdx=1:size(a,2);
			err.U(idx,jdx) = nanmedian(a(idx,jdx).gridN.err.U);
			err.ln_z0(idx,jdx) = nanmedian(a(idx,jdx).gridN.err.ln_z0);
		end;
	end
	e,
subplot(2,2,1)
plot(err.U')
subplot(2,2,2)
plot(err.ln_z0');

subplot(2,2,3)
u = discharge.pseudo.gridN.U;
res = bsxfun(@minus,u,nanmean(u,2));
a = acf_man(res(50:end-50,:));
s=var(a(50:end-50,:));
as = bsxfun(@times,a,s);
errorbar(sign(mean(as')).*sqrt(abs(mean(as'))),sqrt(serr(as')),'.')

