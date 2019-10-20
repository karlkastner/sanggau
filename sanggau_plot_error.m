% 2017-07-22 01:57:58.505944351 +0800

figure(1)
clf;
rmse = [];
Q = [];
field = 'zb';
% TODO integrate ! multiply by width
for idx=1:5;
	subplot(2,2,1);
	plot(-sqrt(calib.cs_(idx).gridNr.msei.(field))./calib.cs_(idx).gridNr.vali.(field));
	ylabel('rmse_q (m^3/s)');

	subplot(2,2,2)
	plot(-sqrt(calib.cs_(idx).gridNr.msei.(field))./calib.cs_(idx).gridNr.vali.(field));
	hold on;
	ylabel('rmse_q/q');
	rmse(idx,1) = sqrt(sum(calib.cs_(idx).gridNr.msei.(field)))
	Q(idx,1) = sum(calib.cs_(idx).gridNr.vali.(field));
end
rmse
Q
rmse_rel = rmse./Q
