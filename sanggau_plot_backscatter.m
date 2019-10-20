% 2015-09-27 14:15:13.236423636 +0200
% Karl Kastner, Berlin

	hadcpname = [ROOTFOLDER,'src/hadcp/mat/hadcp-sanggau-squeezed-1800-preloaded.mat'];
	load(hadcpname);

	E = hadcp.E;
	dE = [abs(E(10,:,1)-E(10,:,2)); abs(E(10,:,1)-E(10,:,3)); abs(E(10,:,2)-E(10,:,3))]';
	legend(num2str([12,13,23]'));

	namedfigure(1,'Backscatter differences at bin 10');
	clf();
	plot(hadcp.time,dE,'-');
	legend(num2str([12,13,23]'));
	datetick('x','dd/mm/yy');

	% find samples where beam was obstructed
	% assumes two beams were not obstructed at least half of the time
%	q = quantile(dE,0.5);
%	fdx = bsxfun(@gt, dE, 2*q);

	figure(10)
	plot(hadcp.time,squeeze(E(10,:,:)))

	figure(2);
	clf();
	plot(sort(dE));
	
	fdx = (dE > 10);
	figure(3);
	clf;
	fdx = max(fdx,[],2);
	plot(fdx,'.');
	

	figure(1);
	pdfprint('img/sanggau-hadcp-backscatter-differences.pdf');
	
	

