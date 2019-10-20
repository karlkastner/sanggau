% Sat 20 May 13:30:08 CEST 2017
function [d hbar h bottom fdx] = sanggau_load_gsd()
	load([ROOTFOLDER,'/src/bed-material/mat/gsd.mat']);

	% select bed material samples at sanggau
	R   = 1500;
	xy0 = [454268    10012714];
	fdx = (bottom.X-xy0(1)).^2 + (bottom.Y-xy0(2)).^2 < R^2;

	d = cvec(2.^bottom.histogram.centre);
	h = bottom.histogram.h(fdx,:);
	hbar = mean(h);
end

