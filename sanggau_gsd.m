% Sun Oct 26 14:34:51 CET 2014

addpath ../bottom/
bottom = load_gsd();

r = 182:185;

ln_z_0 = log10(0.1*2.^bottom.d84(r))
	
