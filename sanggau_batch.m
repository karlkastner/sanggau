% 2015-08-06 14:40:01.596742439 +0200
% 2015-01-09 13:24:05.587222894 +0100

% Karl Kastner, Berlin

% batch script to produce all values and plots used in the thesis

% TODO pflag/reload
% TODO 1d bathy
% todo 2d bathy
% todo vert profile plot
% todo trans profile plot


	% estimate discharge
	sanggau_process_discharge();
	% transversal profile parameter
	sanggau_transversal_profile_parameter();
	% vertical profile parameter
	sanggau_vertical_profile_parameter();
	% transversal-vertical profile error cross correlation
	sanggau_error_correlation()
	% expected error of HADCP discharge estimate
	sanggau_sdm_scale_vs_depth
	% export data to tex file
	sanggau_calib_table()

sanggau_plot_z_0

	jflag = 0;

	fprintf(1,'Starting calibration\n');

	method = {'empirical', 'empirical','theoretical'};
	worder = [0,1,1];
	
	for idx=3:length(method)
		sanggau_hadcp_calibration(method{idx},worder(idx),jflag);
	end
	
