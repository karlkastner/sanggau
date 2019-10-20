% Di 14. Jul 08:32:05 CEST 2015
% Karl Kastner, Berlin

% stat.dir.mean		cross section direction
% stat.dir.std		
% stat.d.sn		delta of each transect in sn-coordinates
% stat.d.xy		delta of each transect in xy coordinate
% stat.d.std;		

function [X Y S N stat] = align_transect(X,Y,id)

% compute principal direction
for each transect individually, to compensate offset between transects
	- covariance
	- median of slopes
- warn if principal direction is not accurate (std of direction vector)
% transform to-n coordinates
% interpolate to constant step size
- get min-max, ngrid
% convolve pairs of transects
	delta a b is maximum / n
% mean delta of each transect
delta = mean(D);
- warn if delta is not accurate (std delta)
% transform delta back to XY
% shift transect coordinates

% test: alignment must be idempotent
% plot original transects and shifted transects
% 

end % function align_transect

