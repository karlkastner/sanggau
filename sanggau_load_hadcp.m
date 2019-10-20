% Sun Nov  2 11:31:42 CET 2014
% Karl Kastner, Berlin

function [hadcp hlevel offset] = sanggau_load_hadcp(meta, reload);
	% squeezed has fast ping samples averaged over 30min
	filename = meta.filename.hadcp;

	if (nargin() > 2 & ~isempty(meta.hadcp.squeeze.dt) & (meta.squeeze.dt > 0))
		pname    = [filename(1:end-4),'-squeezed-',num2str(meta.hadcp.squeeze.dt),'-preloaded.mat'];
		rawname  = [filename(1:end-4) '-squeezed-',num2str(meta.hadcp.squeeze.dt),'.mat'];
	else
		pname     = [filename(1:end-4) '-preloaded.mat'];
		rawname   = filename;
	end
	if (exist(pname) ...
		&& ~system(['test ',pname,' -nt ',rawname]) ... % somehow test returns 0 when second file is newer
		&& ~(nargin() > 1 && ~isempty(reload) && reload))
		load(pname);
	else
		% load hadcp data
		load(rawname);

		% preprocess hadcp data
		hadcp = HADCP(hadcp,'X',meta.hadcp.x0,'Y',meta.hadcp.y0);
		%hadcp = HADCP(hadcp,'X',meta.utm.x,'Y',meta.utm.y);

		% convert to earth velocities
		hadcp.to_earth();

		% correct level for redeployment
		fdx         = hadcp.time > meta.hadcp.redeploy.t;
		hlevel      = hadcp.idepth_m;
		hlevel(fdx) = hlevel(fdx) - meta.hadcp.redeploy.d;
	%	fdx         = find(hadcp.time > meta.hadcp.redeploy.t);
		offset      = zeros(size(hlevel));
		offset(fdx) = meta.hadcp.redeploy.d;

		save(pname,'hadcp','offset','hlevel');
	end
	% project velocity to the cross section
	hadcp.velocity.cs = rotvel(meta.dir, hadcp.velocity.earth);
	% rotate hadcp velocity onto the cross section	
	[hadcp.N, hadcp.T]         = xy2nt(hadcp.bin.X, hadcp.bin.Y, meta.centre(1), meta.centre(2), meta.dir);
        [hadcp.bin.N, hadcp.bin.T] = xy2nt(hadcp.bin.X, hadcp.bin.Y, meta.centre(1), meta.centre(2), meta.dir);
end % sanggau_load_hadcp

