% Fr 19. Jun 11:38:12 CEST 2015

	% maximum velocity vs. cs averaged velocity
	namedfigure(7,'Maximum vs. average velocity');
	meas.u_    = meas.u;
	meas.u     = filtfunc(meas.u,nf);
	[mval mdx] = max(u_mean);
	u_max      = cvec(meas.u(mdx,:));
	bar(meas.U,u_max);	
	c = u_max \ meas.U
	res = u_max*c - meas.U;
	res2 = res'*res/(length(res)-length(c))
	R2 = 1 - res2/var(meas.U)
	% affine - TODO should be the other way around, u_max is less certain then meas.U
	A    = [u_max.^0 u_max];
	c    = A \ meas.U
	res  = A*c - meas.U;
	res2 = res'*res/(length(res)-length(c))
	R2   = 1 - res2/var(meas.U)
