
	meta = sanggau_metadata();
	cc = colormap('lines');

	figure(1);
	clf();
	for idx=1:length(meta.filename.vadcp)	
		load(meta.filename.vadcp{idx});
		vadcp  = VADCP(vadcp);
		A{idx} = vadcp;
		plot(vadcp.utm.X,vadcp.utm.Y,'.','color',cc(idx,:));
		hold on	
	end
	axis equal

