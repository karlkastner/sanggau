% 2014-08-15 19:50:45
% Karl Kastner, Berlin

% instrument depth+redeployment depth
offset = [5.3058    0.5296    2.8458    4.5072] ...
	 - 2.075*[0 1 1 1];
oname_C = { ...
'2013-12-09-sanggau', ...
'2014-02-20-sanggau', ...
'2014-04-18-sanggau', ...
'2014-06-18-sanggau' };
meta = sanggau_metadata();
% TODO make FDM, IDW and RBM strings
%bmethod = [Discharge.FDM, Discharge.IDW, Discharge.RBM];
for idx=1:4
	iname = meta.filename.vadcp{idx};
	oname = oname_C{idx};
	load(iname);
        coord = load([ROOTFOLDER,'/dat/kapuas/metadata-protocols/2013_12_06_sanggau_coordinates.csv']);
        arg   = {'xpmin',coord(1,1),'ypmin',coord(1,2),'xpmax',coord(2,1),'ypmax',coord(2,2)};
	
	tmode   = Discharge.NO_TIME;
	cmethod = 0; % point orientation provided
	hmethod = 0;
	vmethod = Discharge.FDM;
	dw = 1;
	dz = 5;
	dt = [];
	for jdx=1:length(bmethod)
	% TODO this does unnecessarily perform the discharge computation
	discharge = Discharge(adcp, arg{:}, ...
		'cmethod', cmethod, ...
		'bmethod', bmethod(jdx), ...
		'hmethod', hmethod, ...
		'vmethod', vmethod, ...
		'tmode', tmode, ...
		'dt', dt, ...
		'dw', dw, ...
		'dz', dz);
		discharge.calc_discharge();
		bottom(jdx).Z(:,idx) = discharge.bgrid.val(:,1);
		if (jdx==3)
			bottom(jdx).err(:,idx) = discharge.bgrid.val(:,2);
			bottom(jdx).R(:,idx)   = discharge.bgrid.val(:,3);
		end

		% align with respect to the first calibration
		N = discharge.bgrid.cX1*discharge.width;
		bottom(jdx).offset(idx) = nanmedian(bottom(jdx).Z(:,idx) - bottom(jdx).Z(:,1));
		bottom(jdx).Z_(:,idx) = bottom(jdx).Z(:,idx) - offset(idx);
		%bottom(jdx).Z_(:,idx) = bottom(jdx).Z(:,idx) - bottom(jdx).offset(idx);
	end % jdx
end % idx

% error between FDM and IDW
q = quantile(bottom(1).Z-bottom(2).Z,[0.1587 0.5 0.8413]);
bias = q(2,:);
err  = 0.5*(q(3,:)-q(1,:));
namedfigure(1,'FDM vs IDW');
for idx=1:3
	subplot(3,1,idx);
	plot(-bottom(idx).Z_);
end

namedfigure(2,'Individual campaigns');
clf();
for idx=1:4
	for jdx=1:length(bottom)
		b(:,jdx) = bottom(jdx).Z_(:,idx);
	end
	subplot(4,1,idx);
	plot(-b);
end

namedfigure(3,'Bottom profile');
clf();
plot(N,-bottom(3).Z_);
ylabel('m');
xlabel('m');
s0 = 10;
s  = 150;
r  = s*tan(deg2rad(1.3/2));
dd = -2.075;
hold on
patch([10 s0+s s0+s],[0 -r r],[1 1 1],'facecolor','none');
patch([10 s0+s s0+s],dd+[0 -r r],[1 1 1],'facecolor','none');
% plot water level at individual campaigns
plot([0 0 0 0; discharge.width*[1 1 1 1]],[1 1]'*offset);
xlim([0 discharge.width]);

%hold on;plot([10 160],[-2.075 -2.075],'k-');
%hold on;plot([10 160],[0 0],'k-');

if (pflag)
	pdfprint('img/sanggau-bottom-1d');
end

namedfigure(4,'range and error')
subplot(2,1,1)
hist(bottom(3).R)
subplot(2,1,2)
hist(bottom(3).err)

