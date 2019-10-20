% 2015-07-30 15:49:25.160041146 +0200

f=dir('~/phd/src/discharge/mat/*sang*dw-1*1d*dis*');
for idx=1:length(f);
	load(['~/phd/src/discharge/mat/' f(idx).name])
	a(idx)=discharge;
end

f=dir('~/phd/src/discharge/mat/2015-07-30/*sang*dw-1*dis*');
for idx=1:length(f);
	load(['~/phd/src/discharge/mat/2015-07-30/' f(idx).name])
	b(idx)=discharge;
end
for idx=1:length(a)
	d(idx,:) = [a(idx).cs.discharge.total(1) b(idx).cs.discharge.total(1)];
end;
d
100*diff(d,[],2)./sum(d,2)


figure(1);
clf();
for idx=1:length(a)
	subplot(2,3,idx);
	%plot([a(1).gridN.val.U b(1).gridN.val.U])
	U = [nanmean(b(idx).gridNZ.val.U,2) a(idx).gridN.val.U];
	plot(U);
end

figure(2)
clf()
for idx=1:length(a)
	subplot(2,3,idx)
	bin = a(idx).adcp.bin;
	ens = a(idx).adcp.ens;
	vel = bin.velocity;
	%vel = ens.velocity;
	surface(bin.time,bin.Z,vel.cs(:,:,1),'edgecolor','none');
	id     = sub2ind([bin.n ens.n],ens.last,(1:ens.n)');
	Z = bin.Z;
	hold on
	plot(ens.time,Z(id),'k');
%	last = ens.last;
%	time = ens.time;
%	for jdx=1:ens.n
%		plot(time(jdx),Z(last(jdx),jdx),'b')
%	end
end

