% 2014-08-19 15:18:01.502085849 +0200
% Karl Kastner, Berlin

load ../discharge/mat/2013-12-09-sanggau.mat
vel = discharge.bin.velocity;

% todo, this are second derivatives, should be devided by dt and dz

ddU = reshape( vel(2:end-1,2:end-1,2) ...
       -0.25*(   vel(2:end-1,1:end-2,2) ...
		+vel(2:end-1,3:end,2) ...
		+vel(1:end-2,2:end-1,2) ...
		+vel(3:end,2:end-1,2) ),[],1);
std_ddU = nanstd(ddU)

ddZ = reshape(discharge.ens.Z-mean(discharge.ens.Z,2)*[1 1 1 1],[],1);
std_ddZ = nanstd(ddZ)
edges = linspace(-1,1,50);
[hu] = hist(ddU,edges);
hu = hu/sum(hu);
[hz] = hist(ddZ,edges);
hz = hz/sum(hz);
%bar(edges(:),[hu(:),hz(:)]);
%plot(edges(:),[hu(:),hz(:)]);
[ax h1 h2] = plotyy(edges(:),hu(:),edges(:),hz(:));
% ylabel(ax(1),'m/s');
%ylabel(ax(2),'m');
 ylim(ax(1),[0 0.2]);
 ylim(ax(2),[0 0.2]);
 set(ax(1),'ytick',0:0.05:0.2);
  set(ax(2),'ytick',0:0.05:0.2)
legend('\epsilon_u (m/s)', '\epsilon_z (m)')

name = 'img/error_u_z.eps';
print(name,'-depsc');
system(['LD_LIBRARY_PATH= epstopdf ' name]);

