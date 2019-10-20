% Sun Dec 14 14:08:59 CET 2014
% estimationg of standard error of the depth and flow velocity

if (0)
meta  = sanggau_metadata();
load([ROOTFOLDER,'/src/water-level/mat/water-level.mat']);
level.val  = K(nid.sanggau_merged).slot.depth;
level.time = K(nid.sanggau_merged).slot.time;

calib = load_vadcp_discharge(meta.filename_C,level);
%vadcp = sanggau_load_vadcp(meta, level);
hadcp = sanggau_load_hadcp(meta);
end
subplot(2,2,1)
p=repmat(hadcp.pingperens,hadcp.nbins(1),1);
v = hadcp.velocity.earth(:,:,4);
%hist(abs(v(:)),1000);
%v=v(:);
fdx =(v~=0 & isfinite(v) & p == 600);
v=v(fdx);
n = length(v(:));
plot(sort(abs(v(:))),(1:n)/n);
xlim([0 1]);
nanmedian(abs(v))
nanmean(abs(v))

s = std(vadcp.bottom_A');
fdx = find(isfinite(s));
n = length(fdx);
subplot(2,2,3)
plot(sort(s(fdx)),(1:n)/n);
nanmean(s(fdx))

%nom = s(1:end-1).*s(2:end);
%den = s(1:end-1).*s(1:end-1);
subplot(2,2,4)
autocorr(s(fdx),100)

figure(2)
plot([hadcp.pitch_rad(:)*100 hadcp.idepth_m(:)])


% cdf of depth standard deviation (not standard error, does not average out)


