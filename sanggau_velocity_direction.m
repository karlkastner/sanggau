% 2018-01-22 19:53:59.766138035 +0100

load([ROOTFOLDER,'dat/kapuas/2013-12-09-sanggau-transect/mat/vadcp.mat']); v=VADCP(vadcp)
v.depth_average_velocity('earth')
 hist(v.ens.velocity.earth,linspace(-1,1,10))

vel = v.ens.velocity.earth;
fdx=rms(vel(:,1:2),2)>0.2;
hist(vel(fdx,1:2),linspace(-1,1,10))
fdx=1:50:v.nens; quiver(v.X(fdx),v.Y(fdx),v.ens.velocity.earth(fdx,1),v.ens.velocity.earth(fdx,2),10)

