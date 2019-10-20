% 2015-02-26 17:43:27.749544127 +0100

if (~exist('reload','var') || reload)
	load('mat/hadcp-sanggau.mat'); 
	hadcp=HADCP(hadcp) ;
	relaod = 0;
end
figure(1);
clf();
for idx=1:3
	subplot(3,3,idx)
	imagesc(hadcp.vel.beam(:,:,idx)); caxis([-1 1]); colorbar
%	subplot(3,2,2)
%	imagesc(hadcp.vel.earth(:,:,2)); caxis([-1 1]); colorbar
	subplot(3,3,3+idx)
	imagesc(hadcp.E(:,:,idx));  colorbar
end
subplot(3,3,7)
plot(hadcp.time,[cvec(hadcp.heading_rad) cvec(hadcp.pitch_rad) cvec(hadcp.roll_rad)]);
datetick
subplot(3,3,8)
plot(hadcp.time,hadcp.power_W);
datetick

