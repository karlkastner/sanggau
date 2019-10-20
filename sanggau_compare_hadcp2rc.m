% 2017 Feb 27 17:15

m  =   60;
n1 = 1124;
n2 = 4428;


if (~exist('reload','var') || reload)
	meta = sanggau_metadata();
	[hadcp hlevel offset] = sanggau_load_hadcp(meta);
	reload = 0;
end


% correction
v = -hadcp.velocity.earth(1:m,4424:10008,1);
v0 = v;
n=n1;
m1 = median(v(:,n-100:n),2);
m2 = median(v(:,n+1:n+100),2);
c = m1./m2;
v(:,n+1:end) = bsxfun(@times,c,v(:,n+1:end));

n=n2;
m1 = median(v(:,n-100:n),2);
m2 = median(v(:,n+2:n+102),2);
c = m1./m2;
v(:,n+2:end) = bsxfun(@times,c,v(:,n+2:end));
%m1 = median(v(:,n-100:n),2);
%m2 = median(v(:,n+2:n+102),2);
%c = m1./m2;


clf;
subplot(2,2,1)
plot(quantile(v',[0.1:0.1:0.9])');

subplot(2,2,3);
imagesc(v);
caxis([-1 1]) 


subplot(2,2,2)
%v  = quantile(v,[0.25 0.5 0.75])';
%v0 = quantile(v0,[0.25 0.5 0.75])';
v = mean(v(1:m,:))';
v0 = mean(v0(1:m,:))';
%v = quantile(-medfilt1(double(hadcp.velocity.earth(1:m,4424:10008,1))',48)',[0.25 0.5 0.75])';
plot(v);
hold on;
p = hadcp.pressure_bar(4424:10008);
plot(p);
pi = hadcp.pitch_rad(4424:10008);
plot(pi*10);
hold on
plot(hadcp.dat.FileNumber(4424:10008)/10);
vline(n1)
vline(n2)

subplot(2,2,4);
p  = meanfilt1(p(:),48);
v  = meanfilt1(v,48);
v0 = meanfilt1(v0,48);
plot3(p(:),v,(1:length(p))')
hold on
%plot3(p(:),v0,(1:length(p))')
view(0,90)
%plot(p(:),[v(:,2) v0(:,2)])
hold on
plot(p([n1 n2]),v([n1 n2]),'bo');


figure(2)
clf
subplot(2,2,1)
plot(quantile(hadcp.E(:,4424:10008,1),[0.25 0.5 0.75])');
hold on
plot(hadcp.dat.FileNumber(4424:10008)/10)
subplot(2,2,2);
% imagesc(double(max(hadcp.E(:,4424:10008,:),[],3)));
plot(quantile(hadcp.E(:,4424:10008,1)',[0.1:0.2:0.9])');

%-> mean estimate
%-> harmonic mean estimate
%-> lumped estimate
%(vpm)
%(corrected for transversal deviation)
nancorr([v p(:)])

