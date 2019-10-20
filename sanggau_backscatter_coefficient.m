% Tue 16 May 16:34:32 CEST 2017
% theoretical backscatter coefficient for Sanggau based on bed material grain size

% note: this is not in the right ball park,
%       the computed sediment concentration is about a factor 10 too low

if (~exist('reload','var') || reload)
	[d hb h bottom fdx] = sanggau_load_gsd();
	reload = 0;
end

figure(1);
clf();
plot(bottom.X,bottom.Y,'.');
hold on;
plot(bottom.X(fdx),bottom.Y(fdx),'.');

% compute suspended gsd

h_s = suspended_grain_size(d,h)
d'
'h mean of unsuspended'
%hb_ = mean(hb)
hb*d
'd mean of suspended'
d_mu = (h_s*d)'
mean(d_mu)
median(d_mu)

% compute ks
f = 1.2e6;
1e-3*d
ks2_d = backscatter_coefficient(1e-3*d,f);
%ks2 = backscatter_coefficient(1e-3*d,f);

% histrogram of ks2 (not normalised)
hks2  = bsxfun(@times,h,rvec(ks2_d));

% ks^2 for each sample
ks2    = sum(hks2,2);
ks2_mu = mean(ks2)
ks2_sd = std(ks2)
% from mean grain size (mass averaged)
ks2       = backscatter_coefficient(1e-3*d_mu,f);
ks2_mu(2) = mean(ks2)
ks2_sd(2) = std(ks2)

% normalise
hks2 = mean(hks2)';
hk   = hks2/sum(hks2);

hb_mu = (mean(hb))';
h_mu  = mean(h)';

figure(2)
subplot(2,2,1)
%plot(log2(d),[hb_mu h_mu hk])
subplot(2,2,2)
%plot(log2(d),ks2_d);

figure(3);


