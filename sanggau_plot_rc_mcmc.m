% 2015-02-23 21:43:08.716044242 +0100
% Karl Kastner, Berlin

load('mat/sanggau-mcmc.mat');
%load('mat/sanggau-mcmc-2015-02-23-20-08-40.mat');
modelfun = @(h,c) rc.rcfunc(c',[],h);

namedfigure(2,'Parameter chain');
clf();
mcmcplot(chain,[],res,'chainpanel');

namedfigure(30,'Parameter distributions');
p = 0.05;
nh = 100;
for idx=1:size(chain,2)
	q = quantile(chain(:,idx),[p 1-p]);
	subplot(2,2,idx);
	x = linspace(q(1),q(2),100);
	hist(chain(:,idx),x);
	title(idx);
	xlim([q(1) q(2)])
end

figure(3); clf
mcmcplot(chain,[],res,'pairs');
%ylim([0 350]);
%xlim([0 0.4]);

% TODO delaunay mesh the parameter pairs and plot the patch

%%
figure(4);
clf();
mcmcplot(sqrt(s2chain),[],[],'hist')
title('Error std posterior')

%%
% A point estimate of the model parameters can be calculated from the
% mean of the |chain|. Here we plot the fitted model using the
% posterior means of the parameters.

contour_;

figure(1)
clf();
plot(data.xdata,data.ydata,'s');
%xlim([0 400]);
%xlabel('x [mg/L COD]'); ylabel('y [1/h]');
hold on
%plot(x,modelfun(x,mean(chain)),'-k')
x = out.data{1};
plot(x,modelfun(x,mean(chain)),'-k')
plot(x,modelfun(x,median(chain)),'--k')
hold off
legend('data','mean','median',0)

%%
% Instead of just a point estimate of the fit, we should also study
% the predictive posterior distribution of the model. The |mcmcpred|
% and |mcmcpredplot| functions can be used for this purpose. By them
% we can calculate the model fit for a randomly selected subset of the
% chain and calculate the predictive envelope of the model. The grey
% areas in the plot correspond to 50%, 90%, 95%, and 99% posterior
% regions.

figure(5); clf
mcmcpredplot(out);
hold on
plot(data.xdata,data.ydata,'s'); % add data points to the plot
xlabel('x [mg/L COD]'); ylabel('y [1/h]');
hold off
title('Predictive envelopes of the model')

figure(6)
figure(1e3);
clf
r = zeros(size(chain,1),1);
for idx=1:4;
        r = r+(modelfun(rc.h0(idx),chain)'-rc.q0(idx)).^2;
end
r = sqrt(r);
N = nchoosek(1:3,2);
for idx=1:size(N,1)
        [c id] = unique(chain(:,N(idx,:)),'rows');
        r_=(r(id));
        subplot(2,2,idx)
        T = delaunay(c(:,1),c(:,2));
        patch([c(T(:,1),1) c(T(:,2),1) c(T(:,3),1)]',[c(T(:,1),2) c(T(:,2),2) c(T(:,3),2)]',[r_(T(:,1)) r_(T(:,2)) r_(T(:,3))]','edgecolor','none')
        caxis([0 quantile(r,0.95)]);
	colormap gray
end

