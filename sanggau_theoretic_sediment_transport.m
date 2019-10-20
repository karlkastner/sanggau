% Tue  9 May 15:54:35 CEST 2017
% Karl Kastner, Berlin

if (~exist('pflag','var'))
	pflag = 0;
end
fflag = pflag;

ps = 3;

Q = linspace(0,1e4,10)';

site(1).name = 'Sanggau';
site(1).W    = 600;
site(1).C    = 50;
site(1).S    = 5e-5;
site(1).d    = 1; % 1mm
site(1).H    = (Q.^2/(site(1).W^2*site(1).C^2*site(1).S)).^(1/3);
name         = {'Sanggau','Rasau'};

site(2).name = 'Rasau';
site(2).W    = 450;
site(2).C    = 71;
site(2).d    = 0.2; % 0.2mm
site(2).H    = 15*ones(size(Q)); % constant depth


for sdx=1:length(site)

W = site(sdx).W;
H = site(sdx).H;
d = site(sdx).d;
C = site(sdx).C;

d50 = d;
sd  = 1;
d90 = d50*(1 + 1.28*sd);
%exp(1.28*sd)
U   = Q./(W.*H);

Qs = [];
Phi = [];
[Qs(:,1) void Phi(:,1)] = total_transport_bagnold(C,d,U,H,W);
[Qs(:,2) void Phi(:,2)] = total_transport_engelund_hansen(C,d,U,H,W);
[Qs(:,3) void Phi(:,3)] = bed_load_transport_rijn(C,d50,d90,U,H,W);
[Qs(:,4) void Phi(:,4)] = suspended_transport_rijn(C,d50,d90,sd,U,H,W);
[Qs(:,5) void Phi(:,5)] = total_transport_rijn(C,d50,d90,sd,U,H,W);
[Qs(:,6) void Phi(:,6)] = bed_load_engelund_fredsoe(C,d,U,H,W);
[Qs(:,7) void Phi(:,7)] = bed_load_einstein(C,d,U,H,W);

Cs = bsxfun(@times,Qs,1./Q);


%namedfigure(sdx,name{sdx});
splitfigure([2 2],[sdx 1],fflag,name{sdx});
cla();
subplot(2,2,1);
plot(Q,Qs);
ylabel('Q_s (kg/s)');
xlabel('Q_w (m^3/s)');

splitfigure([2 2],[sdx 2],fflag,name{sdx});
plot(Q,Phi);
ylabel('Phi');
xlabel('Q_w (m^3/s)');

splitfigure([2 2],[sdx 3],fflag,name{sdx});
% kg/m^3 to mg/l
% C = C*(kilo/milli)*(m^3/l)
Cs = Cs*(1e3/1e-3)*(10^-3)
plot(Q,Cs);
ylabel('C_s (mg/l) (apparrent, including bed-load)');
xlabel('Q_w (m^3/s)');

splitfigure([2 2],[sdx 4],fflag,name{sdx});
plot(NaN(10,1),NaN(10));
%ox off
axis off
legend('location','northwest','Bagnold total','Engelund-Hansen total','rijn bed','rijn susp','rijn total','bed engelund-fredsoe','bed einstein');

end

% sanngau to sea
L = 3e5;
Q = 5e3;
W = site(1).W;
C = site(1).C;
S = site(1).S;
Hin = (Q.^2/(W^2*C^2*S)).^(1/3);
U = Q/(Hin*W);
Qs = total_transport_engelund_hansen(C,d,U,Hin,W);
H = 15;

% volumetric sediment transport
Qv = Qs/2650;
% volume of kapuas at high flow
V = H*W*L;
% morphological time scale in seconds
T = V/Qv
% in years
T = T/(86400*365.25);
fprintf('Morphological time scale (years): %g\n',T);

if (pflag)
	pdfprint(11,'img/sanggau-sediment-transport-formulae',ps);
	pdfprint(14,'img/legend-sediment-transport-formulae',ps);
	pdfprint(21,'img/rasau-sediment-transport-formulae',ps);
	pflag = 0;
end

