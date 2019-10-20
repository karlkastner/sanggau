
% calculate ... diameter
d50    = [1e-3 3e-2]; % TODO exact
S0     = 1e-5; % TODO exact
%rhos  = 22000;
theta  = Constant.critical_shear_stress_ratio(dstar,'bronwlie-1981')
theta  = Constant.critical_shear_stress_ratio(dstar,'bronwlie-bonneville')
theta  = Constant.critical_shear_stress_ratio(dstar,'soulsby-1997')

% shear stress
tau0 = Constant.g*d50.*theta
	% Yalin and Karahan, 1979, Zeller (1963) for re > 40
	%tau0  = Constant.g*(Constant.density.quartz - Constant.density.water)*d50
us2  = tau0./Constant.density.water
H    = us2./(Constant.g*S0)

%H = (rho-rho)/rho D50 tau / S0

