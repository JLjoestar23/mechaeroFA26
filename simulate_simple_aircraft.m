%% Simulate
system_params.Aw = 24.0; % m^2
system_params.kL = 0.1 * 180/pi; % 1/rad
system_params.kDO = 0.03; % unitless
system_params.kDI = 0.001 * (180/pi)^2; % 1/rad^2
system_params.rho = 1.225; % kg/m^3
system_params.m = 1134; % kg
system_params.g = 9.81; % m/s^2

T = 3400; % N
a = deg2rad(5); % deg -> rad
u = [T, a]; % thrust and AoA input

tspan = [0 300];
X0 = [0; 0; 61.34; 0];
options = odeset('AbsTol', 1e-9);
[t, X] = ode45(@(t, X) simple_aircraft(t, X, u, system_params), tspan, X0, options);

%% Plotting
plot_aircraft(t, X);
