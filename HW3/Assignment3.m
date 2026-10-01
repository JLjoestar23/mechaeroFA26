%% Problem 1

m = 4000; % kg
r = 0.8; % m
k1 = 20000; % N/m
k2 = 20000; % N/m
c1 = 2000; % Ns/m
c2 = 2000; % Ns/m
L1 = 0.9; % m
L2 = 1.4; % m

M = [m  0; 
     0  m*(r^2)];

C = [c1+c2          -c1*L1+c2*L2; 
     -c1*L1+c2*L2   c1*L1^2+c2*L2^2];

K = [k1+k2          -k1*L1+k2*L2; 
     -k1*L1+k2*L2   k1*L1^2+k2*L2^2];

% solve GEVP
[U_mode, wn_sq] = eig(K, M);

U1 = U_mode(:,1);
U2 = U_mode(:,2);

U1 = U1 / U1(1);
U2 = U2 / U2(1);

wn1 = sqrt(wn_sq(1,1));
wn2 = sqrt(wn_sq(2,2));

%% a)

C_diag = U_mode'*C*U_mode;

z1 = C_diag(1,1) / (2*wn1);
z2 = C_diag(2,2) / (2*wn2);

%% b)

M_sqrt = sqrt(M);

Mt = M_sqrt\M/M_sqrt;
Kt = M_sqrt\K/M_sqrt;
Ct = M_sqrt\C/M_sqrt;

[P, ~] = eig(Kt);

K_diag = P'*Kt*P;

C_diag = P'*Ct*P;

M_diag = P'*Mt*P;

%% c) solve IVP for analytical solution

x0 = [0.15; 0];

q0 = U_mode'*x0;

a = U_mode\x0;

a1 = a(1);
a2 = a(2);

t_a = linspace(0,15,500);
s1 = z1*wn1;
s2 = z2*wn2;
wd1 = wn1*sqrt(1-z1^2);
wd2 = wn2*sqrt(1-z2^2);

X_a = a1*U_mode(:,1)*exp(-s1*t).*cos(wd1*t) ... 
      + a2*U_mode(:,2)*exp(-s2*t).*cos(wd2*t);
% X_a = a1*U_mode(:,1)*(exp(-s1*t).*(cos(wd1*t) + (s1/wd1)*sin(wd1*t))) ...
%     + a2*U_mode(:,2)*(exp(-s2*t).*(cos(wd2*t) + (s2/wd2)*sin(wd2*t)));

%% solve numerically

tspan = [0 15];

function dXdt = sus_undamped(t, X)
    
    m = 4000; % kg
    r = 0.8; % m
    k1 = 20000; % N/m
    k2 = 20000; % N/m
    c1 = 2000; % Ns/m
    c2 = 2000; % Ns/m
    L1 = 0.9; % m
    L2 = 1.4; % m
    
    M = [m  0; 
         0  m*(r^2)];
    C = [c1+c2          -c1*L1+c2*L2; 
         -c1*L1+c2*L2   c1*L1^2+c2*L2^2];
    K = [k1+k2          -k1*L1+k2*L2; 
         -k1*L1+k2*L2   k1*L1^2+k2*L2^2];

    q = X(1:2);
    dq = X(3:4);
    
    ddq = (-M\C)*dq + (-M\K)*q;

    dXdt = [dq; ddq];
end

[t_n, X_n] = ode45(@(t, X) sus_undamped(t, X), tspan, [x0; 0; 0]);

% plotting
navy   = [0.00 0.20 0.45];
orange = [0.85 0.40 0.05];
figure('Color','w');
sgtitle('Comparing Analytical and Numerical Solutions')

ax1 = subplot(2,1,1);
h1 = plot(t_a, X_a(1,:), 'Color', navy, 'LineWidth', 2);
hold on;
h2 = plot(t_n, X_n(:,1), '.', 'Color', orange, 'MarkerSize', 8);
grid on;
xlabel('Time (s)');
ylabel('Vertical Displacement');
ylim([min(X_a(1,:))*1.2, max(X_a(1,:))*1.2]);
legend([h1, h2], {'Analytical', 'Numerical'});
hold off;

ax2 = subplot(2,1,2);
h3 = plot(t_a, X_a(2,:), 'Color', navy, 'LineWidth', 2);
hold on;
h4 = plot(t_n, X_n(:,2), '.', 'Color', orange, 'MarkerSize', 10);
grid on;
xlabel('Time (s)');
ylabel('Angular Displacement');
ylim([min(X_a(2,:))*1.2, max(X_a(2,:))*1.2]);
hold off;
%% Problem 2

% a)

Aw = 24.0; % m^2
kL = 0.1 * 180/pi; % 1/rad
kDO = 0.03; % unitless
kDI = 0.001 * (180/pi)^2; % 1/rad^2
rho = 1.225; % kg/m^3
m = 1134; % kg
g = 9.81; % m/s^2

% solve for level flight inputs at 2 deg AoA
a_des = deg2rad(2);
Cl = kL*a_des;
Cd = kDO + kDI*a_des^2;

A1 = [cos(a_des) -0.5*rho*Cd*Aw;
     sin(a_des) 0.5*rho*Cl*Aw];
b1 = [0; m*g];

x1 = A1\b1;

% solve for level flight inputs at 3 deg AoA
a_des = deg2rad(3);
Cl = kL*a_des;
Cd = kDO + kDI*a_des^2;

A2 = [cos(a_des) -0.5*rho*Cd*Aw;
     sin(a_des) 0.5*rho*Cl*Aw];
b2 = [0; m*g];

x2 = A2\b2;

%% b, c)

system_params.Aw = 24.0; % m^2
system_params.kL = 0.1 * 180/pi; % 1/rad
system_params.kDO = 0.03; % unitless
system_params.kDI = 0.001 * (180/pi)^2; % 1/rad^2
system_params.rho = 1.225; % kg/m^3
system_params.m = 1134; % kg
system_params.g = 9.81; % m/s^2

T = 3900; % N
a = deg2rad(3); % deg -> rad
u = [T, a]; % thrust and AoA input

tspan = [0 600];
X0 = [0; 0; 0; 0];
options = odeset('AbsTol', 1e-9);
[t, X] = ode45(@(t, X) simple_aircraft_v2(t, X, u, system_params), tspan, X0, options);

%% plotting
plot_aircraft(t, X);

%% d)

system_params.Aw = 24.0; % m^2
system_params.kL = 0.1 * 180/pi; % 1/rad
system_params.kDO = 0.03; % unitless
system_params.kDI = 0.001 * (180/pi)^2; % 1/rad^2
system_params.rho = 1.225; % kg/m^3
system_params.m = 1134; % kg
system_params.g = 9.81; % m/s^2

T = 3900; % N
a = deg2rad(3); % deg -> rad
u = [T, a]; % thrust and AoA input

tspan = [0 600];
X0 = [0; 0; 0; 0];
options = odeset('AbsTol', 1e-9);
[t, X] = ode45(@(t, X) simple_aircraft_v3(t, X, u, system_params), tspan, X0, options);

%% plotting
plot_aircraft(t, X);