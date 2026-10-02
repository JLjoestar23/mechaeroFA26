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

% U1 = U1 / U1(1);
% U2 = U2 / U2(1);

wn1 = sqrt(wn_sq(1,1));
wn2 = sqrt(wn_sq(2,2));

%% a)

C_diag = U_mode'*C*U_mode;

z1 = C_diag(1,1) / (2*wn1);
z2 = C_diag(2,2) / (2*wn2);

%% c)

M_sqrt = sqrt(M);

Mt = M_sqrt\M/M_sqrt;
Kt = M_sqrt\K/M_sqrt;
Ct = M_sqrt\C/M_sqrt;

[P, ~] = eig(Kt);

K_diag = P'*Kt*P;

C_diag = P'*Ct*P;

M_diag = P'*Mt*P;

%% d) solve IVP for analytical solution

s1 = z1*wn1;
s2 = z2*wn2;
wd1 = wn1*sqrt(1-z1^2);
wd2 = wn2*sqrt(1-z2^2);

x0 = [0.15; 0];

q0 = U_mode\x0;

phi1 = atan2(-s1, wd1);
phi2 = atan2(-s2, wd2);
a1 = q0(1)/cos(phi1);
a2 = q0(2)/cos(phi2);

t_a = linspace(0, 20, 500);

X_a = a1*U_mode(:,1)*exp(-s1*t_a).*cos(wd1*t_a + phi1) ... 
      + a2*U_mode(:,2)*exp(-s2*t_a).*cos(wd2*t_a + phi2);

%% e) solve numerically

tspan = [0 20];

function dXdt = sus_damped(t, X)
    
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

[t_n, X_n] = ode45(@(t, X) sus_damped(t, X), tspan, [x0; 0; 0]);

% plotting
navy   = [0.00 0.20 0.45];
orange = [0.85 0.40 0.05];
figure('Color','w');
sgtitle('Comparing Analytical and Numerical Solutions')

subplot(2,1,1);
plot(t_a, X_a(1,:), 'Color', navy, 'LineWidth', 2, 'DisplayName', 'Analytical');
hold on;
plot(t_n, X_n(:,1), '--', 'Color', orange, 'LineWidth', 2, 'DisplayName', 'Numerical');
grid on;
xlabel('Time (s)');
ylabel('Vertical Displacement');
ylim([min(X_n(:,1))*1.2, max(X_n(:,1))*1.2]);
legend();
hold off;

subplot(2,1,2);
plot(t_a, X_a(2,:), 'Color', navy, 'LineWidth', 2);
hold on;
plot(t_n, X_n(:,2), '--', 'Color', orange, 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Angular Displacement');
ylim([min(X_n(:,2))*1.2, max(X_n(:,2))*1.2]);
hold off;

%% f) solve with forcing input

% spatial distribution of the input in physical coordinates
F_dir = [1; 0];

% spatial distirbution of the input in modal coordinates
Q_dir = U_mode'*F_dir;

% input frequency
Omega = 3.5; % rad/s

% transfer function for each modal coordinate
Q_tfs = Q_dir ./ (-Omega^2 + diag(C_diag)*1i*Omega + diag(K_diag));

% transform back into physical coordinates
X_p = U_mode*Q_tfs;

% get magnitude and phase
X_mags = abs(X_p);
X_phase = angle(X_p);

% particular solutions in physical coordinates
X_p = X_mags*500.*sin(Omega*t_a + X_phase);

%% e) solve numerically

tspan = [0 20];

function dXdt = sus_damped_harmonic(t, X)
    
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

    ddq = M\([500; 0]*sin(3.5*t) - C*dq - K*q);

    dXdt = [dq; ddq];
end

[t_n, X_n] = ode45(@(t, X) sus_damped_harmonic(t, X), tspan, [0; 0; 0; 0]);

% plotting
navy   = [0.00 0.20 0.45];
orange = [0.85 0.40 0.05];
figure('Color','w');
sgtitle('Comparing Analytical and Numerical Solutions')

subplot(2,1,1);
plot(t_a, X_p(1,:), 'Color', navy, 'LineWidth', 2, 'DisplayName', 'Analytical Steady State');
hold on;
plot(t_n, X_n(:,1), '--', 'Color', orange, 'LineWidth', 2, 'DisplayName', 'Numerical');
grid on;
xlabel('Time (s)');
ylabel('Vertical Displacement');
ylim([min(X_n(:,1))*2, max(X_n(:,1))*2]);
legend();
hold off;

subplot(2,1,2);
plot(t_a, X_p(2,:), 'Color', navy, 'LineWidth', 2);
hold on;
plot(t_n, X_n(:,2), '--', 'Color', orange, 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Angular Displacement');
ylim([min(X_n(:,2))*2, max(X_n(:,2))*2]);
hold off;

%% Problem 2

system_params.Aw = 24.0; % m^2
system_params.kL = 0.1 * 180/pi; % 1/rad
system_params.kDO = 0.03; % unitless
system_params.kDI = 0.001 * (180/pi)^2; % 1/rad^2
system_params.rho = 1.225; % kg/m^3
system_params.m = 1134; % kg
system_params.g = 9.81; % m/s^2

% c) circular-loop simulation

T = 1881.157; % N
a = deg2rad(2); % deg -> rad
u = [T, a]; % thrust and AoA input

% Phase 1: simulate level-flight
options = odeset('AbsTol', 1e-9, 'MaxStep', 1e-3);

tspan = [0 2];
X0 = [0; 0; 61.331; 0];
[t1, X1] = ode45(@(t, X) simple_aircraft(t, X, u, system_params), tspan, X0, options);

% Phase 2: simulate circular loop flight
X0 = X1(end, :);
system_params.prev_X0 = X0;
tspan = [0 10];
[t2, X2] = ode45(@(t, X) simple_aircraft_loop(t, X, u, system_params), tspan, X0, options);

% estimate loop completion time
mask = (abs(X2(:,1) - X0(1)) < 1e-2) & (t2 > 5);
t_loop_done = min(t2(mask));

% trim t2 and X2 to only include the loop time interval
t2 = t2(t2<=t_loop_done);
X2 = X2(t2<=t_loop_done, :);

% Phase 3: simulate post-loop level-flight
X0 = X2(t2==t_loop_done, :);
tspan = [0 2];
[t3, X3] = ode45(@(t, X) simple_aircraft(t, X, u, system_params), tspan, X0, options);

t_total_1 = [t1; t2+t1(end); t3+t2(end)+t1(end)];
X_total_1 = [X1; X2; X3];

% Plotting
plot_aircraft(t_total_1, X_total_1);

% commanded AoA over time
V = sqrt(X2(:,3).^2 + X2(:,4).^2);
theta_loop = atan2(X2(:,4), X2(:,3));
AoA = (2/(system_params.rho * system_params.kL * system_params.Aw)) * ...
      (system_params.m/70 + (system_params.m*system_params.g*cos(theta_loop)) ./ V.^2);

AoA_complete = [2*ones(length(t1), 1); rad2deg(AoA); 2*ones(length(t3), 1)];

figure('Color','w');
plot(t_total_1, AoA_complete, 'LineWidth', 2);
grid on; grid minor;
set_ylim_padded(AoA_complete);
xlabel('Time (s)');
ylabel('Commanded AoA (deg)');
title('Commanded AoA vs. Time');

%% d) Constant AoA Loop

T = 1881.157; % N
a = deg2rad(2); % deg -> rad
u = [T, a]; % thrust and AoA input

% Phase 1: simulate level-flight
options = odeset('AbsTol', 1e-9, 'MaxStep', 1e-3);

tspan = [0 2];
X0 = [0; 0; 61.331; 0];
[t1, X1] = ode45(@(t, X) simple_aircraft(t, X, u, system_params), tspan, X0, options);

% Phase 2: simulate circular loop flight w/ constant AoA
a = deg2rad(10.75);
u = [T, a];
X0 = X1(end, :);
system_params.prev_X0 = X0;
tspan = [0 15];
[t_aoa, X_aoa] = ode45(@(t, X) simple_aircraft_loop_const_AoA(t, X, u, system_params), tspan, X0, options);

% estimate completion of the constant-AoA loop
mask = (abs(X_aoa(:,3) - X0(3)) < 1e-2) & (abs(X_aoa(:,4) - X0(4)) < 1e-2) & (t_aoa > 5);
t_loop_done = min(t_aoa(mask));

% trim t2 and X2 to only include the loop time interval
t_aoa = t_aoa(t_aoa<=t_loop_done);
X_aoa = X_aoa(t_aoa<=t_loop_done, :);

% Phase 3: simulate post-loop level-flight
a = deg2rad(2);
u = [T, a];
X0 = X_aoa(t_aoa==t_loop_done, :);
tspan = [0 2];
[t3, X3] = ode45(@(t, X) simple_aircraft(t, X, u, system_params), tspan, X0, options);

t_total_2 = [t1; t_aoa+t1(end); t3+t_aoa(end)+t1(end)];
X_total_2 = [X1; X_aoa; X3];

% Plotting
plot_aircraft(t_total_2, X_total_2);

%% Plot both

figure('Color', 'w');
plot(X_total_1(:, 1), X_total_1(:, 2), 'LineWidth', 2, 'DisplayName', 'Controlled AoA');
hold on;
plot(X_total_2(:, 1), X_total_2(:, 2), '--', 'LineWidth', 2, 'DisplayName', 'Constant AoA');
xlabel('Horizontal Displacement, x (km)');
ylabel('Altitude, y (m)');
title('Comparing Loop Trajectories')
legend();
grid on; grid minor;
axis equal;
hold off;

