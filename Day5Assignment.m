M = [9 0; 0 1];
K = [27 -3; -3 3];

[V, D] = eig(K, M);

wn1 = sqrt(D(1,1));
wn2 = sqrt(D(2,2));

U1 = V(:,1)/V(1,1);
U2 = V(:,2)/V(1,2);

x0 = [2; 3];
dx0 = [0; 0];


fun = @(z) [z(1)*U1*cos(z(3)) + z(2)*U2*cos(z(4)) - x0; ...
            -z(1)*U1*wn1*sin(z(3)) - z(2)*U2*wn2*sin(z(4)) - dx0];

soln = fsolve(fun, [x0; dx0]);

c1 = soln(1);
c2 = soln(2);
phi1 = soln(3);
phi2 = soln(4);

t_ana = linspace(0, 20, 500);
X_ana = c1*U1*cos(wn1*t_ana + phi1) + c2*U2*cos(wn2*t_ana + phi2);

tspan = [0 20];
[t, X] = ode45(@double_mass_spring, tspan, [x0; dx0]);

figure();
plot(t_ana, X_ana(1,:))
hold on;
plot(t_ana, X_ana(2,:));
plot(t, X(:,1), '--', t, X(:,2), '--');
title('State Trajectory of Double Mass-Spring System')
xlabel('Time (s)');
ylabel('Displacement (m)');
legend('Mass 1', 'Mass 2');
hold off;


function dXdt = double_mass_spring(t, X)
    q = X(1:2);
    dq = X(3:4);

    M = [9 0; 0 1];
    K = [27 -3; -3 3];

    ddq = (-M\K)*q;

    dXdt = [dq; ddq];
end