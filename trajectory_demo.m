%% trajectory_demo.m
% Color-mapped 3D trajectories + custom 3D vehicle models (patch meshes).
% Self-contained: generates synthetic loop data, no STL needed.
% Needs MATLAB R2020b+ for 'turbo' (swap for 'jet' on older versions).
clear; close all; clc; rng(0);

%% 1) Synthetic test data: a vertical loop with hover at start/end
T  = 6;                         % duration [s]
N  = 400;
t  = linspace(0, T, N)';
u  = t / T;
s  = 3*u.^2 - 2*u.^3;           % smoothstep -> zero speed at start/end
R  = 1.2;                       % loop radius [m]
th = 2*pi*s;

x = 3*(s - 0.5);                % drift along x
y = R*sin(th);
z = 2 - R*cos(th);              % z up here (flip axis if you use NED)
P = [x y z];                    % N x 3 positions

V     = [gradient(x, t), gradient(y, t), gradient(z, t)];  % velocity [m/s]
speed = vecnorm(V, 2, 2);       % color value for plot (a)

% "Experiment": reference + tracking error
Pexp = P + 0.05*randn(N,3);
err  = 0.5*(1 - cos(2*pi*3*s)) * 0.5 + 0.05*rand(N,1);  % color value for (b)

% Save in a plain CSV format you can swap for your own logs
% writematrix([t P V speed err], 'trajectory_data.csv');
% columns: t, x, y, z, vx, vy, vz, speed, err

%% 2) Simple procedural vehicle mesh (fuselage + wing), body frame x = forward
[V1, F1] = boxMesh([0 0 0],   [0.40 0.10 0.06]);   % fuselage
[V2, F2] = boxMesh([0 0 0.0], [0.10 0.70 0.015]);  % wing
[V3, F3] = boxMesh([-0.2 0 0.06], [0.05 0.01 0.12]); % tail fin
Vm = [V1; V2; V3];
Fm = [F1; F2 + size(V1,1); F3 + size(V1,1) + size(V2,1)];
Vm = 0.6 * Vm;                                     % overall scale

% To use a real CAD model instead:
%   TR = stlread('vehicle.stl'); Vm = TR.Points; Fm = TR.ConnectivityList;

%% 3) Plot
dtPose  = 0.7;                                     % interval between poses [s]
poseIdx = round(interp1(t, 1:N, 0:dtPose:T));
poseIdx = unique(poseIdx(~isnan(poseIdx)));

figure('Color','w');

% (a) reference colored by speed
ax1 = axes; hold on; grid on; box on;
colorLine(P, speed, 3);
drawVehicles(P, V, poseIdx, Vm, Fm);
setup3D(ax1, [0 2.5], 'Speed [m/s]');
title('(a) Reference');

% (b) experiment colored by tracking error
% ax2 = subplot(1,2,2); hold on; grid on; box on;
% colorLine(Pexp, err, 3);
% drawVehicles(Pexp, V, poseIdx, Vm, Fm);
% setup3D(ax2, [0 1], '||x_{ref} - x|| [m]');
% title('(b) Experiment');

%% ---------- local functions (must be at end of script) ----------
function colorLine(P, c, lw)
    % surface() trick: a 2-column "ribbon" with no face, interpolated edge color
    surface([P(:,1) P(:,1)], [P(:,2) P(:,2)], [P(:,3) P(:,3)], [c(:) c(:)], ...
            'FaceColor','none','EdgeColor','interp','LineWidth',lw);
end

function drawVehicles(P, V, idx, Vm, Fm)
    for k = idx(:)'
        f = V(k,:)' / max(norm(V(k,:)), 1e-9);     % forward = velocity direction
        a = [1;0;0];  a = a - (a'*f)*f;  a = a / norm(a);  % side axis
        b = cross(f, a);                           % up axis
        Rb = [f a b];                              % body -> world rotation
        Vt = (Rb * Vm')' + P(k,:);
        patch('Faces',Fm,'Vertices',Vt, ...
              'FaceColor',[0.30 0.32 0.42],'EdgeColor','none', ...
              'FaceLighting','gouraud','AmbientStrength',0.5);
    end
end

function setup3D(ax, lim, label)
    axis(ax,'equal'); view(ax, 35, 22);
    xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
    colormap(ax, turbo); clim(ax, lim);
    cb = colorbar(ax,'southoutside'); cb.Label.String = label;
    camlight(ax,'headlight'); material(ax,'dull');
end

function [V, F] = boxMesh(c, s)
    % axis-aligned box: center c, size s -> 8 vertices, 6 quad faces
    h = s/2;
    V = c + h .* [-1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1; -1 -1 1; 1 -1 1; 1 1 1; -1 1 1];
    F = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
end
