function dXdt = simple_aircraft_plf(t, X, u, system_params)
    Aw = system_params.Aw ; % m^2
    kL = system_params.kL; % 1/rad
    kDO = system_params.kDO; % unitless
    kDI = system_params.kDI; % 1/rad^2
    rho = system_params.rho; % kg/m^3
    m = system_params.m; % kg
    g = system_params.g; % m/s^2
    
    % system inputs
    T = u(1); % thrust input
    a = u(2); % AoA input
    
    % system states
    x = X(1); % x-pos
    y = X(2); % y-pos
    dx = X(3); % x-vel
    dy = X(4); % y-vel

    theta_v = atan2(dy, dx);
    V = sqrt(dx^2 + dy^2);

    pitch_adj_alt = 4100;
    y_des = 4200;
    if y > pitch_adj_alt && y < y_des
        a = deg2rad(2 - (y-pitch_adj_alt)/(y_des-pitch_adj_alt));
        % a = deg2rad(2);
    elseif y >= y_des
        a = deg2rad(2);
    end
    Cl = kL*a;
    Cd = kDO + kDI*a^2;
    
    % aerodynamic forces computed once
    Lx = -0.5*rho*V^2*Cl*Aw*sin(theta_v);
    Dx = -0.5*rho*V^2*Cd*Aw*cos(theta_v);
    Ly =  0.5*rho*V^2*Cl*Aw*cos(theta_v);
    Dy = -0.5*rho*V^2*Cd*Aw*sin(theta_v);
    Fg = -m*g;
    
    if V <= 0
        Lx = 0;
        Dx = 0;
        Ly = 0;
        Dy = 0;
    end
    
    if y > pitch_adj_alt && y < y_des
        % partial linear feedback controller for y-accel
        % desired y-accel as an input is affine to system dynamics
        Kp = 0.1;
        Kd = 1;
        ddy_des = -Kp*(y - y_des) - Kd*dy;
        T = (m*ddy_des - Ly - Dy - Fg) / sin(a + theta_v);
    elseif y >= y_des
        T = 1881.157;
    end
    
    % now compute thrust components once, with final T and a
    Tx = T*cos(a + theta_v);
    Ty = T*sin(a + theta_v);
    
    ddx = (1/m) * (Tx + Lx + Dx);
    ddy = (1/m) * (Ty + Ly + Dy + Fg);
    
    if y <= 0 && ddy < 0
        dy = 0;
        ddy = 0;
    end
    
    dXdt = [dx; dy; ddx; ddy];
end