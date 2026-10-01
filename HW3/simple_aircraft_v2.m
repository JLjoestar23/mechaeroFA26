function dXdt = simple_aircraft_v2(t, X, u, system_params)
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
    
    % insant pitch down at y=4200m
    pitch_adj_alt = 4200;
    if y >= pitch_adj_alt
        a = deg2rad(2);
    end

    Cl = kL*a;
    Cd = kDO + kDI*a^2;
    
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
    
    % switch thrust
    if y >= pitch_adj_alt
        T = 1881.157;
    end

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