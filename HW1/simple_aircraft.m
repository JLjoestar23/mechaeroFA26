function dXdt = simple_aircraft(t, X, u, system_params)
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
    
    theta_v = atan2(dy, dx); % velocity vector angle
    V = sqrt(dx^2 + dy^2); % speed
    Cl = kL*a; % lift coefficient
    Cd = kDO + kDI*a^2; % drag coefficient
    
    % calculate x-accel
    Tx = T*cos(a + theta_v);
    Lx = -0.5*rho*V^2*Cl*Aw*sin(theta_v);
    Dx = -0.5*rho*V^2*Cd*Aw*cos(theta_v);
    ddx = (1/m) * (Tx + Lx + Dx);
    
    % calculate y-accel
    Ty = T*sin(a + theta_v);
    Ly = 0.5*rho*V^2*Cl*Aw*cos(theta_v);
    Dy = -0.5*rho*V^2*Cd*Aw*sin(theta_v);
    Fg = -m*g;
    ddy = (1/m) * (Ty + Ly + Dy + Fg);
    
    % return state vector derivative
    dXdt = [dx; dy; ddx; ddy];
end