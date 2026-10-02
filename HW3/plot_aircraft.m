function plot_aircraft(t, X)
    x  = X(:,1);
    y  = X(:,2);
    vx = X(:,3);
    vy = X(:,4);
    V  = sqrt(vx.^2 + vy.^2);
    
    % set(0, 'DefaultAxesFontSize', 12, 'DefaultAxesFontName', 'Helvetica', ...
    %        'DefaultLineLineWidth', 1.8, 'DefaultAxesBox', 'on', ...
    %        'DefaultAxesLineWidth', 1.0, 'DefaultAxesTickLength', [0.008 0.008]);
     
    navy   = [0.00 0.20 0.45];
    orange = [0.85 0.40 0.05];
    green  = [0.10 0.55 0.30];
    % grey   = [0.45 0.45 0.45];
     
    % Figure 1: Flight path (y vs x)
    figure('Color','w');
    plot(x/1000, y, 'Color', navy, 'LineWidth', 2);
    grid on; grid minor;
    set_ylim_padded(y);
    xlabel('Horizontal Distance, x (km)');
    ylabel('Altitude Change, y (m)');
    title('X-Y Flight Path');
     
    % Figure 2: Velocity components vs time
    figure('Color','w');
     
    subplot(3,1,1);
    plot(t, vx, 'Color', navy, 'LineWidth', 2);
    grid on; grid minor;
    set_ylim_padded(vx);
    ylabel('v_x (m/s)');
    title('Velocity Components vs. Time');
    
    subplot(3,1,2);
    plot(t, vy, 'Color', orange, 'LineWidth', 2);
    grid on; grid minor;
    set_ylim_padded(vy);
    ylabel('v_y (m/s)');
     
    subplot(3,1,3);
    plot(t, V, 'Color', green, 'LineWidth', 2);
    grid on; grid minor;
    set_ylim_padded(V)
    ylabel('|V| (m/s)');
    xlabel('Time (s)');

    % Figure 3: Position components vs time
    figure('Color','w');
     
    subplot(2,1,1);
    plot(t, x/1000, 'Color', navy, 'LineWidth', 2);
    grid on; grid minor;
    set_ylim_padded(x/1000);
    ylabel('x (km)');
    title('Position Components vs. Time');
     
    subplot(2,1,2);
    plot(t, y, 'Color', orange, 'LineWidth', 2);
    grid on; grid minor;
    set_ylim_padded(y)
    ylabel('y (m)');
    xlabel('Time (s)');
end