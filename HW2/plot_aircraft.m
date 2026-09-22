function plot_aircraft(t, X)
    x  = X(:,1);
    y  = X(:,2);
    vx = X(:,3);
    vy = X(:,4);
    V  = sqrt(vx.^2 + vy.^2);
    
    set(0, 'DefaultAxesFontSize', 12, 'DefaultAxesFontName', 'Helvetica', ...
           'DefaultLineLineWidth', 1.8, 'DefaultAxesBox', 'on', ...
           'DefaultAxesLineWidth', 1.0, 'DefaultAxesTickLength', [0.008 0.008]);
     
    navy   = [0.00 0.20 0.45];
    orange = [0.85 0.40 0.05];
    green  = [0.10 0.55 0.30];
    grey   = [0.45 0.45 0.45];
     
    % Figure 1: Flight path (y vs x)
    figure('Color','w','Position',[100 100 850 450]);
    plot(x/1000, y, 'Color', navy);
    grid on; grid minor;
    xlabel('Horizontal Distance, x (km)');
    ylabel('Altitude Change, y (m)');
    title('X-Y Flight Path');
    % ylim_pad = max(5, 0.1*max(abs(y))+1);
    % ylim([-ylim_pad, ylim_pad]);
    set(gca,'GridAlpha',0.15,'MinorGridAlpha',0.08);
     
    % Figure 2: Velocity components vs time
    figure('Color','w','Position',[100 100 900 600]);
     
    subplot(3,1,1);
    plot(t, vx, 'Color', navy); grid on; grid minor;
    ylabel('v_x (m/s)');
    title('Velocity Components vs. Time');
    set(gca,'GridAlpha',0.15,'MinorGridAlpha',0.08);
     
    subplot(3,1,2);
    plot(t, vy, 'Color', orange); grid on; grid minor;
    ylabel('v_y (m/s)');
    set(gca,'GridAlpha',0.15,'MinorGridAlpha',0.08);
     
    subplot(3,1,3);
    plot(t, V, 'Color', green); grid on; grid minor;
    ylabel('|V| (m/s)');
    xlabel('Time (s)');
    set(gca,'GridAlpha',0.15,'MinorGridAlpha',0.08);
     
    % Figure 3: Position components vs time
    figure('Color','w','Position',[100 100 900 450]);
     
    subplot(2,1,1);
    plot(t, x/1000, 'Color', navy); grid on; grid minor;
    ylabel('x (km)');
    title('Position Components vs. Time');
    set(gca,'GridAlpha',0.15,'MinorGridAlpha',0.08);
     
    subplot(2,1,2);
    plot(t, y, 'Color', orange); grid on; grid minor;
    ylabel('y (m)');
    xlabel('Time (s)');
    set(gca,'GridAlpha',0.15,'MinorGridAlpha',0.08);
end