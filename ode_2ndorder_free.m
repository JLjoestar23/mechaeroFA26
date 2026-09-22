% ode_2ndorder_free.m
%
% This code solves ode's by numerical integration using ode45
% e.g., single degree-of-freedom mass-spring-damper oscillator
%
% last modified 1/9/22 CLee
%
function ode_2ndorder_free
% clear out memory 
clear all
close all
clear functions

% Define parameters of single dof system in m-c-k form
% Specify values for  mass,  damping constant, and spring constant,
% m =  1.0;            
% c =  0.0;
% k =  4.0;
% wn = sqrt(k/m);
% zeta = c/(2*wn*m);
% wn2 = k/m;

% Define parameters of single dof system in terms of damping ratio
% and natural frequency
wn =   4;                     % ENTER value for natural frequency
zeta = 2;                     % ENTER value for damping ratio

wn2 = wn*wn;                % natural frequency squared (= k/m)
wd = wn*sqrt(1-zeta^2);     % damped natural frequency 
%
% roots of characteristic eqns (lambda 1 and 2)
r1 = -zeta*wn + sqrt(zeta^2-1)*wn;
r2 = -zeta*wn - sqrt(zeta^2-1)*wn;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify the time span for the simulation
tend = 30;                               % CHANGE end time to scale response 
t_span = [0, tend];                     

% Specify the initial conditions: intial values of the state variables
x0 =    1;           % initial displacement                %ENTER IC value 
v0 =    -100;          % initial velocity                     %ENTER IC value
z_0 = [x0, v0];      % store as a vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define two state variables, Z_1 = x and Z_2 = x_dot,

% call od45 to solve eqns of motion which are defined in the function 
% sdof_fun. Include the time span and intitial conditions
% state values (displacement and velocity) are returned in the vector zout
% (n rows by 2 columns). t is the corresponding time values (n rows by 1
% column
[t, zout] = ode45(@sdof_fun, t_span, z_0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analytical solutions for free response          % PUT ANALYTIC SOLNS HERE
%
% check zeta values to define class of response
    if  zeta == 0                                                 %undamped
xanalytic_free = v0/wn*sin(wn*t) + (x0)*cos(wn*t);
    elseif zeta > 1                                             %overdamped
xanalytic_free = (x0*r2-v0)/(r2-r1)*exp(r1*t)...
                 + (x0*r1-v0)/(r1-r2)*exp(r2*t);    
    elseif zeta == 1                                     %critically damped
xanalytic_free = exp(-wn*t).*( (v0+x0*wn)*t + x0 );   
    elseif  zeta < 1 && zeta > 0                               %underdamped
xanalytic_free = exp(-zeta*wn*t).*( x0*cos(wd*t) ...
                  +  (v0+x0*zeta*wn)/wd*sin(wd*t));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot time responses
% 
figure (1)           % define a new plot
subplot(2,1,1)       % 2 plots on one page: top=displacement, bottom=velocity
plot( t, zout(:,1) ) % zout(:,1) is the column with positions

hold                             % hold plot to overlay a 2nd set of points
plot(t,xanalytic_free,'r+')      % overlay analytic solution 

xlabel('Time')
ylabel('Displacement')
title('SDOF Free Response with ICs')
legend('numerical','analytic')
%
subplot(2,1,2)            % second plot of 2 row and 1 column 
plot( t, zout(:,2) )      % zout(:,2) is the column with velocities
xlabel('Time ')
ylabel('Velocity')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define eqn or motion as a series of first order eqns

function dzdt = sdof_fun(T, z)         % define function 
    
% second order oscillator in first order, state space form
% T is the current time value (within t_span)
% z is the current value of the state variable:
% z(1)= displacement, z(2) = velocity
%
dzdt = zeros(2,1);         % 2x1 vector of first-order derivatives of states
% Eqns of motion in first order form 
dzdt(1) =   z(2);                                          
dzdt(2) =  -wn2*z(1) - 2*zeta*wn*z(2);                     %PUT  EOM's HERE
% 

%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end

