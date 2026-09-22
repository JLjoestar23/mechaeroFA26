% ode_time_2ndorder.m
%
% This code solves a 2nd order ODE by numerical integration using ode45
% e.g., single degree-of-freedom mass-spring-damper oscillator
%
% numerical and analytical(displacement only) solutions for
%   free response
%   step response = F/m
%   harmonic response = F/m*cos(Omega*t)
%
%   last modified 1/9/22 CL
%  
function ode_time_2ndorder
% clear out memory 
clear all
close all
clear functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify system parameters in terms of mass, damping, and spring constants
% m =  1.0;            
% c =  0.0;
% k =  4.0;
% wn = sqrt(k/m);
% zeta = c/(2*wn*m);
% wn2 = k/m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECIFY system parameters in terms of damping ratio and natural frequency
%
wn   = 1.2;                     % ENTER value for natural frequency
zeta = 1;                      % ENTER value for damping ratio
%
wn2  = wn*wn;                   % natural freq squared
wd = wn*sqrt(1-zeta^2);         % damped natural frequency 
coverm = 2*zeta*wn
koverm = wn2

%
% roots of characteristic eqns (lambda 1 and 2)
r1 = -zeta*wn + sqrt(zeta^2-1)*wn;
r2 = -zeta*wn - sqrt(zeta^2-1)*wn;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECIFY excitation case: free=1, step=2, harmonic=3
casenum = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define parameters for excitation terms 
%  ENTER values 
Fmstep  = 2;       % magnitude of the step function = F/m
Fmharm  = 10;       % harmonic, external excitation amplitude = F/m
Omega   = 1;       % harmonic, external excitation frequency in rad/s
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECIFY the time span for the simulation
tend = 40;                             % CHANGE end time to scale response 
t_span = [0, tend];                     

% SPECIFY the initial conditions
x0 =    1.0;           % initial displacement        %ENTER IC value 
v0 =    -3;             % initial velocity             %ENTER IC value
z_0 = [x0, v0];        % store as a vector
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define two state variables, Z_1 = x and Z_2 = x_dot,
%
% call od45 to solve eqns of motion which are defined below in the function 
% sdof_fun. Include the time span and intitial conditions
% state values (displacement and velocity) are returned in the vector zout
% (n rows by 2 columns). t is the corresponding time values (n rows by 1
% column

[t, zout] = ode45(@sdof_fun, t_span, z_0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analytical FREE response solutions                  
%
if casenum == 1
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
xanalytic = xanalytic_free;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analytical STEP response solutions
if casenum == 2
    x0step   = x0 - Fmstep/wn2;  % transient response correction for step
    ss_step  = Fmstep/wn2        % steady-state step value
%   
    if  zeta == 0                                                 %undamped
xanalytic_step = v0/wn*sin(wn*t) + (x0step)*cos(wn*t) + ss_step;
    elseif zeta > 1                                             %overdamped
xanalytic_step = (x0step*r2-v0)/(r2-r1)*exp(r1*t) ...
          + (x0step*r1-v0)/(r1-r2)*exp(r2*t) + ss_step;    
    elseif zeta == 1                                     %critically damped
xanalytic_step = exp(-wn*t).*( (v0+x0step*wn)*t + x0step )+ ss_step;   
    elseif  zeta < 1 && zeta > 0                                %underdamped
xanalytic_step = exp(-zeta*wn*t).*( x0step*cos(wd*t) ...
                 +  (v0+x0step*zeta*wn)/wd*sin(wd*t)) + ss_step;
 end        
  xanalytic = xanalytic_step;  
end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analytical HARMONIC response solutions
if casenum == 3
    droot = sqrt( (wn^2-Omega^2)^2 + (2*zeta*wn*Omega)^2 );
    ph3 = atan2( (2*zeta*wn*Omega),(wn^2-Omega^2) );       %% use atan2!!!!
    ovwd = sqrt(zeta^2-1)*wn;
    zwn = zeta*wn;
    x0A = x0-Fmharm*cos(ph3)/droot;
    v0A = v0 + zwn*x0- Fmharm*(Omega*sin(ph3)+zwn*cos(ph3))/droot;
    xharmss = Fmharm*cos(Omega*t-ph3)/droot;
%   
    if  zeta == 0                                                 %undamped
        if Omega == wn
            xanalytic_harm = x0*cos(wn*t) + v0/wn*sin(wn*t) ...
                             + Fmharm/(2*wn).*t.*sin(wn*t);
        else
            xanalytic_harm = exp(-zeta*wn*t).*( x0A*cos(wd*t)...
                            + 1./wd*v0A.*sin(wd*t) ) + xharmss;
        end
    elseif zeta > 1                                             %overdamped        
xanalytic_harm =  exp(-zwn*t).*(  x0A.*cosh(ovwd*t) ...
                  + 1./ovwd*v0A.*sinh(ovwd*t)  )   + xharmss;     
    elseif zeta == 1                                     %critically damped
xanalytic_harm = exp(-zwn*t).*(x0A + v0A.*t)  + xharmss;    
    elseif  zeta < 1 && zeta > 0                                %underdamped       
xanalytic_harm = exp(-zeta*wn*t).*( x0A*cos(wd*t)...
                  + 1./wd*v0A.*sin(wd*t) ) + xharmss;
 end        
  xanalytic = xanalytic_harm;  
end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if casenum ~= 1 && casenum ~= 2 && casenum ~= 3
   xanalytic = zeros(length(t),1);                  % if case is not 1,2,3
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot time responses
figure (1)  
subplot(2,1,1)
plot( t, zout(:,1) )         % zout(:,1) is the column with positions

% plot analytical solns
hold
plot(t,xanalytic,'r+')       % plot analytical soln values as red + signs

xlabel('Time')
ylabel('Displacement')
title('Total Response')
%
%
subplot(2,1,2)
plot( t, zout(:,2) )          % zout(:,1) is the column with velocities
xlabel('Time ')
ylabel('Velocity')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define eqn of motion as a series of first order eqns

function dzdt = sdof_fun(T, z)
% z = state variables: z(1) = displacement, z(2) = velocity
    
% assign appropriate forcing term 
if casenum == 1                       % free response
    f_in = 0;    
elseif casenum == 2
    f_in = Fmstep;                    % step response      
elseif casenum == 3
    f_in = Fmharm*cos(Omega*T);       % harmonic response (cos) 
else
    f_in = 0;                         % in case of incorrect input
end
%
% define dzdt as a vector with 2 rows and 1 column 
dzdt = zeros(2,1);            % "slope" or first-order derivative of states
%
% Eqns of motion in first order form 
% 
dzdt(1) = z(2);                               
dzdt(2) = -wn2*z(1) - 2*zeta*wn*z(2) + f_in;    
% 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

