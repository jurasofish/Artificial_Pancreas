
c = get_constants();
    
% Get constants out of constants vector.
Tau_i = c(1);
K_i = c(2);
V = c(3);
PEGP = c(4);
p1 = c(5);
p2 = c(6);
p3 = c(7);
G_0_const = c(8);
S_i = c(9);
K_sen = c(10);
Tau_m = c(11);
K_bio = c(12);

% Time span to solve over. In Minutes.
t_start = 0;
t_end = 60*350;
tspan = [t_start, t_end];

% Initial values of the system.
Q_i1_0 = 0.01;
Q_i_0 = 0.01;
I_p_0 = K_i/Tau_i * Q_i_0;
G_0 = 7;
x_0 = 0;
G_s_0 = 7;
Q_m1_0 = 0;
Q_m_0 = 0;
U_m_0 = 0;

sys_0 = [Q_i1_0 Q_i_0 I_p_0 G_0 x_0 G_s_0 Q_m1_0 Q_m_0 U_m_0];

% options for the ode solver
options = odeset('RelTol',1e-7);

% The ODE solver will be run in a loop and solve over the period
% t+dt until it has solved over all of tspan.
%
% After every loop, that is, after every dt, the PID controls will be
% applied.
%
% The initial conditions at each time period are the final output of 
% the previous solution from the ODE. 
% e.g. The initial conditions when solving over [10, 10+dt] are the
% last row of sys that comes from solving over [10-dt, 10]

% dt is minutes of resolution between runs of the ODE solver.
% The PID controls are applied every dt minutes.
% O(n) time with dt, I think.
dt = 10e-0;

% Initialize sys such that sys(end, :) are the initial conditions
% for the next iteration of the loop. t needs a similar thing done,
% just to make sure that it still aligns with sys.
sys = sys_0; 
t = tspan(1);

setpoint = 7;
Kp = -0.00023;
Td = 0.000635;
Ti = 2500;
Ki = Kp/Ti; %0.002;
Kd = Kp/Td; % 0.03;

previous_error = 0;
integral = 0;
for tt = tspan(1):dt:tspan(2)
    
    tt % print to see progress.
    
    sys_old = sys; % So the new ODE outputs can be appended. 
    t_old = t; % So the new time values can be appended. 
    
    % Calculate PID stuff. Thanks wikipedia for the pseudocode.
    error = setpoint - sys(end,6);
    derivative = (error - previous_error) / dt;
    integral = integral + error * dt;
    
    % Insulin to administer -- output of PID.
    pid_insul = Kp * error + Ki * integral + Kd * derivative;
    
    % Can't administer negative insulin.
    if(pid_insul < 0)
        pid_insul = 0;
    end
    
    % For next loop iteration.
    previous_error = error;
    
    % Do the magic. 
    [t,sys] = ode45(@(t,sys) sys_ode(t,sys,c, pid_insul), [tt, tt+dt], sys(end, :), options);
    
    % Vertically concatenate fresh data (from end of period [t, t+dt] to the old
    % data. After all loops, t and sys will contain the data from all loops
    % and can simply be used as though the ODE was solved all at once.
    % ONLY the value from the end of the period [t, t+dt] is recorded.
    sys = [sys_old; sys(end, :)]; 
    t = [t_old; t(end)];
    
end

final_plot(t, sys)


function const = get_constants()
    
    % Time to peak of plasma insulin after bolus.
    Tau_i = 51;
    
    % Ki = 106 /^w KMCRh (10-3 min /L) is a gain inversely proportional to 
    % the metabolic clearance rate KMCR (mL/kg/ min) and the patient 
    % weight w (kg)
    K_i = 2;
    
    % V (mL/kg) is the glucose distribution volume
    V = 160; % average person?
    
    % PEGP (mmol/L/min) describes the rate of endogenous production of glucose
    PEGP = 0.0161;
    
    % p1 (1/min) describes glucose effectiveness (the ability of glucose to
    % promote its own disposal)
    p1 = 0.001;
    
    % p2 (1/min) is a time constant characterizing the delay of the plasma
    % insulin effect on plasma glucose (deactivation rate of insulin effects) 
    p2 = 5;
    
    % p3 (1/min2 per munits/L) describes the activation rate of insulin
    % effects.
    p3 = 0.1;
    
    % G0 (mmol/L) is an equi- librium point for glucose concentration.
    G_0_const = 100; % name to differentiate from function initial value.
    
    % Si ~ G0 Ki P3 /P2 (mmol/L per units) is a positive insulin sen-
    % sitivity factor [the amount of glucose level drop (mmol/L) caused by
    % one unit of insulin]
    S_i = 1800/30; % mg/dL/U
    S_i = S_i / 1000 * 10; % g/L/U
    S_i = S_i / 180; % mol/L/U
    S_i = S_i * 1000; % mmol/L/U
    
    % ksen (1/min) is the transfer-rate constant
    K_sen = 0.066;
    
    % xm (min) is a time constant characterizing the appearance of glucose
    % in the blood circulation from the gut.
    Tau_m = 19;
    
    % KBio is the carbohydrates bioavailability in the meal
    K_bio = 0.8;
    
    const = [Tau_i K_i V PEGP p1 p2 p3 G_0_const S_i K_sen Tau_m K_bio];
end

function final_plot(t, sys)
% Plot the solution

% figure('position', [0, 0, 600, 300]) % new figure Window
% plot(t,sys(:,1),'-o',t,sys(:,2),'-.')
% legend('First Compartment','Second Compartment (Plasma)','Location','northeast')
% title('Insulin in Two Compartments')
% xlabel('Time (minutes)')
% ylabel('Insulin in Compartment (units)')
% 
% figure('position', [0, 0, 600, 300]) % new figure Window
% plot(t,sys(:,3),'-o')
% title('Insulin in Plasma')
% xlabel('Time (minutes)')
% ylabel('Insulin in Plasma (munits/L)')
% 
% figure('position', [0, 0, 600, 300]) % new figure Window
% plot(t,sys(:,4),'-o')
% title('Glucose in Plasma')
% xlabel('Time (minutes)')
% ylabel('Glucose in Plasma (mmol/L)')
% 
% figure('position', [0, 0, 600, 300]) % new figure Window
% plot(t,sys(:,8),'-o')
% title('Glucose Gut Absorption Rate')
% xlabel('Time (minutes)')
% ylabel('Glucose Gut Absorption Rate (?mol/kg/min)')

figure('position', [0, 0, 600, 300]) % new figure Window
plot(t,sys(:,6),'-o')
title('Interstial Glucose')
xlabel('Time (minutes)')
ylabel('Interstial Glucose (mmol/L)')

end