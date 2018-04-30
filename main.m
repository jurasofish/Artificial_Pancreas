
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

% Time span to solve over. In Minutes.
t_start = 0;
t_end = 360;
tspan = [t_start, t_end];

% Initial values of the system.
Q_i1_0 = 0;
Q_i_0 = 0;
I_p_0 = K_i/Tau_i * Q_i_0;
G_0 = 6;
x_0 = 0;
G_s_0 = 6;
Q_m_0 = 13;
U_m_0 = 0;

sys_0 = [Q_i1_0 Q_i_0 I_p_0 G_0 x_0 G_s_0 Q_m_0 U_m_0];

% options for the ode solver
options = odeset('RelTol',1e-7,'Stats','on','OutputFcn',@odeplot);

% solve
[t,sys] = ode45(@(t,sys) sys_ode(t,sys,c), tspan, sys_0, options);

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
    G_0_const = 6; % name to differentiate from function initial value.
    
    % Si ~ G0 Ki P3 /P2 (mmol/L per units) is a positive insulin sen-
    % sitivity factor [the amount of glucose level drop (mmol/L) caused by
    % one unit of insulin]
    P3 = 2; % Unknown
    P2 = 3; % Unknown
    % S_i = G_0_const * K_i * P3 / P2;
    S_i = 1800/40; % guess
    
    % ksen (1/min) is the transfer-rate constant
    K_sen = 0.066;
    
    % xm (min) is a time constant characterizing the appearance of glucose
    % in the blood circulation from the gut.
    Tau_m = 19;
    
    
    const = [Tau_i K_i V PEGP p1 p2 p3 G_0_const S_i K_sen Tau_m];
end

function final_plot(t, sys)
% Plot the solution

figure % new figure Window
plot(t,sys(:,1),'-o',t,sys(:,2),'-.')
legend('First Compartment','Second Compartment (Plasma)','Location','northeast')
title('Insulin in Two Compartments')
xlabel('Time (minutes)')
ylabel('Insulin in Compartment (units)')

figure % new figure Window
plot(t,sys(:,3),'-o')
title('Insulin in Plasma')
xlabel('Time (minutes)')
ylabel('Insulin in Plasma (munits/L)')

figure % new figure Window
plot(t,sys(:,4),'-o')
title('Glucose in Plasma')
xlabel('Time (minutes)')
ylabel('Glucose in Plasma (mmol/L)')

figure % new figure Window
plot(t,sys(:,6),'-o')
title('Interstial Glucose')
xlabel('Time (minutes)')
ylabel('Interstial Glucose (mmol/L)')

figure % new figure Window
plot(t,sys(:,8),'-o')
title('Glucose Gut Absorption Rate')
xlabel('Time (minutes)')
ylabel('Glucose Gut Absorption Rate (?mol/kg/min)')



end