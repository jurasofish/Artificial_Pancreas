
c = get_constants();
    
% Get constants out of constants vector.
Tau_i = c(1);
K_i = c(2);
V = c(3);
PEGP = c(4);
P1 = c(5);
P2 = c(6);
P3 = c(7);

% Time span to solve over. In Minutes.
t_start = 0;
t_end = 180;
tspan = [t_start, t_end];

% Initial values of the system.
Q_i1_0 = 0;
Q_i_0 = 5;
I_p_0 = K_i/Tau_i * Q_i_0;

sys_0 = [Q_i1_0 Q_i_0 I_p_0];

% options for the ode solver
options = odeset('RelTol',1e-7,'Stats','on','OutputFcn',@odeplot);

% solve
[t,sys] = ode45(@(t,sys) sys_ode(t,sys,c), tspan, sys_0, options);

final_plot(t, sys)


function const = get_constants()
    
    % Time to peak of plasma insulin after bolus.
    Tau_i = 60;
    
    % Ki = 106 /^w KMCRh (10-3 min /L) is a gain inversely proportional to 
    % the metabolic clearance rate KMCR (mL/kg/ min) and the patient 
    % weight w (kg)
    K_i = 1;
    
    % V (mL/kg) is the glucose distribution volume
    V = 5000/80; % average person?
    
    % PEGP (mmol/L/min) describes the rate of endogenous production of glucose
    PEGP = 1/30;
    
    % p1 (1/min) describes glucose effectiveness (the ability of glucose to
    % promote its own disposal)
    P1 = 0.5
    
    
    % p2 (1/min) is a time constant characterizing the delay of the plasma
    % insulin effect on plasma glucose (deactivation rate of insulin effects) 
    P2 = 0.5;
    
    % p3 (1/min2 per munits/L) describes the activation rate of insulin
    % effects.
    P3 = 0.1;
    
    const = [Tau_i K_i V PEGP P1 P2 P3];
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

end