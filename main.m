
c = get_constants();
    
% Get constants out of constants vector.
Tau_i = c(1);
K_i = c(2);

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
    
    const = [Tau_i K_i];
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