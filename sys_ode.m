function sys_deriv = sys_ode(t, sys, c)
    % t is the current time.
    % sys is a vector containing the values of the of the
    % functions in the system.
    % sys_deriv is a vector containing the values of the derivatives of the
    % functions in the system.
    % c is array of constants.
    % This function is passed to the ode solver.
    
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
    
    % Get constants defined by functions.
    U_i = get_U_i(t);
    D_m = get_D_m(t);

    % Get values of functions from sys vector.
    Q_i1 = sys(1);
    Q_i  = sys(2);
    I_p = sys(3);
    G = sys(4);
    x = sys(5);
    G_s = sys(6);
    Q_m1 = sys(7);
    Q_m = sys(8);
    U_m = sys(9);

    % Evaluate derivatives.
    d_Q_i1 = -1*Q_i1/Tau_i + U_i;
    d_Q_i = Q_i1/Tau_i - Q_i/Tau_i;
    d_I_p = K_i/Tau_i * d_Q_i;
    d_G = -p1*G - S_i/Tau_i*Q_i + PEGP + U_m/V;
    d_x = -p2*x + p3 * I_p;
    d_G_s = K_sen*(G-G_s);
    d_Q_m1 = -1/Tau_m*Q_m1 + K_bio*D_m;
    d_Q_m = 1/Tau_m*Q_m1 - 1/Tau_m*Q_m;
    d_U_m = 1/Tau_m * d_Q_m;
    
    sys_deriv = zeros(size(sys)); % Initialize memory.
    % populate sys_deriv
    sys_deriv(1) = d_Q_i1;
    sys_deriv(2) = d_Q_i;
    sys_deriv(3) = d_I_p;
    sys_deriv(4) = d_G;
    sys_deriv(5) = d_x;
    sys_deriv(6) = d_G_s;
    sys_deriv(7) = d_Q_m1;
    sys_deriv(8) = d_Q_m;
    sys_deriv(9) = d_U_m;
    
end

function U_i = get_U_i(t)
% Ui(t) (unit/min) is the external insulin infusion rate
    if(t >= 1 && t < 10)
        U_i = 0.75;
    else
        U_i = 0;
    end
end

function D_m = get_D_m(t)
    % Dm (t) (?mol/kg/min) is the rate of glucose ingestion
    
    % QCHO (g) is the quantity of carbohydrates ingested
    Q_cho = 60;
    
    % MCHO = 180.156 (g/mol) is the molar mass of glucose
    M_cho = 180.156;
    
    % w (kg) is the patient weight
    w = 80;
    if(t >= 1 && t < 10)
        D_m = 1e6 * Q_cho/(w*M_cho)/10;
        
    else
        D_m = 0;
    end
end









