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
    P1 = c(5);
    P2 = c(6);
    P3 = c(7);
    G_0_const = c(8);
    S_i = c(9);
    K_sen = c(10);
    Tau_m = c(11);
    
    % Get constants defined by functions.
    U_i = get_U_i(t);
    G_emp = get_G_emp(t);

    % Get values of functions from sys vector.
    Q_i1 = sys(1);
    Q_i  = sys(2);
    I_p = sys(3);
    G = sys(4);
    x = sys(5);
    G_s = sys(6);
    Q_m = sys(7);
    U_m = sys(8)

    % Evaluate derivatives.
    d_Q_i1 = -1*Q_i1/Tau_i + U_i;
    d_Q_i = Q_i1/Tau_i - Q_i/Tau_i;
    d_I_p = K_i/Tau_i * d_Q_i;
    d_G = -P1*G - S_i/Tau_i*Q_i + PEGP + U_m/V;
    d_x = -P2*x + P3 * I_p;
    d_G_s = K_sen*(G-G_s);
    d_Q_m = -1/Tau_m * Q_m + G_emp;
    d_U_m = 1/Tau_m * d_Q_m;
    
    sys_deriv = zeros(size(sys)); % Initialize memory.
    % populate sys_deriv
    sys_deriv(1) = d_Q_i1;
    sys_deriv(2) = d_Q_i;
    sys_deriv(3) = d_I_p;
    sys_deriv(4) = d_G;
    sys_deriv(5) = d_x;
    sys_deriv(6) = d_G_s;
    sys_deriv(7) = d_Q_m;
    sys_deriv(8) = d_U_m;
    
end

function U_i = get_U_i(t)
% Ui(t) (unit/min) is the external insulin infusion rate
    if(t >= 60 && t < 61)
        U_i = 5;
    else
        U_i = 0;
    end
end


function G_emp = get_G_emp(t)
    % gastric emptying.
    % does this need to be adjusted for a meal to start at any time?
    %  Dm (?mol/kg) represents the total amount of ingested glucose
    D_m = 10;

    % KBio is the carbohydrates bioavailability in the meal
    K_bio = 1.5;

    % Tasc , Tdes , and T (min) are the durations of ascend- ing and
    % descending rates and the total duration of gastric emptying,
    % respectively.
    T_asc = 60;
    T_des = 45;
    T = T_asc + T_des;
    const = K_bio*D_m/(T - (T_des+T_asc)/2);

    if(t < T_asc)
        G_emp = const * t/T_asc;
    elseif(T_asc <= t <= T-T_des)
        G_emp = const * 1;
    elseif(T-T_des<t)
        G_emp = const * (T-t)/T_des;
    elseif(T<t)
        G_emp = const * 0;
    end
end









