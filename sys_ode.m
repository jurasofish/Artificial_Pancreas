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
    
    % Get constants defined by functions.
    U_i = get_U_i(t);

    % Get values of functions from sys vector.
    Q_i1 = sys(1);
    Q_i  = sys(2);

    % Evaluate derivatives.
    d_Q_i1 = -1*Q_i1/Tau_i + U_i;
    d_Q_i = Q_i1/Tau_i - Q_i/Tau_i;
    d_K_i = K_i/Tau_i * d_Q_i;
    
    sys_deriv = zeros(size(sys)); % Initialize memory.
    % populate sys_deriv
    sys_deriv(1) = d_Q_i1;
    sys_deriv(2) = d_Q_i;
    sys_deriv(3) = d_K_i;
end

function U_i = get_U_i(t)
% Ui(t) (unit/min) is the external insulin infusion rate
    if(t >= 60 && t < 61)
        U_i = 5;
    else
        U_i = 0;
    end
end












