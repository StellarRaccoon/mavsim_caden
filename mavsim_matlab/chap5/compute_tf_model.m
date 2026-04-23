%  Create a function that computes the transfer function models
% described in this chapter, linearized about the trim state and
% trim inputs.

    % - Implment "compute_tf_model()" and the dependent functions "dT_dVa()",
    % "dT_ddelta_t()" in "compute_models.py".

function [T_phi_delta_a,...
            T_chi_phi,...
            T_theta_delta_e,...
            T_h_theta,...
            T_h_Va,...
            T_Va_delta_t,...
            T_Va_theta,...
            T_v_delta_r]...
            = compute_tf_model(x_trim,u_trim,MAV)
%I assume the original arg P is actually MAV
pn  = x_trim(1);
pe  = x_trim(2);
pd  = x_trim(3);
u   = x_trim(4);
v   = x_trim(5);
w   = x_trim(6);
phi     = x_trim(7);
theta   = x_trim(8);
psi     = x_trim(9);
p   = x_trim(10);
q   = x_trim(11);
r   = x_trim(12);
delta_e_trim = u_trim(1);
delta_a_trim = u_trim(2);
delta_r_trim = u_trim(3);
delta_t_trim = u_trim(4);

%get Va
Va_trim = sqrt(u^2+v^2+w^2);
    % d_phi_1 = q*sin(phi)*tan(theta)+r*cos(phi)*tan(theta);
    a_phi1 = -1/2*MAV.rho*Va_trim^2*MAV.S_wing*MAV.C_Y_p*b/(2*Va_trim);
    a_phi2 = 1/2*MAV.rho*Va_trim^2*MAV.S_wing*MAV.C_Y_delta_a;

    % d_phi2 = Gamma1*p*q-Gamma2*q*r...
    %         +/2*rho*Va^2*MAV.S_wing*b(C_Y_0+C_Y_beta*beta- MAV.C_Y_p*b/(2*Va)*d_phi_1+C_Y_r*b*r/(2*Va)+C_Y_delta_r*del_r)...
    %         - Gamma3*Q_p+d_phi_1;

    % tf num and  denominator - vector of the coiffecents of s in decending order
    % a_pi2/(s^2+a_pi1*s+0)
    T_phi_delta_a = tf([a_phi2], [1 a_phi1 0]);

    T_chi_phi = tf([MAV.gravity/Va_trim], [1 0]);
    
    a_theta1 = (-MAV.rho*Va_trim^2*c*MAV.S_wing/(2*Jy))*MAV.C_m_q*c/(2*Va_trim);
    a_theta2 = (-MAV.rho*Va_trim^2*c*MAV.S_wing/(2*Jy))*MAV.C_m_alpha;
    a_theta3 = (-MAV.rho*Va_trim^2*c*MAV.S_wing/(2*Jy))*MAV.C_m_delta_e;
    %a_theta3 *s /s^2+atheta1*s+a_theta2
    T_theta_delta_e = tf([a_theta3], [1 a_theta1 a_theta2]);

    %Va/s
    T_h_theta = tf([Va_trim], [1 0]);
    %theta/s
    T_h_Va = tf([theta_trim], [1 0]);

% velocities
dT_dVa = dT_dVa(mav, Va_trim, delta_t_trim, MAV);
dT_ddelta_t = dT_ddelta_t(mav, Va_trim, delta_t_trim, MAV);

a_V1 = (-MAV.rho*Va_trim^2*c*MAV.S_wing/MAV.mass)*(MAV.C_D_0+MAV.C_D_alpha*alpha_trim+MAV.C_D_delta_e*delta_e_trim)...
        -(1/MAV.mass)*dT_dVa;
a_V2 = (1/MAV.mass)*dT_ddelta_t;
a_V3 = MAV.gravity*cos(theta_trim-alpha_trim);

T_Va_delta_t    = tf([a_V2],[1,a_V1]);
T_Va_theta      = tf([-a_V3],[1,a_V1]);
% TO DO FIND TV DELTA R
% T_v_delta_r     = tf([a_beta2],[1,a_beta1]);
save transfer_function_coef.mat a_phi1 a_phi1 a_theta1 a_theta2 a_theta3 a_V1 a_V2 a_V3;
%apptox values
% a_phi_2 = 130.6
% a_phi_1 = 22.6
% a_theta_1 = 5.288
% a_theta_2 = 99.7
% a_theta_3 = -36.02
% a_v_1 = 0.6607
% a_v_2 = 47.02
end

% returns the derivative of motor thrust with respect to Va
% use finite central difference
function dThurst = dT_dVa(mav, Va, delta_t, MAV)
    % set step very small
    dVa = 10^-8;
    % find Tp using Va+dVa
    [T_p1, Q_prop] = prop_force_and_torque(Va+dVa, delta_t, MAV);
    % find Tp using Va-dVa
    [T_m1, Q_prop] = prop_force_and_torque(Va-dVa, delta_t, MAV);
    % apply central difference formula
    dThrust = (T_p1-T_m1)/(2*dVa);
end

% returns the derivative of motor thrust with respect to delta_t
function dThurst = dT_ddelta_t(mav, Va, delta_t)
    % set step very small
    ddetlta_t = 10^-8;
    % find Tp using Va+dVa
    [T_p1, Q_prop] = prop_force_and_torque(Va, delta_t+ddetlta_t, MAV);
    % find Tp using Va-dVa
    [T_m1, Q_prop] = prop_force_and_torque(Va, delta_t-ddetlta_t, MAV);
    % apply central difference formula
    dThrust = (T_p1-T_m1)/(2*ddetlta_t);
end