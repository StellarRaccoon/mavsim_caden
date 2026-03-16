%  Create a function that computes the transfer function models
% described in this chapter, linearized about the trim state and
% trim inputs.


function [T_phi_delta_a,...
          T_chi_phi,...
          T_theta_delta_e,...
          T_h_theta,...
          T_h_Va,...
          T_Va_delta_t,...
          T_Va_theta,...
          T_v_delta_r]...
          = compute_tf_model(x_trim,u_trim,P)
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
del_e = u_trim(1);
del_a = u_trim(2);
del_r = u_trim(3);
del_t = u_trim(4);
%u= del e, del a, del r, del t
%TO DO WHAT IS VA HOW TO GET VA
    % i think C_p is C_Y??
    % TO DO WHAT IS Q_p
    d_phi_1 = q*sin(phi)*tan(theta)+r*cos(phi)*tan(theta);
    a_phi1 = -1/2*rho*Va^2*S_wing*C_Y_p*b/(2*Va);
    a_phi2 = 1/2*rho*Va^2*S_wing*C_Y_delta_a;
    d_phi2 = Gamma1*p*q-Gamma2*q*r+1/2*rho*Va^2*S_wing*b(C_Y_0+C_Y_beta*beta- C_Y_p*b/(2*Va)*d_phi_1+C_Y_r*b*r/(2*Va)+C_Y_delta_r*del_r)- Gamma3*Q_p+d_phi_1;
    % tf num and  denominator - vector of the coiffecents of s in decending order
    % a_pi2/(s^2+a_pi1*s+0)
    T_phi_delta_a = tf([a_phi2], [1 a_phi1 0]);

    %g/Vg / s
    %TO DO HOW TO GET VG
    T_chi_phi = tf([gravity/Vg], [1 0]);
    
    a_theta1 = (-rho*Va^2*c*S_wing/(2*Jy))*C_m_q*c/(2*Va);
    a_theta2 = (-rho*Va^2*c*S_wing/(2*Jy))*C_m_alpha;
    a_theta3 = (-rho*Va^2*c*S_wing/(2*Jy))*C_m_delta_e;
    %a_theta3 *s /s^2+atheta1*s+a_theta2
    T_theta_delta_e = tf([a_theta3 0], [1 a_theta1 a_theta2]);

    %Va/s
    T_h_theta = tf([Va], [1 0]);
    %theta/s
    T_h_Va = tf([theta], [1 0]);

    %TODO hmmmm??
    T_Va_delta_t
    T_Va_theta
    T_v_delta_r
    
end

% returns the derivative of motor thrust with respect to Va
function dThurst = dT_dVa(mav, Va, delta_t)

end

% returns the derivative of motor thrust with respect to delta_t
function dThurst = dT_ddelta_t(mav, Va, delta_t)
end