% % 5.3 Create a function that computes the longitudinal and lateral state
% % space models described in this chapter, linearized around trim.

% %     - Implment "compute_ss_model()" and the dependent functions "df_dx()",
% %       "df_du()", "f_euler()", "euler_state()", "quaternion_state()" in 
% %       "compute_models.py".

function [A_lon,B_lon,A_lat,B_lat] = compute_ss_model(mav, x_trim, u_trim, MAV)
    x_euler = euler_state(trim_state)
    
    A = df_dx(mav, x_euler, trim_input)
    B = df_du(mav, x_euler, trim_input)

    % matrcies tp convert to extract longitudinal states in correct order
    % u, w, q, theta, pd
    E1 = [
    0 0 0 1 0 0 0 0 0 0 0 0;
    0 0 0 0 1 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 1 0;
    0 0 0 0 0 0 0 1 0 0 0 0;
    0 0 -1 0 0 0 0 0 0 0 0 0
    ];
    % delta t and delta e
    E2 = [
    1 0 0 0;
    0 1 0 0
    ];
    % matrcies tp convert to extract lateral states in correct order
    %v, p, r, phi, psi
    E3 = [
    0 0 0 0 1 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 1 0 0;
    0 0 0 0 0 0 0 0 0 0 0 1;
    0 0 0 0 0 0 1 0 0 0 0 0;
    0 0 0 0 0 0 0 0 1 0 0 0
    ];
    %delta_a, delta_r
    E4 = [
    0 0 1 0;
    0 0 0 1
    ];

    %  extract longitudinal states (u, w, q, theta, pd)
    % iput them in the order of original indices
    A_lon = E1*A*E1'
    B_lon = E1*B*E2';
    %  extract lateral states (v, p, r, phi, psi)
    A_lat = E3*A*E3';
    B_lat = E3*B*E4';

end

% convert state x with attitude represented by quaternion
% to x_euler with attitude represented by Euler angles
function x_euler = euler_state(x_quat)
p=x_euler(1:3);
v=x_euler(4:6);
q=x_euler(7:10);
w=x_euler(11:13);
x_quant=[p;v;Euler2Quaternion(q);w];

end

% convert state x_euler with attitude represented by Euler angles
% to x_quat with attitude represented by quaternions
% euler state to quant state
function x_quat = quaternion_state(x_euler)
p=x_euler(1:3);
v=x_euler(4:6);
E=x_euler(7:9);
w=x_euler(10:12);
x_quant=[p;v;Euler2Quaternion(E);w];
end
% return 12x1 dynamics (as if state were Euler state)
% compute f at euler_state
function xdot = f_euler(mav, x_euler, input, MAV)
    %%FIND F(X,U)
    xdot = zeros(12,1);
    % get the variables from x
    pn      = x_euler(1);
    pe      = x_euler(2);
    h       = x_euler(3);
    u       = x_euler(4);
    v       = x_euler(5);
    w       = x_euler(6);
    phi     = x_euler(7);
    theta   = x_euler(8);
    psi     = x_euler(9);
    p       = x_euler(10);
    q       = x_euler(11);
    r       = x_euler(12);
    delta_e = input(1);
    delta_t = input(2);
    delta_a = input(3);
    delta_r = input(4);

    [CD, CL] = drag_and_lift_coiff(mav.alpha);
    % X-axis force coefficients
    CX      = -CD*cos(mav.alpha) + CL*sin(mav.alpha);
    CXq     = -CDq*cos(mav.alpha) + CLq*sin(mav.alpha);
    CXde    = -CDde*cos(mav.alpha) + CLde*sin(mav.alpha);

    % Z-axis force coefficients
    CZ      = -CD*sin(mav.alpha) - CL*cos(mav.alpha);
    CZq     = -CDq*sin(mav.alpha) - CLq*cos(mav.alpha);
    CZde    = -CDde*sin(mav.alpha) - CLde*cos(mav.alpha);

    % get propeller force and torque
    [T_p, Q_p] = prop_force_and_torque(Va, delta_t);

        % Position kinematics
    pndot =  cos(theta)*cos(psi)*u ...
            + (sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi))*v ...
            + (cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi))*w;
    
    pedot =  cos(theta)*sin(psi)*u ...
            + (sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi))*v ...
            + (cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi))*w;

    hdot  =  u*sin(theta) ...
            - v*sin(phi)*cos(theta) ...
            - w*cos(phi)*cos(theta);

    % Translational dynamics
    udot =  r*v - q*w ...
            - MAV.gravity*sin(theta) ...
            + (MAV.rho*Va^2*MAV.S_wing)/(2*m) * ( CX + CXq*(MAV.c*q/(2*Va)) + CXde*delta_e ) ...
            + T_p;

    vdot =  p*w - r*u ...
            + MAV.gravity*cos(theta)*sin(phi) ...
            + (MAV.rho*Va^2*MAV.S_wing)/(2*m) * ( CY0 + CYb*mav.beta + CYp*(MAV.b*p/(2*Va)) ...
            + CYr*(MAV.b*r/(2*Va)) + CYda*delta_a + CYdr*delta_r );

    wdot =  q*u - p*v ...
            + MAV.gravity*cos(theta)*cos(phi) ...
            + (MAV.rho*Va^2*MAV.S_wing)/(2*m) * ( CZ + CZq*(MAV.c*q/(2*Va)) + CZde*delta_e );

    % Euler angle kinematics
    phidot   = p + q*sin(phi)*tan(theta) + r*cos(phi)*tan(theta);
    thetadot = q*cos(phi) - r*sin(phi);
    psidot   = q*sin(phi)/cos(theta) + r*cos(phi)/cos(theta);

    % Rotational dynamics
    pdot =  MAV.Gamma1*p*q - MAV.Gamma2*q*r ...
            + 0.5*MAV.rho*Va^2*MAV.S_wing*MAV.b * ( MAV.C_p_0 + MAV.C_p_beta*mav.beta + MAV.C_p_p*(MAV.b*p/(2*Va)) ...
            + MAV.C_p_r*(MAV.b*r/(2*Va)) + MAV.C_p_delta_a*delta_a + MAV.C_p_delta_r*delta_r ) ...
            -MAV.Gamma3*Q_p;
    qdot = MAV.Gamma5*p*r ...
            - MAV.Gamma6*(p^2 - r^2) ...
            + (MAV.rho*Va^2*MAV.S_wing*MAV.c)/(2*MAV.Jy) * ( MAV.C_m_0 + MAV.C_m_alpha*mav.alpha + MAV.C_m_q*(MAV.c*q/(2*Va)) + C_m_delta_e* delta_e);

    rdot = MAV.Gamma7*p*q ...
            - MAV.Gamma1*q*r ...
            + 0.5*MAV.rho*Va^2*MAV.S_wing*MAV.b * ( MAV.C_r_0 + MAV.C_r_beta*mav.beta ...
            + MAV.C_r_p*(MAV.b*p/(2*Va)) + MAV.C_r_r*(MAV.b*r/(2*Va)) ...
            + MAV.C_r_delta_a*delta_a + MAV.C_r_delta_r*delta_r ) ...
            - Gamma4*Q_p;

    xdot(1,1)   = pndot;
    xdot(2,1)   = pedot;
    xdot(3,1)   = hdot;
    xdot(4,1)   = udot;
    xdot(5,1)   = vdot;
    xdot(6,1)   = wdot;
    xdot(7,1)   = phidot;
    xdot(8,1)   = thetadot;
    xdot(9,1)   = psidot;
    xdot(10,1)  = pdot;
    xdot(11,1)  = qdot;
    xdot(12,1)  = rdot;

end

% take partial of f_euler with respect to x_euler
%get Jacobian of f(x,u) wrt state variables x
function A = df_dx(mav, x_euler, input, MAV)
    A = zeros(12, 12);

    % get f at x_euler
    fx = f_euler(mav, x_euler, input, MAV)

    % go column by column computing at each x
    for i =1:12    
        %compute x+eps
        x_eps(i,0) = x_euler(i,0) + eps;
        % compute dynamics at x+eps
        fx_eps = f_euler(mav, x_eps, input, MAV);
        % find df/dx
        df_dxi(fx_eps-fx)/eps;
        % place into ith column of A
        A(:,i) = df_dxi[:,0];
    end
end

% take partial of f_euler with respect to delta
function B = df_du(mav, x_euler, delta, MAV)
    B = zeros(12, 4)  # Jacobian of f wrt u

    % get f at u
    fu = f_euler(mav, x_euler, input, MAV)

    % go column by column computing at each f
    for i =1:4 
        %compute u+eps
        u_eps(i,0) = input(i,0) + eps;
        % compute dynamics at x+eps
        fu_eps = f_euler(mav, x_euler, u_eps, MAV);
        % find df/dx
        df_dui(fu_eps-fu)/eps;
        % place into ith column of A
        B(:,i) = df_dui[:,0];
    end
end





% function [A_lon,B_lon,A_lat,B_lat] = compute_ss_model(mav, x_trim, u_trim, MAV)
    
% end

% % convert state x with attitude represented by quaternion
% % to x_euler with attitude represented by Euler angles
% function x_euler = euler_state(x_quat)
% end

% % convert state x_euler with attitude represented by Euler angles
% % to x_quat with attitude represented by quaternions
% function x_quat = quaternion_state(x_euler)
% end

% % return 12x1 dynamics (as if state were Euler state)
% % compute f at euler_state
% function xdot = f_euler(mav, x_euler, input, MAV)
%             pn    = state(1);
%             pe    = state(2);
%             pe    = state(3);
%             u     = state(4);
%             v     = state(5);
%             w     = state(6);
%             e0    = state(7);
%             e1    = state(8);
%             e2    = state(9);
%             e3    = state(10);
%             p     = state(11);
%             q     = state(12);
%             r     = state(13);
% end

% % take partial of f_euler with respect to x_euler
% function A = df_dx(mav, x_euler, input, MAV)
%     % where does C_X_alpha come from
%     Xu = (u_star*rho*S_wing/m)*(C_X_0 + C_X_alpha*P.alpha_star + C_X_delta_e*delta_e_star) ...
%         - (rho*S_wing*w_star*C_X_alpha)/(2*m) ...
%         + (rho*S_wing*C_X_q*u_star*q_star)/(4*m*Va_star) ...
%         + dTp_du;

%     % Xw
%     Xw = -q_star ...
%         + (w_star*rho*S_wing/m)*(C_X_0 + C_X_alpha*P.alpha_star + C_X_delta_e*delta_e_star) ...
%         + (rho*S_wing*C_X_alpha*w_star*q_star)/(4*m*Va_star) ...
%         + (rho*S_wing*C_X_alpha*u_star)/(2*m) ...
%         + dTp_dw;

%     % Xq
%     Xq = -w_star + (rho*Va_star*S_wing*C_X_q*c)/(4*m);

%     % X_delta_e
%     X_de = (rho*Va_star^2*S_wing*C_X_delta_e)/(2*m);

%     % X_delta_t
%     X_delta_t = dTp_ddt;

%     % Zu
%     Zu = q_star ...
%         + (u_star*rho*S_wing/m)*(C_Z_0 + C_Z_alpha*P.alpha_star + C_Z_delta_e*delta_e_star) ...
%         - (rho*S_wing*C_Z_alpha*w_star)/(2*m) ...
%         + (u_star*rho*S_wing*C_Z_q*q_star)/(4*m*Va_star);

%     % Zw
%     Zw = (w_star*rho*S_wing/m)*(C_Z_0 + C_Z_alpha*P.alpha_star + C_Z_delta_e*delta_e_star) ...
%         + (rho*S_wing*C_Z_alpha*u_star)/(2*m) ...
%         + (rho*w_star*S_wing*C_Z_q*q_star)/(4*m*Va_star);

%     % Zq
%     Zq = u_star + (rho*Va_star*S_wing*C_Z_q*c)/(4*m);

%     % Z_delta_e
%     Z_delta_e = (rho*Va_star^2*S_wing*C_Z_delta_e)/(2*m);

%     % Mu
%     Mu = (u_star*rho*S_wing*c/Jy)*(C_m_0 + C_m_alpha*P.alpha_star + Cmde*delta_e_star) ...
%         - (rho*S_wing*c*C_m_alpha*w_star)/(2*Jy) ...
%         + (rho*S_wing*c^2*C_m_q*q_star*u_star)/(4*Jy*Va_star);

%     % Mw
%     Mw = (w_star*rho*S_wing*c/Jy)*(C_m_0 + C_m_alpha*P.alpha_star + Cmde*delta_e_star) ...
%         + (rho*S_wing*c*C_m_alpha*u_star)/(2*Jy) ...
%         + (rho*S_wing*c^2*C_m_q*q_star*w_star)/(4*Jy*Va_star);

%     % Mq
%     Mq = (rho*Va_star*S_wing*c^2*C_m_q)/(4*Jy);


%     % Yv
%     Yv = (rho*S_wing*b*v_star)/(4*m*Va_star)*(C_Y_p*p_star + C_Y_r*r_star) ...
%         + (rho*S_wing*C_Y_beta)/(2*m)*sqrt(u_star^2 + w_star^2) ...
%         + (rho*S_wing*v_star/m)*(CY0 + C_Y_beta*beta_star + C_Y_delta_a*delta_a_star + C_Y_delta_r*delta_r_star);

%     % Yp
%     Yp = w_star + (rho*Va_star*S_wing*b)/(4*m)*C_Y_p;

%     % Yr
%     Yr = -u_star + (rho*Va_star*S_wing*b)/(4*m)*C_Y_r;

%     % Y_delta_a
%     Y_delta_a = (rho*Va_star^2*S_wing)/(2*m)*C_Y_delta_a;

%     % Y_delta_r
%     Y_delta_r = (rho*Va_star^2*S_wing)/(2*m)*C_Y_delta_r;

%     % Lv
%     Lv = (rho*S_wing*b^2*v_star)/(4*Va_star)*(C_p_p*p_star + C_p_r*r_star) ...
%         + (rho*S_wing*b*C_p_b)/2*sqrt(u_star^2 + w_star^2) ...
%         + rho*S_wing*b*v_star*(C_p_0 + C_p_b*beta_star + C_p_delta_a*delta_a_star + C_p_delta_r*delta_r_star);

%     % Lp
%     Lp = Gamma1*q_star + (rho*Va_star*S_wing*b^2)/4*C_p_p;

%     % Lr
%     Lr = -Gamma2*q_star + (rho*Va_star*S_wing*b^2)/4*C_p_r;

%     % L_delta_a
%     L_da = (rho*Va_star^2*S_wing*b)/2*C_p_delta_a;

%     % L_delta_r
%     L_dr = (rho*Va_star^2*S_wing*b)/2*C_p_delta_r;

%     % Nv
%     Nv = (rho*S_wing*b^2*v_star)/(4*Va_star)*(C_r_p*p_star + C_r_r*r_star) ...
%         + (rho*S_wing*b*Crb)/2*sqrt(u_star^2 + w_star^2) ...
%         + rho*S_wing*b*v_star*(Cr0 + Crb*beta_star + C_r_delta_a*delta_a_star + C_r_delta_r*delta_r_star);

%     % Np
%     Np = Gamma7*q_star + (rho*Va_star*S_wing*b^2)/4*C_r_p;

%     % Nr
%     Nr = -Gamma1*q_star + (rho*Va_star*S_wing*b^2)/4*C_r_r;

%     % N_delta_a
%     N_delta_a = (rho*Va_star^2*S_wing*b)/2*C_r_delta_a;

%     % N_delta_r
%     N_delta_r = (rho*Va_star^2*S_wing*b)/2*C_r_delta_r;

%     % A matrix terms
%     A14 = P.g * cos(theta_star) * cos(phi_star);
%     A43 = cos(phi_star) * tan(theta_star);
%     A44 = q_star * cos(phi_star) * tan(theta_star) - r_star * sin(phi_star) * tan(theta_star);
%     A53 = cos(phi_star) / cos(theta_star);
%     A54 = p_star * cos(phi_star) / cos(theta_star) - r_star * sin(phi_star) / cos(theta_star);
%     % M_delta_e
%     M_delta_e = (rho*Va_star^2*S*c*Cmde)/(2*Iy);
%     A = [...
%         X_u X_w X_q -g*cos(theta) 0; %u
%         Y_v Y_p Y_r A_14 0; %v
%         Z_u Z_w Z_q -g*sin(theta) 0; %w
%         L_v L_p L_r 0 0;    %p
%         M_u M_w M_q 0 0;    %q
%         N_v N_p N_r 0 0;    %r
%         0 0 1 0 0;          %theta
%         0 1 A_43 A_44 0;    %phi
%         0 0 A_53 A_54 0;    %psi
%         sin(theta) -cos(theta) 0 u_star*cos(theta_star)+w_star*sin(theta_star) 0; %h
%     ];

% end

% % take partial of f_euler with respect to delta
% function B = df_du(mav, x_euler, delta, MAV)
%     B = [
%         0 0 X_del_e X_del_t;%u
%         0 0 Y_del_a Y_del_r;%v
%         Z_del_e 0 0 0;      %w
%         L_del_a L_del_r 0 0;%p
%         M_del_a 0 0 0;      %q
%         N_del_a N_del_r 0 0;%r
%         0 0 0 0;            %theta
%         0 0 0 0;            %phi
%         0 0 0 0;            %psi
%         0 0 0 0;            %h
%     ];
% end