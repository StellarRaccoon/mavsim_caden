%  MAV.Create a function that computes the trim state and the trim
% inputs for a desired airspeed of ∗ and a desired flight path
% angle of ∗. Set the initial condition in your simulation to the
% trim state and trim input, and verify that the airMAV.craft maintains
% trim until numerical errors cause it to drift

%unknowns are angle of attack, pitch angle theta, elevator deflection delta e, throttle delta e
function [x_trim, u_trim] = compute_trim(mav, Va, gamma, MAV)
    % Va is the desired airspeed (m/s)
    % gamma is the desired flight path angle (radians)
    % R is the desired radius (m) - use (+) for right handed orbit, 
%                                       (-) for left handed orbit

    % define initial state and input
    addpath('../tools');
    state0 = [
        MAV.pn0;
        MAV.pe0;
        MAV.pd0;
        MAV.u0;
        MAV.v0;
        MAV.w0;
        MAV.phi0;
        MAV.theta0;
        MAV.psi0;
        MAV.p0;
        MAV.q0;
        MAV.r0
    ];

    delta_e = 0;
    delta_t = 0.67;
    delta_a = 0;  
    delta_r = 0;
    delta0 = [delta_e; delta_t; delta_a; delta_r];

    x0 = [ state0; delta0 ];
    xstar = fmincon(@trim_objective, x0, [], [],...
                    [], [], [], [], @trim_constraints, [],...
                    mav, Va, gamma, MAV);
    x_trim = xstar(1:13);
    u_trim = xstar(14:17);
    J = trim_objective(xstar, mav, Va, gamma, MAV);
end

% objective function to be minimized f(x,u) where x is the states and u is the control surface deltas
function J = trim_objective(x, mav, Va, gamma, MAV)
    x = x(1:13);
    u = x(14:17);
    
    % get the nonlinear xdot
    forces_moments = mav.forces_moments(u, MAV);
    xdot = mav.derivatives(mav, x, forces_moments, MAV)
    % desired trim
    xdot_star = [
        0.1;
        0.1;
        Va * sin(gamma);
        0;
        0;
        0;
        0;
        0;
        0;
        0;
        0;
        0
    ];
    % we want 0=f(x,u)-xdot* where f(x,u) = xdot and f(x*,u*) = xdot*
    r = xdot - xdot_star;
    % Squared norm - a scalar
    J = r'*r;


end

% nonlinear constraints for trim optimization
function [c, ceq] = trim_constraints(x, mav, Va, gamma, MAV)
u       = x(4);
v       = x(5);
w       = x(6);

% Va must be constant
Va_current = sqrt(u^2+v^2+w^2);
%gamma must be constant
gamma_current = atan2(-w,sqrt(u^2+v^2));
ceq = [Va_current - Va; gamma_current - gamma];
c = [];
end



% %%FIND F(X,U)
% % get the variables from x
% pn      = x(1);
% pe      = x(2);
% h       = x(3);
% u       = x(4);
% v       = x(5);
% w       = x(6);
% phi     = x(7);
% theta   = x(8);
% psi     = x(9);
% p       = x(10);
% q       = x(11);
% r       = x(12);
% delta_e = u(1);
% delta_t = u(2);
% delta_a = u(3);
% delta_r = u(4);

% [CD, CL] = drag_and_lift_coiff(mav.alpha);
% % X-axis force coefficients
% CX      = -CD*cos(mav.alpha) + CL*sin(mav.alpha);
% CXq     = -CDq*cos(mav.alpha) + CLq*sin(mav.alpha);
% CXde    = -CDde*cos(mav.alpha) + CLde*sin(mav.alpha);

% % Z-axis force coefficients
% CZ      = -CD*sin(mav.alpha) - CL*cos(mav.alpha);
% CZq     = -CDq*sin(mav.alpha) - CLq*cos(mav.alpha);
% CZde    = -CDde*sin(mav.alpha) - CLde*cos(mav.alpha);

% % get propeller force and torque
% [T_p, Q_p] = prop_force_and_torque(Va, delta_t);

%     % Position kinematics
% pndot =  cos(theta)*cos(psi)*u ...
%         + (sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi))*v ...
%         + (cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi))*w;

% pedot =  cos(theta)*sin(psi)*u ...
%         + (sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi))*v ...
%         + (cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi))*w;

% hdot  =  u*sin(theta) ...
%         - v*sin(phi)*cos(theta) ...
%         - w*cos(phi)*cos(theta);

% % Translational dynamics
% udot =  r*v - q*w ...
%         - MAV.gravity*sin(theta) ...
%         + (MAV.rho*Va^2*MAV.S_wing)/(2*m) * ( CX + CXq*(MAV.c*q/(2*Va)) + CXde*delta_e ) ...
%         + T_p;

% vdot =  p*w - r*u ...
%         + MAV.gravity*cos(theta)*sin(phi) ...
%         + (MAV.rho*Va^2*MAV.S_wing)/(2*m) * ( CY0 + CYb*mav.beta + CYp*(MAV.b*p/(2*Va)) ...
%         + CYr*(MAV.b*r/(2*Va)) + CYda*delta_a + CYdr*delta_r );

% wdot =  q*u - p*v ...
%         + MAV.gravity*cos(theta)*cos(phi) ...
%         + (MAV.rho*Va^2*MAV.S_wing)/(2*m) * ( CZ + CZq*(MAV.c*q/(2*Va)) + CZde*delta_e );

% % Euler angle kinematics
% phidot   = p + q*sin(phi)*tan(theta) + r*cos(phi)*tan(theta);
% thetadot = q*cos(phi) - r*sin(phi);
% psidot   = q*sin(phi)/cos(theta) + r*cos(phi)/cos(theta);

% % Rotational dynamics
% pdot =  MAV.Gamma1*p*q - MAV.Gamma2*q*r ...
%         + 0.5*MAV.rho*Va^2*MAV.S_wing*MAV.b * ( MAV.C_p_0 + MAV.C_p_beta*mav.beta + MAV.C_p_p*(MAV.b*p/(2*Va)) ...
%         + MAV.C_p_r*(MAV.b*r/(2*Va)) + MAV.C_p_delta_a*delta_a + MAV.C_p_delta_r*delta_r ) ...
%         -MAV.Gamma3*Q_p;
% qdot = MAV.Gamma5*p*r ...
%         - MAV.Gamma6*(p^2 - r^2) ...
%         + (MAV.rho*Va^2*MAV.S_wing*MAV.c)/(2*MAV.Jy) * ( MAV.C_m_0 + MAV.C_m_alpha*mav.alpha + MAV.C_m_q*(MAV.c*q/(2*Va)) + C_m_delta_e* delta_e);

% rdot = MAV.Gamma7*p*q ...
%         - MAV.Gamma1*q*r ...
%         + 0.5*MAV.rho*Va^2*MAV.S_wing*MAV.b * ( MAV.C_r_0 + MAV.C_r_beta*mav.beta ...
%         + MAV.C_r_p*(MAV.b*p/(2*Va)) + MAV.C_r_r*(MAV.b*r/(2*Va)) ...
%         + MAV.C_r_delta_a*delta_a + MAV.C_r_delta_r*delta_r ) ...
%         - Gamma4*Q_p;
