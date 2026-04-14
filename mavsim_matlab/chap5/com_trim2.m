%  MAV.Create a function that computes the trim state and the trim
% inputs for a desired airspeed of Va and a desired flight path
% angle of gamma. Set the initial condition in your simulation to the
% trim state and trim input, and verify that the airMAV.craft maintains
% trim until numerical errors cause it to drift

%unknowns are angle of attack, pitch angle theta, elevator deflection delta e, throttle delta e
function [x_trim, u_trim] = compute_trim(mav, Va, gamma, MAV)
    gamma
    e = Euler2Quaternion(0, gamma,0);
    state0 = [
        0; % pn
        0; % pe
        mav.state(3); % pd
        Va; % u
        0; % v
        0; % w
        e(1); % e0
        e(2); % e1
        e(3); % e2
        e(4); % e3
        0; % p
        0; % q
        0 % r
    ];

    delta0 = [
            0; % delta_e 
            0.5; % delta_t
            0; % delta_a
            0  % delta_r
        ];
        % bounds for states (13 states)
lb_x = -inf(13,1);
ub_x =  inf(13,1);

% bounds for controls (match Python)
lb_delta = -1.0 * ones(4,1);
ub_delta =  1.0 * ones(4,1);

lb = [lb_x; lb_delta];
ub = [ub_x; ub_delta];

    x0 = [ state0; delta0 ];

    to =  @(xu) trim_objective(xu, mav, Va, gamma, MAV);
    tc = @(xu) trim_constraints(xu, mav, Va, gamma, MAV);
    xstar = fmincon(to, x0, [], [],...
                    [], [], lb, ub, tc);
    % xstar = fmincon(to, x0, [], [],...
    %                 [], [], [], [], tc);
    [c, ceq] = trim_constraints(xstar, mav, Va, gamma, MAV);
    disp('J:');
    % disp(max(abs(ceq)));  % Be near 0
    fprintf("J\n");
    % J = trim_objective(xstar, mav, Va, gamma, MAV);
    x_trim = xstar(1:13);
    u_trim = xstar(14:17);
end

% objective function to be minimized f(x,u) where x is the states and u is the control surface deltas
%we want all the derivatives to be 0
function J = trim_objective(xu, mav, Va, gamma, MAV)
    % x     = xu(1:13);
    % delta = xu(14:17);

    % % small controls, small pitch
    % theta = quat2euler(x(7:10));  % or your own quaternion→Euler

    % w_delta = 1.0;
    % w_theta = 0.1;

    % J = w_delta * (delta.'*delta) + w_theta * theta^2;
    x = xu(1:13);
    delta = xu(14:17); %extract the deltas

    %goal derivatives are at 0
    xdot_target = [0;0;-Va*sin(gamma);0;0;0;0;0;0;0;0;0;0];
    
    %find the forces and moments needed given the current delta
    mav.state = x;
    mav.update_velocity_data(zeros(6,1)); %set wind as 0
    forces_moments = mav.forces_moments(delta, MAV);
    xdot = mav.derivatives(x, forces_moments, MAV);

    error = xdot_target-xdot;

    % J = norm(error(3:13).^2);
    J = error(3:13).' * error(3:13);   % ||error(3:13)||^2

end

% nonlinear constraints for trim optimization
function [c, ceq] = trim_constraints(xu, mav, Va, gamma, MAV)
    x = xu(1:13);
    delta = xu(14:17); %extract the deltas
    mav.state = x;
    fm = mav.forces_moments(delta, MAV);
xdot = mav.derivatives(x, fm, MAV);
    ceq = [                             
        xu(4)^2 + xu(5)^2 + xu(6)^2 - Va^2,  # Va=Vg
        xu(5),  # v=0
        xu(7)^2 + xu(8)^2 + xu(9)^2 + xu(10)^2 - 1,  # Unit length quant
        xu(8),  # e1=0
        xu(10),  # e3=0
        xu(11),  # p=0 
        xu(12),  # q=0
        xu(13),  # r=0
    %         xdot(3) + Va*sin(gamma);   % pd_dot = -Va*sin(gamma)
    % xdot(4);                   % u_dot = 0
    % xdot(5);                   % v_dot = 0
    % xdot(6);                   % w_dot = 0
    % xdot(11);                  % p_dot = 0
    % xdot(12);                  % q_dot = 0
    % xdot(13)                % r_dot = 0
        ];
    c =[];

end
% function [x_trim, u_trim] = compute_trim(mav, Va, gamma, MAV)

%     theta0 = gamma;
%     e = Euler2Quaternion(0, theta0, 0);

%     x0 = [
%         0; 0; mav.state(3);
%         Va; 0; 0;
%         e(1); e(2); e(3); e(4);
%         0; 0; 0
%     ];

%     delta0 = [0.5; 1; 0; 0];

%     xu0 = [x0; delta0];

%     lb = [-inf(13,1); -1*ones(4,1)];
%     ub = [ inf(13,1);  1*ones(4,1)];

%     options = optimset("Display","iter","MaxIter",5000,"MaxFunEvals",5000);

%     to = @(xu) trim_objective(xu, mav, Va, gamma, MAV);
%     tc = @(xu) trim_constraints(xu, mav, Va, gamma, MAV);

%     xu_star = fmincon(to, xu0, [], [], [], [], lb, ub, tc, options);

%     x_trim = xu_star(1:13);
%     u_trim = xu_star(14:17);
% end
% function [c, ceq] = trim_constraints(xu, mav, Va, gamma, MAV)

%     x = xu(1:13);

%     u = x(4); v = x(5); w = x(6);

%     ceq = [
%         u^2 + v^2 + w^2 - Va^2;    % Va constraint
%         v;                         % no sideslip
%         x(7)^2 + x(8)^2 + x(9)^2 + x(10)^2 - 1;  % unit quaternion
%         x(8);                      % e1 = 0
%         x(10);                     % e3 = 0
%         x(11);                     % p = 0
%         x(12);                     % q = 0
%         x(13);                     % r = 0
%     ];

%     c = [];
% end
% function J = trim_objective(xu, mav, Va, gamma, MAV)

%     x     = xu(1:13);
%     delta = xu(14:17);

%     mav.state = x;
%     mav.update_velocity_data(zeros(6,1));
%     fm   = mav.forces_moments(delta, MAV);
%     xdot = mav.derivatives(x, fm, MAV);

%     % desired xdot
%     xdot_target = [0;0;-Va*sin(gamma);0;0;0;0;0;0;0;0;0;0];

%     err = xdot - xdot_target;

%     % minimize dynamics residual (same as Python)
%     J = err(3:13).' * err(3:13);
% end
