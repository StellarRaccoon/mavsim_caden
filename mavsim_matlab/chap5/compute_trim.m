%  MAV.Create a function that computes the trim state and the trim
% inputs for a desired airspeed of Va and a desired flight path
% angle of gamma. Set the initial condition in your simulation to the
% trim state and trim input, and verify that the airMAV.craft maintains
% trim until numerical errors cause it to drift

%unknowns are angle of attack, pitch angle theta, elevator deflection delta e, throttle delta e
function [x_trim, u_trim] = compute_trim(mav, Va, gamma, MAV)
    gamma
    alpha0 = 5*pi/180;


    theta = alpha0;
    e = Euler2Quaternion(0, theta,0);
    state0 = [
        0; % pn
        0; % pe
        mav.state(3); % pd
        Va*cos(theta);
        0; % v
        Va*sin(theta);; % w
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
lb_delta = [-1; 0; -1; -1];   % throttle >= 0
ub_delta = [ 1; 1;  1;  1];

lb = [lb_x; lb_delta];
ub = [ub_x; ub_delta];

    x0 = [ state0; delta0 ];

    to =  @(xu) trim_objective(xu, mav, Va, gamma, MAV);
    tc = @(xu) trim_constraints(xu, mav, Va, gamma, MAV);
options = optimset('Display','iter','MaxIter',1000,'TolFun',1e-9,'TolX',1e-9);
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

% % nonlinear constraints for trim optimization
% function [c, ceq] = trim_constraints(xu, mav, Va, gamma, MAV)
%     ceq = [                             
%         xu(4)^2 + xu(5)^2 + xu(6)^2 - Va^2,  # Va=Vg
%         xu(5),  # v=0
%         xu(7)^2 + xu(8)^2 + xu(9)^2 + xu(10)^2 - 1,  # Unit length quant
%         xu(8),  # e1=0
%         xu(10),  # e3=0
%         xu(11),  # p=0 
%         xu(12),  # q=0
%         xu(13)  # r=0
%         ];
%     c =[];

% end
function [c, ceq] = trim_constraints(xu, mav, Va, gamma, MAV)
    x     = xu(1:13);
    delta = xu(14:17);

    mav.state = x;
    mav.update_velocity_data(zeros(6,1));
    fm   = mav.forces_moments(delta, MAV);
    xdot = mav.derivatives(x, fm, MAV);
mav.state = x;
mav.update_velocity_data(zeros(6,1));
fm   = mav.forces_moments(delta, MAV);
xdot = mav.derivatives(x, fm, MAV);


    ceq = [
        x(4)^2 + x(5)^2 + x(6)^2 - Va^2;
        x(5);
        x(7)^2 + x(8)^2 + x(9)^2 + x(10)^2 - 1;
        x(8);
        x(10);
        x(11);
        x(12);
        x(13);
        % xdot(3) + Va*sin(gamma);   % enforce vertical trim
    ];
    c = [];
    ceq(end+1) = xdot(3) + Va*sin(gamma);   % enforce vertical trim
end
