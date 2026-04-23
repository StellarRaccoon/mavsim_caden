% % 5.3 Create a function that computes the longitudinal and lateral state
        % % space models described in this chapter, linearized around trim.

% %     - Implment "compute_ss_model()" and the dependent functions "df_dx()",
        % %       "df_du()", "f_euler()", "euler_state()", "quaternion_state()" in 
        % %       "compute_models.py".

function [A_lon,B_lon,A_lat,B_lat] = compute_ss_model(mav, x_trim, u_trim, MAV)
        x_euler = euler_state(x_trim);

        A = df_dx(mav, x_euler, u_trim, MAV);
        B = df_du(mav, x_euler, u_trim, MAV);

        % matrcies tp convert to extract longitudinal states in correct order
        % u, w, q, theta, pd
        E1 = [
        0 0 0 1 0 0 0 0 0 0 0 0;
        0 0 0 0 0 1 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 1 0;
        0 0 0 0 0 0 0 1 0 0 0 0;
        0 0 -1 0 0 0 0 0 0 0 0 0;
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
        A_lon = E1*A*E1';
        B_lon = E1*B*E2';
        %  extract lateral states (v, p, r, phi, psi)
        A_lat = E3*A*E3';
        B_lat = E3*B*E4';

end

        % convert state x with attitude represented by quaternion
function x_euler = euler_state(x_quat)
        p=x_quat(1:3);
        v=x_quat(4:6);
        q=x_quat(7:10);
        w=x_quat(11:13);
        e=Quaternion2Euler(q);
        x_euler=[p;v;e;w];

end

        % convert state from euler to quant
function x_quat = quaternion_state(x_euler)
        p=x_euler(1:3);
        v=x_euler(4:6);
        E=x_euler(7:9);
        w=x_euler(10:12);
        phi=E(1);
        theta = E(2);
        psi= E(3);
        
        Eq= Euler2Quaternion(phi, theta, psi);
        x_quat=[p;v;Eq(1); Eq(2);Eq(3); Eq(4);w];
end
        % get euler xdot dynamics 
        % compute f at euler_state
function xdot = f_euler(mav, x_euler, input, MAV)

        mav.state = quaternion_state(x_euler);
        mav.update_velocity_data(zeros(6,1));

        forces_moments = mav.forces_moments(input, MAV);
        xdot = mav.derivatives(mav.state, forces_moments, MAV);
        xdot = euler_state(xdot);
end

        % take partial of f_euler with respect to x_euler
        %get Jacobian of f(x,u) wrt state variables x
function A = df_dx(mav, x_euler, input, MAV)
        A = zeros(12, 12);

        % get f at x_euler
        fx = f_euler(mav, x_euler, input, MAV)
        eps= 10^-8;
        % go column by row computing at each x
        for i =1:12    
        %compute x+eps
        x_eps = x_euler;
        x_eps(i,1) = x_euler(i,1) + eps;
        % compute dynamics at x+eps
        fx_eps = f_euler(mav, x_eps, input, MAV);
        % find df/dx
        df_dxi=(fx_eps-fx)/eps;
        % place into ith column of A
        A(:,i) = df_dxi(:,1);
        end

end

        % take partial of f_euler with respect to delta
function B = df_du(mav, x_euler, input,  MAV)
        B = zeros(12, 4);  # Jacobian of f wrt u

        % get f at u
        fu = f_euler(mav, x_euler, input, MAV);
        eps= 10^-8;
        % go column by column computing at each f
        for i =1:4 
        u_eps = input;
        %compute u+eps
        u_eps(i,1) = input(i,1) + eps;
        % compute dynamics at x+eps
        fu_eps = f_euler(mav, x_euler, u_eps, MAV);
        % find df/dx
        df_dui = (fu_eps-fu)/eps;
        % place into ith column of A
        B(:,i) = df_dui(:,1);
        end
end
