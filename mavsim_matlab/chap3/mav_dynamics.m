% mav dynamics - implement rigid body dynamics for mav
%
% mavsim
%     - Beard & McLain, PUP, 2012
%     - Update history:  
%         12/19/2018 - RWB
classdef mav_dynamics < handle
    %--------------------------------
    properties
        ts_simulation
        state
        true_state
    end
    %--------------------------------
    methods
        %------constructor-----------
        function self = mav_dynamics(Ts, MAV)
            self.ts_simulation = Ts; % time step between function calls
            self.state = [MAV.pn0; MAV.pe0; MAV.pd0; MAV.u0; MAV.v0; MAV.w0;...
                MAV.e0; MAV.e1; MAV.e2; MAV.e3; MAV.p0; MAV.q0; MAV.r0];
            addpath('../message_types'); self.true_state = msg_state();
        end
        %---------------------------
        function self=update_state(self, forces_moments, MAV)
            %
            % Integrate the differential equations defining dynamics
            % forces_moments are the forces and moments on the MAV.
            % 
            % Integrate ODE using Runge-Kutta RK4 algorithm
            k1 = self.derivatives(self.state, forces_moments, MAV);
            k2 = self.derivatives(self.state + self.ts_simulation/2*k1, forces_moments, MAV);
            k3 = self.derivatives(self.state + self.ts_simulation/2*k2, forces_moments, MAV);
            k4 = self.derivatives(self.state + self.ts_simulation*k3, forces_moments, MAV);
            self.state = self.state + self.ts_simulation/6 * (k1 + 2*k2 + 2*k3 + k4);
            
            % normalize the quaternion
            self.state(7:10) = self.state(7:10)/norm(self.state(7:10));
            self.update_true_state();
        end
        %----------------------------
        function xdot = derivatives(self, state, forces_moments, MAV)
            pn    = state(1);
            pe    = state(2);
            pe    = state(3);
            u     = state(4);
            v     = state(5);
            w     = state(6);
            e0    = state(7);
            e1    = state(8);
            e2    = state(9);
            e3    = state(10);
            p     = state(11);
            q     = state(12);
            r     = state(13);
            fx    = forces_moments(1);
            fy    = forces_moments(2);
            fz    = forces_moments(3);
            ell   = forces_moments(4);
            m     = forces_moments(5);
            n     = forces_moments(6);
            
            % quanternion rotation matrix

            % position kinematics
            pos_dot = Quaternion2Rotation([e0, e1, e2, e3]) * [u; v; w];
            pn_dot = pos_dot(1);
            pe_dot = pos_dot(2);
            pd_dot = pos_dot(3);
            fprintf("pn_dot: %d pe_dot: %d pd_dot: %d\n", pn_dot, pe_dot, pd_dot);

            % position dynamics
            vel_dot = [r*v-q*w; p*w-r*u; q*u-p*v;] + m^(-1)*[fx; fy; fz]
            u_dot = vel_dot(1); 
            v_dot = vel_dot(2);
            w_dot = vel_dot(3);
            fprintf("u_dot: %d v_dot: %d w_dot: %d\n", u_dot, v_dot, w_dot);
            % rotational kinematics
            temp = [
                    0, -p, -q, -r;
                    p, 0, r, -q;
                    q, -r, 0, p;
                    r, q, -p, 0
                    ]
            e_dot = temp * [e0; e1; e2; e3];
            e0_dot = e_dot(1);
            e1_dot = e_dot(2);
            e2_dot = e_dot(3);
            e3_dot = e_dot(4);
            % rotational dynamics
            fprintf('Js\n')
            MAV.Jx
            MAV.Jy
            MAV.Jz
            MAV.Jxz
            p_dot = 0;
            q_dot = 0;
            r_dot = 0;
        
            % collect all the derivaties of the states
            xdot = [pn_dot; pe_dot; pd_dot; u_dot; v_dot; w_dot;...
                    e0_dot; e1_dot; e2_dot; e3_dot; p_dot; q_dot; r_dot];
        end
        %----------------------------
        function self=update_true_state(self)
            [phi, theta, psi] = Quaternion2Euler(self.state(7:10));
            self.true_state.pn = self.state(1);  % pn
            self.true_state.pe = self.state(2);  % pd
            self.true_state.h = -self.state(3);  % h
            self.true_state.phi = phi; % phi
            self.true_state.theta = theta; % theta
            self.true_state.psi = psi; % psi
            self.true_state.p = self.state(11); % p
            self.true_state.q = self.state(12); % q
            self.true_state.r = self.state(13); % r
        end
    end
end