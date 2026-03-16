% mav dynamics - implement rigid body dynamics for mav
%
% mavMatSim 
%     - Beard & McLain, PUP, 2012
%     - Update history:  
%         1/18/2019 - RWB
classdef mav_dynamics < handle
%--------------------------------
    properties
        ts_simulation
        state
        Va
        alpha
        beta
        wind
        true_state
    end
    %--------------------------------
    methods
        %------constructor-----------
        function self = mav_dynamics(Ts, MAV)
            self.ts_simulation = Ts; % time step between function calls
            self.state = [MAV.pn0; MAV.pe0; MAV.pd0; MAV.u0; MAV.v0; MAV.w0;...
                MAV.e0; MAV.e1; MAV.e2; MAV.e3; MAV.p0; MAV.q0; MAV.r0];
            self.Va = 0;
            self.alpha = 0; 
            self.beta = 0;
            self.wind = [0;0;0];
            addpath('../message_types'); self.true_state = msg_state();
        end
        %---------------------------
        function self=update_state(self, delta, wind, MAV)
            %
            % Integrate the differential equations defining dynamics
            % forces_moments are the forces and moments on the MAV.
            % 
            
            % get forces and moments acting on rigid body
            forces_moments = self.forces_moments(delta, MAV);
            
            % Integrate ODE using Runge-Kutta RK4 algorithm
            k1 = self.derivatives(self.state, forces_moments, MAV);
            k2 = self.derivatives(self.state + self.ts_simulation/2*k1, forces_moments, MAV);
            k3 = self.derivatives(self.state + self.ts_simulation/2*k2, forces_moments, MAV);
            k4 = self.derivatives(self.state + self.ts_simulation*k3, forces_moments, MAV);
            self.state = self.state + self.ts_simulation/6 * (k1 + 2*k2 + 2*k3 + k4);
            
            % normalize the quaternion
            self.state(7:10) = self.state(7:10)/norm(self.state(7:10));
            
            % update the airspeed, angle of attack, and side slip angles
            self.update_velocity_data(wind);
            
            % update the message class for the true state
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
            %fprintf("pn_dot: %d pe_dot: %d pd_dot: %d\n", pn_dot, pe_dot, pd_dot);

            % position dynamics
            vel_dot = [r*v-q*w; p*w-r*u; q*u-p*v;] + MAV.mass^(-1)*[fx; fy; fz];
            u_dot = vel_dot(1); 
            v_dot = vel_dot(2);
            w_dot = vel_dot(3);
            %fprintf("u_dot: %d v_dot: %d w_dot: %d\n", u_dot, v_dot, w_dot);
            % rotational kinematics
            temp = [
                    0, -p, -q, -r;
                    p, 0, r, -q;
                    q, -r, 0, p;
                    r, q, -p, 0
                    ];
            e_dot = temp * [e0; e1; e2; e3];
            e0_dot = e_dot(1);
            e1_dot = e_dot(2);
            e2_dot = e_dot(3);
            e3_dot = e_dot(4);
            % rotational dynamics
            
            Gamma1 = MAV.Gamma1;
            Gamma2 = MAV.Gamma2;
            Gamma3 = MAV.Gamma3;
            Gamma4 = MAV.Gamma4;
            Gamma5 = MAV.Gamma5;
            Gamma6 = MAV.Gamma6;
            Gamma7 = MAV.Gamma7;
            Gamma8 = MAV.Gamma8;
            Jy = MAV.Jy;
            
            rate_dot= [Gamma1*p*q-Gamma2*q*r; Gamma5*p*r-Gamma6*(p^2-r^2); Gamma7*p*q-Gamma1*q*r]+...
            [Gamma3, 0, Gamma4; 0, 1/Jy, 0; Gamma4, 0, Gamma8;]* [ell; m; n];
            p_dot = rate_dot(1);
            q_dot = rate_dot(2);
            r_dot = rate_dot(3);
            % collect all the derivaties of the states
            xdot = [pn_dot; pe_dot; pd_dot; u_dot; v_dot; w_dot;...
                    e0_dot; e1_dot; e2_dot; e3_dot; p_dot; q_dot; r_dot];
        end
        %----------------------------

        % Compute the wind corrected body velocities
        function self=update_velocity_data(self, wind)

            % find wind in body velocity
            %TODO ACCESS GUST??
            wind_vel_body = Quaternion2Rotation(self.state(7:10))*wind + gust
            % R*steady wind + wind gust
            % wind velocity in body frame
            self.wind =  Quaternion2Rotation(self.state(7:10))*wind + gust
            % find relative wind

            % find airspeed in body frame. Va = [ur;vr;wr] ur = u-uw
            VaBody = [...
                        self.state(4) - self.wind(1);...
                        self.state(5) - self.wind(2);...
                        self.state(6) - self.wind(3)...
                    ];
            self.Va = sqrt(VaBody(1)^2+VaBody(2)^2+VaBody(3)^2);
            % find angle of attack
            self.alpha = tan^-1(VaBody(3)/VaBody(1));
            % find side slip
            self.beta = sin^-1(VaBody(2)/self.Va);
        end
        %----------------------------
        % compute for control surfaces
        function out=forces_moments(self, delta, MAV)
            [delta_e; delta_t; delta_a; delta_r] = delta;
            %% PROPELLER FORCE AND TORQUE Algorithim 4.1
            % Compute Vin , input voltage, using eq 4.22
                V_in = MAV.V_max*delta_t;
            % Compute Omega_prop, propeller speed, using the quadratic formula from eq 4.21
                % get quadratic coifficents
                
                a = (MAV.rho*MAV.D_prop^5/(2*pi)^2)*MAV.C_Q0;
                b = (MAV.rho*MAV.D_prop^4/(2*pi))*MAV.C_Q1*self.Va+(MAV.K_Q*MAV.K_V/MAV.R_motor);
                c = MAV.rho*MAV.D_prop^3*MAV.C_Q2*self.Va^2-(MAV.K_Q*V_in/MAV.R_motor)+MAV.K_Q*MAV.i0;
                % solve eq 4.21
                r = roots([a, b, c]);
                Omega_prop = r(r > 0);
            % Compute T_prop, thrust produced by propeller, using eq 4.17
                % find J(Omega, Va)
                J = 2*pi*self.Va/(Omega_prop*MAV.D_prop);
                % find CT(J)
                C_T_J = MAV.C_T0+MAV.C_T1*J+C_T2*J^2;
                T_prop = MAV.rho*MAV.D_prop^4/(4*pi^2)*Omega_prop^2*(MAV.C_T2*J^2+MAV.C_T1*J+MAV.C_T0);
            % Compute Q_prop, torque produced by propeller, using eq 4.18
                % Find C_Q(J)
                C_Q_J = MAV.C_Q0+MAV.C_Q1*J+C_Q2*J^2;
                % Find Q_Prop from eq 4.18
                Q_prop = MAV.rho*MAV.D_prop^5/(4*pi^2)*Omega_prop^2*(MAV.C_Q2*J^2+MAV.C_Q1*J+MAV.C_Q0);
            
            %% AERODYNAMIC FORCES [fx; fy; fz] eq 4.24

                % find coifficents eq 4.25                
                C_X = -C_D*cos(alpha) + C_L*sin(alpha);
                C_X_q = -C_D_q*cos(alpha) + C_L_q*sin(alpha);
                C_X_delta_e = -C_D_delta_e*cos(alpha) + C_L_delta_e*sin(alpha);
                C_Z = -C_D*sin(alpha) + C_L*cos(alpha);
                C_Z_q = -C_D_q*sin(alpha) + C_L_q*cos(alpha);
                C_Z_delta_e = -C_D_delta_e*sin(alpha) + C_L_delta_e*cos(alpha);
                % find force due to gravity
                F_g = [...
                        -m*g*sin(theta);
                        m*g*cos(theta)*sin(phi);
                        m*g*cos(theta)*cos(phi)...
                    ];
                % find aerodynamic force due to Lift and pitch rate
                F_aero = [...
                            C_X+C_X_q*MAV.c/(2*self.Va)*self.state(12);
                            MAV.C_Y_0+MAV.C_Y_beta*self.beta+MAV.C_Y_p*MAV.b/(2*self.Va)*state(11)+MAV.C_Y_r*MAV.b/(2*self.Va)*state(13);
                            C_Z+C_Z_q*MAV.c/(2*self.Va)*state(12)...
                        ];

                % find force due to control surfaces
                F_cs = [...
                        C_X_delta_e*delta_e;
                        C_Y_delta_a*delta_a+C_Y_delta_r*delta_r;
                        C_z_delta_e*delta_e...
                    ];
                % sum forces together as in eq 4.24
                Force = F_g+[T_prop;0;0]+1/2*MAV.rho*self.Va^2*MAV.S_wing*F_aero+1/2*MAV.rho*self.Va^2*MAV.S_wing*F_cs;
            %% AERODYNAMIC MOMENTS
                % Get first matrix of moment
                M_aero = [
                            MAV.b*(MAV.C_ell_0+C_ell_beta*self.beta+C_ell_p*(MAV.b/(2*self.Va)*p)+C_ell_r*(MAV.b/(2*self.Va)*r));
                            MAV.c*(MAV.C_m_0+C_m_alpha*self.alpha+C_m_q*(MAV.c/(2*self.Va)*q));
                            MAV.b*(MAV.C_n_0+C_n_beta*self.beta+C_n_p*(MAV.b/(2*self.Va)*p)+C_n_r*(MAV.b/(2*self.Va)*r))
                        ];
                
                % get moments from control surfaces
                M_cs = [...
                        MAV.b*(MAV.C_ell_delta_a*delta_a+MAV.C_ell_delta_r*delta_r);
                        MAV.c*(MAV.C_m_delta_e*delta_e);
                        MAV.b*(MAV.C_n_delta_a*delta_a+MAV.C_n_delta_r*delta_r)
                    ];
                % sum forces together as in eq 4.24
                Torque = M_g+1/2*MAV.rho*self.Va^2*MAV.S_wing*M_aero+1/2*MAV.rho*self.Va^2*MAV.S_wing*M_cs+[-Q_prop;0;0];

            % output total force and torque
            out = [Force'; Torque'];
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
            self.true_state.Va = self.Va;
            self.true_state.alpha = self.alpha;
            self.true_state.beta = self.beta;
            % ground speed
            self.true_state.Vg = 
            % find course angle
            self.true_state.chi = 
            % flight path angle
            self.true_state.gamma = 
            self.true_state.wn = self.wind(1);
            self.true_state.we = self.wind(2);
        end
    end
end