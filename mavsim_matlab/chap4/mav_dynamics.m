% mav dynamics - implement rigid body dynamics for mav
%
% mavMatSim 
%     - BeMAV.ard & McLain, PUP, 2012
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
        Va_body
        true_state
    end
    %--------------------------------
    methods
        %------constructor-----------
        % initiate with constant wind and perfectly aligned using inital values
        function self = mav_dynamics(Ts, MAV, WIND)
            self.ts_simulation = Ts; % time step between function calls
            self.state = [MAV.pn0; MAV.pe0; MAV.pd0; MAV.u0; MAV.v0; MAV.w0;
                MAV.e0; MAV.e1; MAV.e2; MAV.e3; MAV.p0; MAV.q0; MAV.r0];
            self.Va = MAV.Va0;
            % alpha is angle of attack. angle of up/down tilt on wings
            self.alpha = 0; 
            % beta is side slip angle. angle that plane is facing off straight path
            self.beta = 0;
            self.wind = [WIND.wind_n;WIND.wind_e;WIND.wind_d];
            addpath('../message_types'); self.true_state = msg_state();
        end
        %---------------------------
        function self=update_state(self, delta, wind, MAV)
            % determine current state of aircraft from the derivatives
            % find current location and orientation of craft based on the force and moments which are time derivatives
            % we solve the time derivatives for the current time
            %
            % Integrate the differential equations defining dynamics
            % forces_moments MAV.are the forces and moments on the MAV.
            % 

            % Determine the current force and moments (rotation and location derivatives) acting on the body
            forces_moments = self.forces_moments(delta, MAV);
            
            % Integrate ODE using Runge-Kutta RK4 algorithm
            % undo the derivative to get the current loc, orientation of the craft
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
            % find the true state by adjusting for wind
            self.update_true_state(wind);
        end
        %----------------------------
        function xdot = derivatives(self, state, forces_moments, MAV)
            %get the current location of aircraft
            pn    = state(1);
            pe    = state(2);
            pd    = state(3);
            u     = state(4);
            v     = state(5);
            w     = state(6);
            e0    = state(7);
            e1    = state(8);
            e2    = state(9);
            e3    = state(10);
            %angular rates in body frame
            p     = state(11);
            q     = state(12);
            r     = state(13);
            %get the current forces and moments of the aircraft (inital values)
            fx    = forces_moments(1);
            fy    = forces_moments(2);
            fz    = forces_moments(3);
            ell   = forces_moments(4);
            m     = forces_moments(5);
            n     = forces_moments(6);
            

            % position kinematics - Velocity of aircraft in north, east, down direction wrt to ground
            % found by rotating the ground velocity measured from the aircraft to the ground frame (straigting it out)
            % dpos/dt = velocity rotated to align w/ ground
            pos_dot = Quaternion2Rotation([e0, e1, e2, e3]) * [u; v; w];
            pn_dot = pos_dot(1);
            pe_dot = pos_dot(2);
            pd_dot = pos_dot(3);
            %fprintf("pn_dot: %d pe_dot: %d pd_dot: %d\n", pn_dot, pe_dot, pd_dot);

            % position dynamics - Acceleration a=f/m
            % Forces acting on vehicle: gravity and external
            % project the angular velocities to the ground velocity in the body frame (cross [pq,r] x [uvw])
            % 1/m * f is acceleration in from applied forces (fx,fy,fz are found later)
            vel_dot = [r*v-q*w; p*w-r*u; q*u-p*v;] + MAV.mass^(-1)*[fx; fy; fz];
            u_dot = vel_dot(1); 
            v_dot = vel_dot(2);
            w_dot = vel_dot(3);
            %fprintf("u_dot: %d v_dot: %d w_dot: %d\n", u_dot, v_dot, w_dot);

            % rotational kinematics
            % velocity of spin - rotate the angular rates to ground frame
            temp = [
                    0, -p, -q, -r;
                    p, 0, r, -q;
                    q, -r, 0, p;
                    r, q, -p, 0
                    ];
            % rotation matrix * body rate
            e_dot = 0.5 * temp * [e0; e1; e2; e3];
            e0_dot = e_dot(1);
            e1_dot = e_dot(2);
            e2_dot = e_dot(3);
            e3_dot = e_dot(4);

            % rotational dynamics - 
            % angular acceleration momentum
            % uses the tendency to appose acceleration of rotation
            rate_dot= [
                MAV.Gamma1*p*q-MAV.Gamma2*q*r; 
                MAV.Gamma5*p*r-MAV.Gamma6*(p^2-r^2); 
                MAV.Gamma7*p*q-MAV.Gamma1*q*r]+...
                [MAV.Gamma3, 0, MAV.Gamma4; 
                0, 1/MAV.Jy, 0; 
                MAV.Gamma4, 0, MAV.Gamma8;]* [ell; m; n];

            p_dot = rate_dot(1);
            q_dot = rate_dot(2);
            r_dot = rate_dot(3);
            % collect all the derivaties of the states
            xdot = [pn_dot; pe_dot; pd_dot; u_dot; v_dot; w_dot;...
                    e0_dot; e1_dot; e2_dot; e3_dot; p_dot; q_dot; r_dot];
        end
        %----------------------------

        % Compute the wind corrected body velocities
        % update the mav data to include effects of wind 
        function self=update_velocity_data(self, wind)

            % find wind in body velocity
            steady_wind = wind(1:3);
            gust_wind = wind(4:6);
            
            % R*steady wind + wind gust
            % wind velocity in body frame
            V_wind =  Quaternion2Rotation(self.state(7:10))'*steady_wind +gust_wind;
            % find relative wind
            % find airspeed in body frame. Va = [ur;vr;wr] ur = u-uw
                Va_body= [...
                        self.state(4) - V_wind(1);...
                        self.state(5) - V_wind(2);...
                        self.state(6) - V_wind(3)...
                    ];
            self.Va_body=Va_body;
            
            %update state
            self.wind = V_wind;
            self.Va = sqrt(Va_body(1)^2+Va_body(2)^2+Va_body(3)^2);
            % find angle of attack
            self.alpha = atan2(Va_body(3),Va_body(1));
            % find side slip
            self.beta = asin(Va_body(2)/self.Va);
        end
        %----------------------------
        % Determine the forces and moments effecting the aircraft to use when caluclating derivatives
        % forces and momentum rely on the control surfaces, velocity of aircraft, and properties of aircraft
        function out=forces_moments(self, delta, MAV)
        
            e0 = self.state(7);
            e1 = self.state(8);
            e2 = self.state(9);
            e3 = self.state(10);
            p = self.state(11);
            q = self.state(12);
            r = self.state(13);

            delta_e = delta(1);
            delta_t = delta(2);
            delta_a = delta(3);
            delta_r = delta(4);
            %% PROPELLER FORCE AND TORQUE Algorithim 4.1
            [~, Q_prop] = prop_force_and_torque(self.Va, delta_t,MAV);
            % This uses a simpliefied propeller force 
            F_prop = 0.5*MAV.rho*MAV.S_prop*MAV.C_prop*((MAV.k_motor*delta_t)^2-self.Va^2);
            
            %% AERODYNAMIC FORCES [fx; fy; fz] eq 4.24
            % effect of lift and drag on aircraft
                [C_D, C_L] = drag_and_lift_coiff(self.alpha, MAV);
                
                % find coifficents eq 4.25                
                C_X = -C_D*cos(self.alpha) + C_L*sin(self.alpha);
                C_X_q = -MAV.C_D_q*cos(self.alpha) + MAV.C_L_q*sin(self.alpha);
                C_X_delta_e = -MAV.C_D_delta_e*cos(self.alpha) + MAV.C_L_delta_e*sin(self.alpha);
                C_Z = -C_D*sin(self.alpha) - C_L*cos(self.alpha);
                C_Z_q = -MAV.C_D_q*sin(self.alpha) - MAV.C_L_q*cos(self.alpha);
                C_Z_delta_e = -MAV.C_D_delta_e*sin(self.alpha) - MAV.C_L_delta_e*cos(self.alpha);
                
            %% GRAVITATIONAL FORCE
            % F= mg, reminder gravity points striaght down
            F_g = MAV.mass*MAV.gravity*[
                                    2*(e1*e3-e2*e0);
                                    2*(e2*e3 + e1*e0);
                                    e3^2+e0^2-e1^2-e2^2
                                    ];
            e = Quaternion2Euler(self.state(7:10));
            phi=e(1); theta=e(2); psi=e(3);
                % find aerodynamic force due to Lift and pitch rate
            c = MAV.c/(2*self.Va);
            b = MAV.b/(2*self.Va);

            Fa = 0.5*MAV.rho*self.Va^2*MAV.S_wing * [
                C_X + C_X_q*c*q + C_X_delta_e*delta_e;
                MAV.C_Y_0 + MAV.C_Y_beta*self.beta + MAV.C_Y_p*b*p + MAV.C_Y_r*b*r + MAV.C_Y_delta_a*delta_a + MAV.C_Y_delta_r*delta_r;
                C_Z + C_Z_q*c*q + C_Z_delta_e*delta_e
            ];
                
                % sum forces together as in eq 4.24
                Force = F_g+Fa+[F_prop;0;0];
            %% AERODYNAMIC MOMENTS
                    p = self.state(11); % p
                    q = self.state(12); % q
                    r = self.state(13); % r
                % Get first matrix of moment
                M_aero = [
                            MAV.b*(MAV.C_ell_0+MAV.C_ell_beta*self.beta+MAV.C_ell_p*(MAV.b/(2*self.Va)*p)+MAV.C_ell_r*(MAV.b/(2*self.Va)*r));
                            MAV.c*(MAV.C_m_0+MAV.C_m_alpha*self.alpha+MAV.C_m_q*(MAV.c/(2*self.Va)*q));
                            MAV.b*(MAV.C_n_0+MAV.C_n_beta*self.beta+MAV.C_n_p*(MAV.b/(2*self.Va)*p)+MAV.C_n_r*(MAV.b/(2*self.Va)*r))
                        ];
                % get moments from control surfaces
                M_cs = [...
                        MAV.b*(MAV.C_ell_delta_a*delta_a+MAV.C_ell_delta_r*delta_r);
                        MAV.c*(MAV.C_m_delta_e*delta_e);
                        MAV.b*(MAV.C_n_delta_a*delta_a+MAV.C_n_delta_r*delta_r)
                    ];
                    
                    
                % sum forces together as in eq 4.24
                Torque = 1/2*MAV.rho*self.Va^2*MAV.S_wing*M_aero+1/2*MAV.rho*self.Va^2*MAV.S_wing*M_cs+[-Q_prop;0;0];
            % output total force and torque
            out = [Force; Torque];
        end
        %----------------------------
        function self=update_true_state(self, wind)
            e = Quaternion2Euler(self.state(7:10));
            self.true_state.pn = self.state(1);  % pn
            self.true_state.pe = self.state(2);  % pd
            self.true_state.h = -self.state(3);  % h
            self.true_state.phi = e(1); % phi
            self.true_state.theta = e(2); % theta
            self.true_state.psi = e(3); % psi
            self.true_state.p = self.state(11); % p
            self.true_state.q = self.state(12); % q
            self.true_state.r = self.state(13); % r
            self.true_state.Va = self.Va;
            self.true_state.alpha = self.alpha;
            self.true_state.beta = self.beta;
            % ground speed
            % inertial wind from wind_simulation
            wn = wind(1);
            we = wind(2);
            wd = wind(3);
            
            % Vg = self.state(4:6);
            % Vg_M = norm(Vg);
            % Vg_h = [self.state(4),self.state(5),0];
            % Vg_h_M = norm(Vg_h);
            % %horizontal airspeed
            % Va_h = [self.Va_body(1),self.Va_body(2),0];
            % Va_h_M = norm(Va_h);
            % %determine gamma
            % self.true_state.gamma = acos((Vg'*Vg_h')/(Vg_M*Vg_h_M));
            % fprintf("Gamma %.5f\n",self.true_state.gamma);
            % %determine course angle chi
            % num = Vg_h*Va_h';
            % den = (Vg_h_M*Va_h_M);
            % frac = num/den; 
            % frac = max(min(frac,1),-1);% clamp prevent going outside domain of -1 to 1
            % self.true_state.chi = self.true_state.psi + self.beta + acos(frac);
            % fprintf("Chi %.5f\n",self.true_state.chi)
            % % transform velocity back to interial fram to get ground velocity
            % R_vb = Quaternion2Rotation(self.state(7:10));
            % V_body = [self.state(4); self.state(5); self.state(6)];
            % V_inertial = R_vb * V_body;
R_vb = Quaternion2Rotation(self.state(7:10));  % body → inertial
V_body = [self.state(4); self.state(5); self.state(6)];
V_inertial = R_vb * V_body;                    % body → inertial

            % ground velocity
            Vg_vec = V_inertial + [wn; we; wd];

            % Update remaining states
            self.true_state.Vg = norm(Vg_vec);
            self.true_state.chi = atan2(Vg_vec(2), Vg_vec(1));
            %current flight path angle - angle between velocity and horizonntal plane (pitch angle to aoa)
            % self.true_state.gamma = self.true_state.theta - self.alpha;%-asin(Vg_vec(3)/self.true_state.Vg);
            self.true_state.gamma = atan2(-Vg_vec(3), sqrt(Vg_vec(1)^2 + Vg_vec(2)^2));
% self.msg_true_state.gamma = accos(Vg.dot(Vg_h)/(Vg_M*Vg_h_M))
            self.true_state.wn = wind(1);
            self.true_state.we = wind(2);

        end
    end
end