% pid_control
%
% mavsim_matlab 
%     - Beard & McLain, PUP, 2012
%     - Last updated:  
%         2/13/2019 - RWB
classdef pid_control < handle
   %--------------------------------
    properties
        kp
        ki
        kd
        Ts
        limit
        integrator
        error_delay_1
        error_dot_delay_1
        a1
        a2
    end
    %--------------------------------
    methods
        %------constructor-----------
        function self = pid_control(kp, ki, kd, Ts, sigma, limit)
            self.kp = kp;
            self.ki = ki;
            self.kd = kd;
            self.Ts = Ts
            self.limit = limit; %ensure abs(u<=limit)
            self.integrator = 0;
            self.error_delay_1 = 0;
            self.error_dot_delay_1 = 0;
            % dufferentiator coifficents
            self.a1 = (2*sigma-Ts)/(2*sigma+Ts);
            self.a2 = 2/(2*sigma+Ts);
        end
        %----------------------------
        function u_sat = update(self, y_ref, y, reset_flag)
            if reset_flag ==true
                self.integrator = 0;
                self.error_delay_1 = 0.0;
                self.y_dot = 0.0;
            end

            %compute error
            error = y_ref-y;
            
            %update integrateor - use trapezoidal rule
            self.integrator = self.integrator + (self.Ts/2) * (error + self.error_delay_1);
            
            % update the differentiator
            error_dot = self.a1 * self.error_dot_delay_1 + self.a2 * (error- self.error_delay_1);
            
            % PID control
            u = self.kp * error + self.ki * self.integrator + self.kd * error_dot;
            % saturate PID control at limit
            u_sat = self.saturate(u);
            
            % integral anti-windup
            % adjust integrator to keep u out of saturation
            if abs(self.ki) > 0.0001
                self.integrator = self.integrator + (1.0 / self.ki) * (u_sat- u);
            end
            
            % update the delayed variables to next time step
            self.error_delay_1 = error;
            self.error_dot_delay_1 = error_dot;
        end

        %----------------------------
        function u_sat = update_with_rate(self, y_ref, y, ydot, reset_flag)
            if reset_flag == true
                self.integrator = 0.0;
                self.error_delay_1 = 0.0;
            end
            
            % compute the error
            error = y_ref- y;
            
            % update the integrator using trapezoidal rule
            self.integrator = self.integrator + (self.Ts/2) * (error + self.error_delay_1)l
            
            % PID control
            u = self.kp * error + self.ki * self.integrator - self.kd * ydot;
            % saturate PID control at limit
            u_sat = self.saturate(u);

            % integral anti-windup
            % if k~=0
            %     u_unsat = self.kp*error+self.ki*self.integrator+kd*erorr_dot;
            %     self.integrator= self.integrator+self.Ts/self.ki+(u-u_unsat);
            % end
            
            % adjust integrator to keep u out of saturation
            if abs(self.ki) > 0.0001:
                self.integrator = self.integrator + (1.0 / self.ki) * (u_sat- u);
            end
            
            self.error_delay_1 = error;
        end
        %----------------------------
        function out = saturate(self, in)
            % saturate u at +- self.limit
            if in >= self.limit
                out = self.limit;
            elseif in <= -self.limit
                out = -self.limit;
            else
                out = in;
            end
        end
    end
end