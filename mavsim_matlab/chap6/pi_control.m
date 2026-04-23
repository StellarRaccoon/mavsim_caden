% pi_control
%
% mavsim_matlab 
%     - Beard & McLain, PUP, 2012
%     - Last updated:  
%         2/13/2019 - RWB
classdef pi_control < handle
%--------------------------------
    properties
        kp
        ki
        Ts
        limit
        integrator
        error_delay_1
    end
    %--------------------------------
    methods
        %------constructor-----------
        function self = pi_control(kp, ki, Ts, limit)
            self.kp = kp;
            self.ki = ki;
            self.Ts = Ts;
            self.limit = limit;
            self.integrator = 0;
            self.error_delay_1 = 0;
        end
        %---------------------------- Similar ro PID, but remove differentiato edot and kd
        function u_sat = update(self, y_ref, y)
            if reset_flag ==true
                self.integrator = 0;
                self.error_delay_1 = 0.0;
            end
            %compute error
            error = y_ref-y;
            %update integrateor - use trapezoidal rule
            self.integrator = self.integrator + (self.Ts/2) * (error + self.error_delay_1);

            % PI control
            u = self.kp * error + self.ki * self.integrator;

            % saturate PID control at limit
            u_sat = self.saturate(u);
            % integral anti-windup
            % adjust integrator to keep u out of saturation
            if abs(self.ki) > 0.0001
                self.integrator = self.integrator + (1.0 / self.ki) * (u_sat- u);
            end
            % update the delayed variables to next time step
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