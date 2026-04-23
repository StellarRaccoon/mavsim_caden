% pd_control_with_rate
%
% mavsim_matlab 
%     - Beard & McLain, PUP, 2012
%     - Last updated:  
%         2/13/2019 - RWB
classdef pd_control_with_rate < handle
   %--------------------------------
    properties
        kp
        kd
        limit
    end
    %--------------------------------
    methods
        %------constructor-----------
        function self = pd_control_with_rate(kp, kd, limit)
            self.kp = kp;
            self.kd = kd;
            self.limit = limit;
        end
        %----------------------------
        function u_sat = update(self, y_ref, y, ydot)
            % compute the error
            error = y_ref- y;
            % PD control
            u = self.kp * error  - self.kd * ydot;
            % saturate PD control at limit
            u_sat = self.saturate(u);
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