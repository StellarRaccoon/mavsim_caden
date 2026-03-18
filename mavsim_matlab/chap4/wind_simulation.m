% wind Simulation - simulates steady state and wind gusts
%
% mavMatSim 
%     - Beard & McLain, PUP, 2012
%     - Update history:  
%         12/27/2018 - RWB
classdef wind_simulation < handle
   %--------------------------------
    
    properties
        steady_state
        A
        B
        C
        gust_state
        _gust
        Ts
    end
    %--------------------------------
    methods
        %------constructor-----------
        function self = wind_simulation(Ts)
            run('../parameters/wind_parameters')  % load WIND
            self.steady_state = [0;0;0];%todo
            %H(u)
            Ku = WIND.sigma_u * sqrt(2*WIND.Va0/(pi*WIND.L_u));
            Hu = Ku * tf(1, [WIND.L_u/WIND.Va0 1]);

            %H(v)
            Kv = WIND.sigma_v * sqrt(3*WIND.Va0/(pi*WIND.L_v));
            num = [2*sqrt(3)*WIND.L_v/WIND.Va0 1];
            den = conv([WIND.L_v/WIND.Va0 1], [WIND.L_v/WIND.Va0 1]);   % (1 + L/V s)^2
            Hv = Kv * tf(num, den);

            %H(w)
            Kw = WIND.sigma_w * sqrt(3*WIND.Va0/(pi*WIND.L_w));
            num = [2*sqrt(3)*WIND.L_w/WIND.Va0 1];
            den = conv([WIND.L_w/WIND.Va0 1], [WIND.L_w/WIND.Va0 1]);
            Hw = Kw * tf(num, den);

            H = blkdiag(Hu, Hv, Hw);
            stateSpace = ss(H);
            self.A = stateSpace.A;
            self.B = stateSpace.B;
            self.C = stateSpace.C;
            self.gust_state = zeros(size(self.A,1),1)
            self._gust = [0; 0; 0];
            self.Ts = Ts;
        end
        %---------------------------
        function wind=update(self)
            wind = [self.steady_state; self._gust];
        end
        %----------------------------
        function self = gust(self)
            w = randn(3,1);
            %state space equation
            % x = Ax +Bw
            self.gust_state = self.gust_state + self.Ts*(self.A*self.gust_state + self.B*w);
            % gust velocities = C*x
            self._gust = self.C*self.gust_state;
        end
    end
end