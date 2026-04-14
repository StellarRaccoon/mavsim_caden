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
        gust
        Ts
    end
    %--------------------------------
    methods
        %------constructor-----------
        function self = wind_simulation(Ts, WIND)
            self.steady_state = [WIND.wind_n;WIND.wind_e;WIND.wind_d];

            %apply Dryden filters  and transform functions
            %H(u)
            outer = WIND.sigma_u * sqrt(2*WIND.Va0/(pi*WIND.L_u));
            Hu = tf(outer, [WIND.L_u/WIND.Va0, 1]);

            %H(v)
            outer = WIND.sigma_v * sqrt(3*WIND.Va0/(pi*WIND.L_v));
            num = [outer, outer*(WIND.Va0/(sqrt(3)*WIND.L_v))];
            den = [1, 2*WIND.Va0/WIND.L_v, (WIND.Va0/WIND.L_v)^2];   % (1 + L/V s)^2
            Hv = tf(num, den);

            %H(w)
            outer = WIND.sigma_w * sqrt(3*WIND.Va0/(pi*WIND.L_w));
            num = [outer, outer*(WIND.Va0/(sqrt(3)*WIND.L_w))];
            den = [1, 2*WIND.Va0/WIND.L_w, (WIND.Va0/WIND.L_w)^2];   % (1 + L/V s)^2
            Hw = tf(num, den);

            % create diagnoal matrix
            H = blkdiag(Hu, Hv, Hw);
            stateSpace = ss(H);
            self.A = stateSpace.A;
            self.B = stateSpace.B;
            self.C = stateSpace.C;
            self.gust_state = zeros(size(self.A,1),1);
            self.gust = [0; 0; 0];
            self.Ts = Ts;
        end
        %---------------------------
        function wind=update(self)
            % self.gust();
            wind = [self.steady_state; self.gust];
        end
        %----------------------------
        function self = gust(self)
            w = randn(3,1);
            %state space equation
            % x = Ax +Bw
            self.gust_state = self.gust_state + self.Ts*(self.A*self.gust_state + self.B*w);
            % gust velocities = C*x
            self.gust = self.C*self.gust_state;
        end
    end
end