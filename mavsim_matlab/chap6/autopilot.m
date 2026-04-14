% autopilot block for mavsim_matlab
%
% mavsim_matlab 
%     - Beard & McLain, PUP, 2012
%     - Last updated:  
%         2/13/2019 - RWB
classdef autopilot < handle
   %--------------------------------
    properties
        ts_control
        roll_from_aileron
        course_from_roll
        sideslip_from_rudder
        pitch_from_elevator
        altitude_from_pitch
        airspeed_from_pitch
        airspeed_from_throttle
        altitude_zone
        commanded_state
    end
    %--------------------------------
    methods
        %------constructor-----------
        function self = autopilot(ts_control)
             % load AP: control gains/parameters
            run('../parameters/control_parameters') 
            addpath('../chap6')
            self.roll_from_aileron = pd_control_with_rate(...
                                        AP.roll_kp,...
                                        AP.roll_kd,...
                                        45 * (pi/180));
            self.course_from_roll = pi_control(...
                                        AP.course_kp,...
                                        AP.course_ki,...
                                        ts_control,...
                                        30 * (pi/180));
            self.sideslip_from_rudder = pi_control(...
                                        AP.sideslip_kp,...
                                        AP.sideslip_ki,...
                                        ts_control,...
                                        45 * (pi/180));
            self.pitch_from_elevator = pd_control_with_rate(...
                                        AP.pitch_kp,...
                                        AP.pitch_kd,...
                                        45 * (pi/180));
            self.altitude_from_pitch = pi_control(...
                                        AP.altitude_kp,...
                                        AP.altitude_ki,...
                                        ts_control,...
                                        30 * (pi/180));
            self.airspeed_from_throttle = pi_control(...
                                        AP.airspeed_throttle_kp,...
                                        AP.airspeed_throttle_ki,...
                                        ts_control,...
                                        1);     
            self.altitude_zone = AP.altitude_zone;
            addpath('../message_types'); 
            self.commanded_state = msg_state();
        end
        %------methods-----------
        function [delta, commanded_state] = update(self, cmd, state)
            % update using the defined functions to get the commanded states as reswuested by the course command. 
            % this puts the close the loop steps into order
            % lateral  autopilot
            phi_c = self.roll_from_roll.update(cmd.course_command,state.chi, reset_flag=true);
            delta_a = self.roll_from_aileron.update(phi_c, state.phi,state.p);
            delta_r = self.sideslip_from_rudder.update(0,state.beta);

            % longitudinal autopilot
            h_c = cmd.altitude_command;
            theta_c = pi/16;
            delta_e = self.altitude_from_pitch.update(h_c,state.h);
            delta_t = self.pitch_from_elevator.update(theta_c, state.theta, state.q);
            delta_t = self.airspeed_from_throttle.update(cmd.airspeed_command, state.Va);

            % construct output and commanded states
            delta = [delta_e; delta_a; delta_r; delta_t];
            self.commanded_state.h = cmd.altitude_command;
            self.commanded_state.Va = cmd.airspeed_command;
            self.commanded_state.phi = phi_c;
            self.commanded_state.theta = theta_c;
            self.commanded_state.chi = cmd.course_command;
        	commanded_state = self.commanded_state;
        end
        %---------------------------
        function output = saturate(self, input, low_limit, up_limit)
            if input <= low_limit
                output = low_limit;
            elseif input >= up_limit
                output = up_limit;
            else
                output = input;
            end
        end
    end
end
% Implement the longitudinal autopilot and tune the gains. The
% input the longitudinal autopilot is the commanded airspeed
% and the commanded altitude. To make the process easier, you
% may want to initialize the system to trim inputs and then tune
% airspeed control first, followed by the pitch loop, and then the
% altitude loop.