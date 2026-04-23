% control_parameters
%
% mavsim_matlab 
%     - Beard & McLain, PUP, 2012
%     - Last updated:  
%         2/13/2019 - RWB
addpath('../chap5')
load transfer_function_coef.mat

% AP stands for autopilot
AP.gravity = MAV.gravity;
AP.sigma = 0.05;
AP.Va0 = MAV.Va0;
%TODO what is omega_n in the equations
%----------roll loop-------------
% these need to be determined from trial and error
wn_phi = 11;
zeta_phi = 0.707;
%eq 6.1 and 6.2
AP.roll_kp = omega_nphi^2/a_phi2;
AP.roll_kd = ((2*zeta_phi*omega_nphi)-a_phi1)/a_phi2;

%----------course loop-------------
zeta_chi = 0.707
omega_nchi = 0.5
%eq 6.5
AP.course_kp =  2*zeta_chi*omega_nchi*AP.Va0/AP.gravity;
AP.course_ki = (omega_nchi^2* AP.Va0)/AP.gravity;

% %----------sideslip loop-------------
% AP.sideslip_ki = 
% AP.sideslip_kp = 

%----------pitch loop-------------
omega_ntheta = 15;
zeta_theta = 0.8;
%eq 6.9-6.11
AP.pitch_kp = (omega_ntheta^2-a_theta2)/a_theta3;
AP.pitch_kd = (2*zeta_theta*omega_ntheta-a_theta1)/a_theta3;
K_theta_DC = AP.pitch_kp*a_theta3/omega_ntheta^2;

%----------altitude loop-------------
zeta_h = 0.8;
wn_h = 0.25;
AP.altitude_kp = 2*zeta_h*omega_nh/(K_theta_DC*AP.Va0); 
AP.altitude_ki = omega_nh^2/(K_theta_DC*AP.Va0); %eq 6.12
AP.altitude_zone = 10; 

%---------airspeed hold using throttle---------------
omega_nv = 0.5;
zeta_v = 0.8;
AP.airspeed_throttle_kp = (2*zeta_v*omega_nv-a_v1)/a_v2; %eq 6.15
AP.airspeed_throttle_ki = omega_nv^2/a_v2; %eq 6.15
