%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mavSimMatlab 
%     - Chapter 5 assignment for Beard & McLain, PUP, 2012
%     - Update history:  
%         2/5/2019 - RWB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run('../parameters/simulation_parameters')  % load SIM: simulation parameters
run('../parameters/aerosonde_parameters')  % load MAV: aircraft parameters
run('../parameters/wind_parameters')  % load WIND

% initialize the mav viewer
addpath('../chap2'); mav_view = mav_viewer();  
addpath('../chap3'); data_view = data_viewer();

% initialize the video writer
VIDEO = 0;  % 1 means write video, 0 means don't write video
if VIDEO==1
    video=video_writer('chap5_video.avi', SIM.ts_video);
end

% initialize elements of the architecture
addpath('../chap4'); 
addpath('../chap4'); wind = wind_simulation(SIM.ts_simulation, WIND);
mav = mav_dynamics(SIM.ts_simulation, MAV, WIND);

% compute trim
addpath('../chap5');
Va = MAV.Va0;
gamma = 0*pi/180;
[trim_state, trim_input] = compute_trim(mav, Va, gamma, MAV);
mav.state = trim_state;
delta = trim_input;

mav.update_velocity_data(zeros(6,1));   % no wind
fm   = mav.forces_moments(delta, MAV);
xdot = mav.derivatives(trim_state, fm, MAV)

xdot_target = [0; 0; -Va*sin(0); 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
err = xdot - xdot_target
temp1=Quaternion2Euler([trim_state(7);trim_state(8);trim_state(9);trim_state(10)]);

fprintf("computed trim values\n")
fprintf("TRIM STATE - Straight and level should have no:\n Side speed(v), roll angle (phi), yaw angle(psi), angular momentum (p,q,r)\n\
        \tpn = %.2f, pe - %.2f, pd = %.2f\n\
        \tu = %.4f, v = %.4f, w = %.4f\n\
        \tphi = %.4f, theta = %.4f, psi = %.4f\n\
        \tp = %.4f, q = %.4f, r = %.4f\n",... 
        trim_state(1),...
        trim_state(2),...
        trim_state(3),...
        trim_state(4),...
        trim_state(5),...
        trim_state(6),...
        temp1(1),...
        temp1(2),...
        temp1(3),...
        trim_state(11),...
        trim_state(12),...
        trim_state(13));
fprintf("TRIM INPUT- straight and level should onlt have elevator and throttle:\n\
        \tElevator: %.4f\n\
        \tThrottle: %.4f\n\
        \tAileron: %.4f\n\
        \tRudder: %.4f\n",...
        delta(1), delta(2), delta(3), delta(4));

gamma1 = atan(-trim_state(6)/trim_state(4)) ;
gamma_true = mav.true_state.gamma;
fprintf("Gamma calculated from trim - steady flight should be 0\n\tGamma: %.4f\n", gamma_true);
% compute linearized models
[A_lon, B_lon, A_lat, B_lat] = compute_ss_model(mav, trim_state, trim_input, MAV);

% %compute eigen values
% lambda = eig(A_lon)

% % find complex conjugate pair
% idx = find(imag(lambda) ~= 0);
% lambda_complex = lambda(idx);

% % extract real and imaginary parts
% sigma_1 = real(lambda_complex(1));
% omega_1 = abs(imag(lambda_complex(1)));

% omega_n_1 = sqrt(sigma_1^2 + omega_1^2)
% zeta_1 = -sigma_1 / omega_n_1

% % extract real and imaginary parts
% sigma_2 = real(lambda_complex(2));
% omega_2 = abs(imag(lambda_complex(2)));

% omega_n_2 = sqrt(sigma_2^2 + omega_2^2)
% zeta_1 = -sigma_2 / omega_n_2
% initialize the simulation time
sim_time = SIM.start_time;

lambda = eig(A_lat);
    d_trim = delta(1)
% main simulation loop
disp('Type CTRL-C to exit');
while sim_time < SIM.end_time

    %-------physical system-------------
    % current_wind = wind.update();
    current_wind = zeros(6,1);
    % fprintf("Chap 5 delta\n")
    % delta
    if abs(sim_time - 3.0) < 1e-6
        fprintf("phougoid impulese")
        delta(1) = d_trim+0.4;
    else
    delta(1) = d_trim;
    end
    % elseif abs(sim_time - (0.6+SIM.ts_simulation)) < 1e-6
    % delta(3) = -0.2;
    % end
    mav.update_state(delta, current_wind, MAV);
    mav.true_state;
    %-------update viewer-------------
    mav_view.update(mav.true_state);  % plot body of MAV
    data_view.update(mav.true_state,...  % true states
                     mav.true_state,...  % estimated states
                     mav.true_state,...  % commmanded states
                     SIM.ts_simulation); 
    if VIDEO==1
        video.update(sim_time);  
    end

    %-------increment time-------------
    sim_time = sim_time + SIM.ts_simulation;
end

if VIDEO==1
    video.close(); 
end

