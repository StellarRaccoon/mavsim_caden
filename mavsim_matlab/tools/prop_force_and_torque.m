

function [T_prop, Q_prop] = prop_force_and_torque(Va, delta_t, MAV)
        % run('../parameters/aerosonde_parameters');  % load MAV: aircraft parameters
        V_in = MAV.V_max*delta_t;
        % Compute Omega_prop, propeller speed, using the quadratic formula from eq 4.21
        % get quadratic coifficents
        a = (MAV.rho*MAV.D_prop^5/((2*pi)^2))*MAV.C_Q0;
        b = (MAV.rho*MAV.D_prop^4/(2*pi))*MAV.C_Q1*Va+(MAV.K_Q*MAV.K_V/MAV.R_motor);
        c = MAV.rho*MAV.D_prop^3*MAV.C_Q2*Va^2-(MAV.K_Q*V_in/MAV.R_motor)+MAV.K_Q*MAV.i0;
        % I get nasty  root cannot be inf or NaN errors so prevent those
        if isnan(a) || isinf(a) || a == 0
            a = 0;
        end
        if isnan(b) || isinf(b)
            b = 0;
        end
        if isnan(c) || isinf(c)
            c = 0;
        end
        % solve eq 4.21
        r = roots([a, b, c]);

        % keep only real, positive roots
        r = r(imag(r)==0 & r > 0);
        if isempty(r)
            Omega_prop = 1e-3;
        else
            Omega_prop = max(r);
        end

    % Compute T_prop, thrust produced by propeller, using eq 4.17
        % find J(Omega, Va)
        J = 2*pi*Va/(Omega_prop*MAV.D_prop);
        J = min(max(J, 0), 5);   % clamp to reasonable range
        % find CT(J)
        C_T_J = MAV.C_T0+MAV.C_T1*J+MAV.C_T2*J^2;
        % T_prop = (MAV.rho*MAV.D_prop^4/(4*pi^2))*Omega_prop^2*(MAV.C_T2*J^2+MAV.C_T1*J+MAV.C_T0);
        n = Omega_prop / (2*pi);
        T_prop = MAV.rho * n^2 * MAV.D_prop^4 * C_T_J;
    % Compute Q_prop, torque produced by propeller, using eq 4.18
        % Find C_Q(J)
        C_Q_J = MAV.C_Q0+MAV.C_Q1*J+MAV.C_Q2*J^2;
        % Find Q_Prop from eq 4.18
        Q_prop = MAV.rho*MAV.D_prop^5/(4*pi^2)*Omega_prop^2*(MAV.C_Q2*J^2+MAV.C_Q1*J+MAV.C_Q0);
end

