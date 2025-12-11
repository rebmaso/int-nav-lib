#include "navigation.h"

namespace intnavlib {

template<typename T>
template<int n_x>
void Navigation<T>::predictKF(const Eigen::Matrix<T, n_x, n_x> & Phi_matrix, 
                                const Eigen::Matrix<T, n_x, n_x> & Q_matrix,
                                const Eigen::Matrix<T, n_x, n_x> & P_matrix_old, 
                                Eigen::Matrix<T, n_x, n_x> & P_matrix) {

    P_matrix = Phi_matrix * (P_matrix_old + 0.5 * Q_matrix) * Phi_matrix.transpose() + 0.5 * Q_matrix;

}

template<typename T>
template<int n_x, int n_z, int max_n_z>
bool Navigation<T>::updateKF(const Eigen::Matrix<T, n_z, 1, 0, max_n_z, 1> & delta_z,
            const Eigen::Matrix<T, n_x, n_x> & P_matrix, 
            const Eigen::Matrix<T, n_z, n_x, 0, max_n_z, n_x> & H_matrix,
            const Eigen::Matrix<T, n_z, n_z> & R_matrix,
            const T & p_value,
            Eigen::Matrix<T, n_z, n_z, 0, max_n_z, max_n_z> & S_matrix,
            Eigen::Matrix<T, n_x, 1> & x_est_new,
            Eigen::Matrix<T, n_x, n_x> & P_matrix_post) {

    // Calculate Kalman gain using (3.21)
    S_matrix = H_matrix * P_matrix * H_matrix.transpose() + R_matrix;
    
    // auto S_matrix_inv = S_matrix.inverse();

    Eigen::Matrix<T, n_z, n_z, 0, max_n_z, max_n_z> S_matrix_inv = S_matrix.template selfadjointView<Eigen::Upper>().llt().solve(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(R_matrix.rows(), R_matrix.cols()));

    Eigen::Matrix<T, n_x, n_z, 0, n_x, max_n_z> K_matrix = P_matrix * H_matrix.transpose() * S_matrix_inv;

    // Update error state estimates using (3.24)
    // A priori error state is always zero in closed loop filter
    x_est_new = /*x_est_propagated + */ K_matrix * delta_z;

    // Update state estimation error covariance matrix using (3.25)
    P_matrix_post = (Eigen::Matrix<T, n_x, n_x>::Identity() - K_matrix * H_matrix) * P_matrix;

    // Real-time consistency check
    // See: Estimation with Applications to Tracking and Navigation -- Yaakov Bar-Shalom et al. p.237
    // The quadratic form res' * inv(S) * res is a determination of the standard Chi squared distribution with one dof.
    // We check it's inside the desired p_value
    // Simplified version: no history, just current measurement, 1 * n_res dof's.
    T chi2 = delta_z.dot(S_matrix_inv * delta_z);
    unsigned int n_dofs = delta_z.rows();
    boost::math::chi_squared chi2_dist(n_dofs);
    T chi2_thresh = boost::math::quantile(chi2_dist, p_value);
    if (chi2 > chi2_thresh)
        return false;
    else return true;

}

template<typename T>
typename Types<T>::NavSolutionEcef 
Navigation<T>::navEquationsEcef(const NavSolutionEcef& old_nav, 
                                const ImuMeasurements& imu_meas, 
                                const T& tor_i) {
    NavSolutionEcef new_nav;
    new_nav.time = imu_meas.time;
    // ATTITUDE UPDATE
    // From (2.145) determine the Earth rotation over the update interval
    // C_Earth = C_e_i' * old_C_e_i
    T alpha_ie = kOmega_ie * tor_i;
    T cos_alpha_ie = cos(alpha_ie);
    T sin_alpha_ie = sin(alpha_ie);
    Matrix3 C_Earth;
    C_Earth << cos_alpha_ie, sin_alpha_ie, 0.0,
                -sin_alpha_ie, cos_alpha_ie, 0.0,
                            0.0,             0.0,  1.0;
    // Calculate attitude increment, magnitude, and skew-symmetric matrix
    Vector3 alpha_ib_b = imu_meas.omega * tor_i;
    T mag_alpha = alpha_ib_b.norm();
    Matrix3 Alpha_ib_b = skewSymmetric(alpha_ib_b);  
    // Obtain coordinate transformation matrix from the new attitude w.r.t. an
    // inertial frame to the old using Rodrigues' formula, (5.73)
    Matrix3 C_new_old;
    if (mag_alpha>1.0e-8) {
        C_new_old = Matrix3::Identity() + 
                    ((sin(mag_alpha) / mag_alpha) * Alpha_ib_b) +
                    ((1.0 - cos(mag_alpha)) / pow(mag_alpha,2.0)) * Alpha_ib_b * Alpha_ib_b;
    }
    else {
        C_new_old = Matrix3::Identity() + Alpha_ib_b;
    }
    // Update attitude using (5.75)
    new_nav.C_b_e = C_Earth * old_nav.C_b_e * C_new_old;
    // SPECIFIC FORCE FRAME TRANSFORMATION
    // Calculate the average body-to-ECEF-frame coordinate transformation
    // matrix over the update interval using (5.84) and (5.85)
    Vector3 alpha_ie_vec;
    alpha_ie_vec << 0.0 , 0.0 , alpha_ie;   
    Matrix3 ave_C_b_e;
    if (mag_alpha>1.0e-8) {
        ave_C_b_e = old_nav.C_b_e * 
            (Matrix3::Identity() + 
            ((1.0 - cos(mag_alpha)) / pow(mag_alpha,2.0)) *
            Alpha_ib_b + 
            ((1.0 - sin(mag_alpha) / mag_alpha) / pow(mag_alpha,2.0)) * 
            Alpha_ib_b * Alpha_ib_b) - 
            0.5 * skewSymmetric(alpha_ie_vec) * old_nav.C_b_e;
    }
    else { // Approximate if angle small enough (sum not multiply)
        ave_C_b_e = old_nav.C_b_e -
            0.5 * skewSymmetric(alpha_ie_vec) * old_nav.C_b_e;
    }
    // Transform specific force to ECEF-frame resolving axes using (5.85)
    Vector3 f_ib_e = ave_C_b_e * imu_meas.f;
    // UPDATE VELOCITY
    // From (5.36)
    Vector3 kOmega_ie_vec;
    kOmega_ie_vec << 0.0 , 0.0 , kOmega_ie;   
    new_nav.v_eb_e = old_nav.v_eb_e + tor_i * (f_ib_e + gravityEcef(old_nav.r_eb_e) -
                            2.0 * skewSymmetric(kOmega_ie_vec) * old_nav.v_eb_e);
    // UPDATE CARTESIAN POSITION
    // From (5.38),
    new_nav.r_eb_e = old_nav.r_eb_e + (new_nav.v_eb_e + old_nav.v_eb_e) * 0.5 * tor_i; 
    return new_nav;
}

template<typename T>
typename Types<T>::StateEstEcef 
Navigation<T>::lcPredictKF(const StateEstEcef & state_est_old, 
                            const ImuMeasurements & imu_meas,
                            const KfConfig & lc_kf_config,
                            const T & tor_i) {
    // Compensate IMU measurements
    ImuMeasurements imu_meas_comp = imu_meas;
    imu_meas_comp.f -= state_est_old.acc_bias;
    imu_meas_comp.omega -= state_est_old.gyro_bias;
    // Prop uncertainty
    StateEstEcef state_est_ecef;
    state_est_ecef.valid = state_est_old.valid;
    state_est_ecef.acc_bias = state_est_old.acc_bias;
    state_est_ecef.gyro_bias = state_est_old.gyro_bias;
    state_est_ecef.innovations_sigmas = state_est_old.innovations_sigmas;
    state_est_ecef.P_matrix.template block<15,15>(0,0) = lcPropUnc(state_est_old.P_matrix.template block<15,15>(0,0), 
                                                            state_est_old.nav_sol,
                                                            ecefToNed(state_est_old.nav_sol),
                                                            imu_meas_comp,
                                                            lc_kf_config,
                                                            tor_i);
    // Predict state
    state_est_ecef.nav_sol = navEquationsEcef(state_est_old.nav_sol, imu_meas_comp, tor_i);
    return state_est_ecef;
}

template<typename T>
typename Types<T>::StateEstEcef 
Navigation<T>::tcPredictKF(const StateEstEcef & state_est_old, 
                            const ImuMeasurements & imu_meas,
                            const KfConfig & kf_config,
                            const T & tor_i) {
    // Compensate IMU measurements
    ImuMeasurements imu_meas_comp = imu_meas;
    imu_meas_comp.f -= state_est_old.acc_bias;
    imu_meas_comp.omega -= state_est_old.gyro_bias;
    // Prop uncertainty
    StateEstEcef state_est_ecef;
    state_est_ecef.valid = state_est_old.valid;
    state_est_ecef.acc_bias = state_est_old.acc_bias;
    state_est_ecef.gyro_bias = state_est_old.gyro_bias;
    state_est_ecef.clock_offset = state_est_old.clock_offset;
    state_est_ecef.clock_drift = state_est_old.clock_drift;
    state_est_ecef.innovations_sigmas = state_est_old.innovations_sigmas;
    state_est_ecef.P_matrix = tcPropUnc(state_est_old.P_matrix, 
                                        state_est_old.nav_sol,
                                        ecefToNed(state_est_old.nav_sol),
                                        imu_meas_comp,
                                        kf_config,
                                        tor_i);
    // Predict state
    state_est_ecef.nav_sol = navEquationsEcef(state_est_old.nav_sol, imu_meas_comp, tor_i);
    return state_est_ecef;
}

template<typename T>
Eigen::Matrix<T,15,15> 
Navigation<T>::lcPropUnc(const Eigen::Matrix<T,15,15> & P_matrix_old, 
                                const NavSolutionEcef & old_nav_est_ecef,
                                const NavSolutionNed & old_nav_est_ned,
                                const ImuMeasurements & imu_meas,
                                const KfConfig & lc_kf_config,
                                const T & tor_i) {

    T geocentric_radius = kR0 / sqrt(1.0 - pow(kEccentricity * sin(old_nav_est_ned.latitude),2.0)) *
        sqrt(pow(cos(old_nav_est_ned.latitude), 2.0) + pow(1.0 - kEccentricity*kEccentricity, 2.0) * pow(sin(old_nav_est_ned.latitude), 2.0)); // from (2.137)

    // Skew symmetric matrix of Earth rate
    Vector3 kOmega_ie_vec;
    kOmega_ie_vec << 0.0,0.0,kOmega_ie;
    Matrix3 Omega_ie = skewSymmetric(kOmega_ie_vec);

    // Determine error-state transition matrix using (14.50) (first-order approx)
    Eigen::Matrix<T,15,15> Phi_matrix = Eigen::Matrix<T,15,15>::Identity();
    Phi_matrix.template block<3,3>(0,0) -= Omega_ie * tor_i;
    Phi_matrix.template block<3,3>(0,12) = old_nav_est_ecef.C_b_e * tor_i;
    Phi_matrix.template block<3,3>(3,0) = - tor_i * skewSymmetric(Vector3(old_nav_est_ecef.C_b_e * imu_meas.f));
    Phi_matrix.template block<3,3>(3,3) -= 2.0 * Omega_ie * tor_i;
    Phi_matrix.template block<3,3>(3,6) = -tor_i * 2 * gravityEcef(old_nav_est_ecef.r_eb_e) /
        geocentric_radius * old_nav_est_ecef.r_eb_e.transpose() / old_nav_est_ecef.r_eb_e.norm();
    Phi_matrix.template block<3,3>(3,9) = old_nav_est_ecef.C_b_e * tor_i;
    Phi_matrix.template block<3,3>(6,3) = Matrix3::Identity() * tor_i;

    // Determine approximate system noise covariance matrix using (14.82)
    Eigen::Matrix<T,15,15> Q_prime_matrix = Eigen::Matrix<T,15,15>::Zero();
    Q_prime_matrix.template block<3,3>(0,0) = Matrix3::Identity() * lc_kf_config.gyro_noise_psd * tor_i;
    Q_prime_matrix.template block<3,3>(3,3) = Matrix3::Identity() * lc_kf_config.accel_noise_psd * tor_i;
    Q_prime_matrix.template block<3,3>(9,9) = Matrix3::Identity() * lc_kf_config.accel_bias_psd * tor_i;
    Q_prime_matrix.template block<3,3>(12,12) = Matrix3::Identity() * lc_kf_config.gyro_bias_psd * tor_i;

    // Propagate state estimation error covariance matrix using (3.46)
    Eigen::Matrix<T,15,15> P_matrix;
    predictKF<15>(Phi_matrix, Q_prime_matrix, P_matrix_old, P_matrix);
    
    return P_matrix;
}

template<typename T>
Eigen::Matrix<T,17,17> 
Navigation<T>::tcPropUnc(const Eigen::Matrix<T,17,17> & P_matrix_old, 
                                const NavSolutionEcef & old_nav_est_ecef,
                                const NavSolutionNed & old_nav_est_ned,
                                const ImuMeasurements & imu_meas,
                                const KfConfig & tc_kf_config,
                                const T & tor_i) {

    T geocentric_radius = kR0 / sqrt(1.0 - pow(kEccentricity * sin(old_nav_est_ned.latitude),2.0)) *
        sqrt(pow(cos(old_nav_est_ned.latitude), 2.0) + pow(1.0 - kEccentricity*kEccentricity, 2.0) * pow(sin(old_nav_est_ned.latitude), 2.0)); // from (2.137)

    // Skew symmetric matrix of Earth rate
    Vector3 kOmega_ie_vec;
    kOmega_ie_vec << 0.0,0.0,kOmega_ie;
    Matrix3 Omega_ie = skewSymmetric(kOmega_ie_vec);

    // Determine error-state transition matrix using (14.50) (first-order approx)
    Eigen::Matrix<T,17,17> Phi_matrix = Eigen::Matrix<T,17,17>::Identity();
    Phi_matrix.template block<3,3>(0,0) -= Omega_ie * tor_i;
    Phi_matrix.template block<3,3>(0,12) = old_nav_est_ecef.C_b_e * tor_i;
    Phi_matrix.template block<3,3>(3,0) = - tor_i * skewSymmetric(Vector3(old_nav_est_ecef.C_b_e * imu_meas.f));
    Phi_matrix.template block<3,3>(3,3) -= 2.0 * Omega_ie * tor_i;
    Phi_matrix.template block<3,3>(3,6) = -tor_i * 2 * gravityEcef(old_nav_est_ecef.r_eb_e) /
        geocentric_radius * old_nav_est_ecef.r_eb_e.transpose() / old_nav_est_ecef.r_eb_e.norm();
    Phi_matrix.template block<3,3>(3,9) = old_nav_est_ecef.C_b_e * tor_i;
    Phi_matrix.template block<3,3>(6,3) = Matrix3::Identity() * tor_i;

    // Clock offset
    Phi_matrix(15,16) = tor_i;

    // Determine approximate system noise covariance matrix using (14.82)
    Eigen::Matrix<T,17,17> Q_prime_matrix = Eigen::Matrix<T,17,17>::Zero();
    Q_prime_matrix.template block<3,3>(0,0) = Matrix3::Identity() * tc_kf_config.gyro_noise_psd * tor_i;
    Q_prime_matrix.template block<3,3>(3,3) = Matrix3::Identity() * tc_kf_config.accel_noise_psd * tor_i;
    Q_prime_matrix.template block<3,3>(9,9) = Matrix3::Identity() * tc_kf_config.accel_bias_psd * tor_i;
    Q_prime_matrix.template block<3,3>(12,12) = Matrix3::Identity() * tc_kf_config.gyro_bias_psd * tor_i;

    // Clock offset, drift
    Q_prime_matrix(15,15) = tc_kf_config.clock_phase_psd * tor_i;
    Q_prime_matrix(16,16) = tc_kf_config.clock_freq_psd * tor_i;

    // Propagate state estimation error covariance matrix using (3.46)

    Eigen::Matrix<T,17,17> P_matrix;
    predictKF<17>(Phi_matrix, Q_prime_matrix, P_matrix_old, P_matrix);

    return P_matrix;
}

template<typename T>
typename Types<T>::StateEstEcef 
Navigation<T>::lcUpdateKFPosEcef (const PosMeasEcef & pos_meas, 
                                const StateEstEcef & state_est_prior,
                                const T & p_value) {
        
    // Set-up measurement matrix using (14.115)
    Eigen::Matrix<T,3,15> H_matrix = Eigen::Matrix<T,3,15>::Zero();
    H_matrix.template block<3,3>(0,6) = - Matrix3::Identity(); // Position

    // Set-up measurement noise covariance matrix
    Matrix3 R_matrix = pos_meas.cov_mat;

    // Formulate measurement innovations using (14.102), noting that zero
    // lever arm is assumed here
    Vector3 delta_z = pos_meas.r_eb_e - state_est_prior.nav_sol.r_eb_e;

    // Do error-state Kalman filter update
    Eigen::Matrix<T, 3, 3> S_matrix;
    Eigen::Matrix<T, 15, 1> x_est_new;
    Eigen::Matrix<T, 15, 15> P_matrix_post;
    bool valid_update = updateKF<15, 3>(delta_z,
                                        state_est_prior.P_matrix.template block<15,15>(0,0),
                                        H_matrix,
                                        R_matrix,
                                        p_value,
                                        S_matrix,
                                        x_est_new,
                                        P_matrix_post);

    // CLOSED-LOOP CORRECTION

    // Correct attitude, velocity, and position using (14.7-9)
    StateEstEcef state_est_post;
    state_est_post.valid = valid_update;
    state_est_post.nav_sol.time = state_est_prior.nav_sol.time;
    state_est_post.nav_sol.C_b_e = (Matrix3::Identity() - 
                                skewSymmetric(Vector3(x_est_new.template block<3,1>(0,0)))) * 
                                state_est_prior.nav_sol.C_b_e;
    state_est_post.nav_sol.v_eb_e = state_est_prior.nav_sol.v_eb_e - x_est_new.template block<3,1>(3,0);
    state_est_post.nav_sol.r_eb_e = state_est_prior.nav_sol.r_eb_e - x_est_new.template block<3,1>(6,0);
    state_est_post.P_matrix.template block<15,15>(0,0) = P_matrix_post;

    // Update IMU bias estimates
    state_est_post.acc_bias = state_est_prior.acc_bias + x_est_new.template block<3,1>(9,0);
    state_est_post.gyro_bias = state_est_prior.gyro_bias + x_est_new.template block<3,1>(12,0);

    // Update innovations and sigmas
    state_est_post.innovations_sigmas.clear();
    for(int i=0; i<delta_z.rows(); i++){
        state_est_post.innovations_sigmas.push_back(std::make_pair(delta_z(i), sqrt(S_matrix(i,i))));
    }
    
    // Return 
    return state_est_post;
}

template<typename T>
typename Types<T>::StateEstEcef 
Navigation<T>::lcUpdateKFGnssEcef (const GnssMeasurements & gnss_meas, 
                                    const StateEstEcef & state_est_prior,
                                    const GnssConfig & gnss_config,
                                    const T & p_value) {

    // Get position + velocity estimate with nlls
    GnssPosVelMeasEcef pos_vel_gnss_meas = Navigation<T>::gnssLsPositionVelocity(gnss_meas, 
                                                                state_est_prior.nav_sol.r_eb_e, 
                                                                state_est_prior.nav_sol.v_eb_e,
                                                                gnss_config);
        
    // Set-up measurement matrix using (14.115)
    Eigen::Matrix<T,6,15> H_matrix = Eigen::Matrix<T,6,15>::Zero();
    H_matrix.template block<3,3>(0,6) = - Matrix3::Identity(); // Position
    H_matrix.template block<3,3>(3,3) = - Matrix3::Identity(); // Velocity

    // Set-up measurement noise covariance matrix
    Eigen::Matrix<T,6,6> R_matrix = pos_vel_gnss_meas.cov_mat;

    // Formulate measurement innovations using (14.102), noting that zero
    // lever arm is assumed here
    Eigen::Matrix<T,6,1> delta_z;
    delta_z.template block<3,1>(0,0) = pos_vel_gnss_meas.r_ea_e - state_est_prior.nav_sol.r_eb_e;
    delta_z.template block<3,1>(3,0) = pos_vel_gnss_meas.v_ea_e - state_est_prior.nav_sol.v_eb_e;

    // Do error-state Kalman filter update
    Eigen::Matrix<T, 6, 6> S_matrix;
    Eigen::Matrix<T, 15, 1> x_est_new;
    Eigen::Matrix<T, 15, 15> P_matrix_post;
    bool valid_update = updateKF<15, 6>(delta_z,
                                        state_est_prior.P_matrix.template block<15,15>(0,0),
                                        H_matrix,
                                        R_matrix,
                                        p_value,
                                        S_matrix,
                                        x_est_new,
                                        P_matrix_post);

    // CLOSED-LOOP CORRECTION

    // Correct attitude, velocity, and position using (14.7-9)
    StateEstEcef state_est_post;
    state_est_post.valid = valid_update;
    state_est_post.nav_sol.time = state_est_prior.nav_sol.time;
    state_est_post.nav_sol.C_b_e = (Matrix3::Identity() - 
                                skewSymmetric(Vector3(x_est_new.template block<3,1>(0,0)))) * 
                                state_est_prior.nav_sol.C_b_e;
    state_est_post.nav_sol.v_eb_e = state_est_prior.nav_sol.v_eb_e - x_est_new.template block<3,1>(3,0);
    state_est_post.nav_sol.r_eb_e = state_est_prior.nav_sol.r_eb_e - x_est_new.template block<3,1>(6,0);
    state_est_post.P_matrix.template block<15,15>(0,0) = P_matrix_post;

    // Update IMU bias estimates
    state_est_post.acc_bias = state_est_prior.acc_bias + x_est_new.template block<3,1>(9,0);
    state_est_post.gyro_bias = state_est_prior.gyro_bias + x_est_new.template block<3,1>(12,0);

    // Update innovations and sigmas
    state_est_post.innovations_sigmas.clear();
    for(int i=0; i<delta_z.rows(); i++){
        state_est_post.innovations_sigmas.push_back(std::make_pair(delta_z(i), sqrt(S_matrix(i,i))));
    }

    return state_est_post;
}

template<typename T>
typename Types<T>::StateEstEcef 
Navigation<T>::tcUpdateKFGnssEcef (const GnssMeasurements & gnss_meas, 
                                const StateEstEcef & state_est_prior,
                                const T & tor_s,
                                const T & p_value) {

    // Compute predicted range and range rate measurements from prior state est
    Eigen::Matrix<T, Eigen::Dynamic, 2, 0, kMaxGnssSatellites, 2> pred_meas(gnss_meas.no_meas,2);
    // Line of sight unit vectors in ecef
    Eigen::Matrix<T, Eigen::Dynamic, 3, 0, kMaxGnssSatellites, 3> u_as_e_T(gnss_meas.no_meas,3);

    // For each satellite
    for(int j=0; j<gnss_meas.no_meas; j++) {

        // Predicted range, approximated
        Vector3 delta_r = gnss_meas.meas.template block<1,3>(j,2).transpose() - state_est_prior.nav_sol.r_eb_e;
        T range = delta_r.norm();

        // Calculate ecef rotation during signal transit time using (8.36)
        Matrix3 C_e_I;
        C_e_I << 1, kOmega_ie * range / kC, 0,
                -kOmega_ie * range / kC, 1, 0,
                0, 0, 1;

        // Predicted range, corrected
        // Convert satellite position in ecef frame at signal reception time, 
        // by taking into account earth rotation. This way you get the true range
        delta_r = C_e_I *  gnss_meas.meas.template block<1,3>(j,2).transpose() - state_est_prior.nav_sol.r_eb_e;
        range = delta_r.norm();
        // Also consider prior on clock offset (meters) 
        pred_meas(j,0) = range + state_est_prior.clock_offset + state_est_prior.clock_drift * tor_s;

        // Predict pseudo-range rate using (9.165)
        // As before, get satellite velocity at signal reception time
        // Skew symmetric matrix of Earth rate
        Vector3 kOmega_ie_vec;
        kOmega_ie_vec << 0.0,0.0,kOmega_ie;
        Matrix3 Omega_ie = skewSymmetric(kOmega_ie_vec);
        // Get line of sight unit vector
        u_as_e_T.template block<1,3>(j,0) = delta_r / range;
        T range_rate = u_as_e_T.template block<1,3>(j,0) * (C_e_I * (gnss_meas.meas.template block<1,3>(j,5).transpose() +
                                                        Omega_ie * gnss_meas.meas.template block<1,3>(j,2).transpose()) - 

                                                        (state_est_prior.nav_sol.v_eb_e + 
                                                        Omega_ie * state_est_prior.nav_sol.r_eb_e)); 
        // Also consider prior on clock drift (meters/seconds) 
        pred_meas(j,1) = range_rate + state_est_prior.clock_drift;
    }

    // Set-up measurement matrix using (14.126) - simplified Jacobian
    Eigen::Matrix<T, Eigen::Dynamic, 17, 0, 2 * kMaxGnssSatellites, 17> H_matrix = 
        Eigen::Matrix<T, Eigen::Dynamic, 17, 0, 2 * kMaxGnssSatellites, 17>::Zero(2*gnss_meas.no_meas, 17);

    // Ranges
    H_matrix.block(0,6,gnss_meas.no_meas,3) = u_as_e_T.block(0,0,gnss_meas.no_meas,3);
    H_matrix.block(0,15,gnss_meas.no_meas,1) = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Ones(gnss_meas.no_meas,1);
    // Range rates
    H_matrix.block(gnss_meas.no_meas,3, gnss_meas.no_meas,3) = u_as_e_T.block(0,0,gnss_meas.no_meas, 3);
    H_matrix.block(gnss_meas.no_meas,16, gnss_meas.no_meas,1) = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Ones(gnss_meas.no_meas,1);

    
    // Set-up measurement noise covariance matrix
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, 0, 2* kMaxGnssSatellites, 2* kMaxGnssSatellites> R_matrix = gnss_meas.cov_mat;

    // Formulate measurement innovations using (14.119)
    Eigen::Matrix<T, Eigen::Dynamic, 1, 0, 2* kMaxGnssSatellites, 1> delta_z = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(2*gnss_meas.no_meas,1);

    // Range innovations
    delta_z.block(0,0,gnss_meas.no_meas,1) = gnss_meas.meas.block(0,0,gnss_meas.no_meas,1) -
                                                pred_meas.block(0,0,gnss_meas.no_meas,1);
    // Range rates innovations
    delta_z.block(gnss_meas.no_meas,0,gnss_meas.no_meas,1) = gnss_meas.meas.block(0,1,gnss_meas.no_meas,1) -
                                                                pred_meas.block(0,1,gnss_meas.no_meas,1);


    // Do error-state Kalman filter update
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, 0, 2* kMaxGnssSatellites, 2* kMaxGnssSatellites> S_matrix;
    Eigen::Matrix<T, 17, 1> x_est_new;
    Eigen::Matrix<T, 17, 17> P_matrix_post;
    bool valid_update = updateKF<17, Eigen::Dynamic, 2* kMaxGnssSatellites>(delta_z,
                                                                            state_est_prior.P_matrix,
                                                                            H_matrix,
                                                                            R_matrix,
                                                                            p_value,
                                                                            S_matrix,
                                                                            x_est_new,
                                                                            P_matrix_post);

    // CLOSED-LOOP CORRECTION

    // Correct attitude, velocity, and position using (14.7-9)
    StateEstEcef state_est_post;
    state_est_post.valid = valid_update;
    state_est_post.nav_sol.time = state_est_prior.nav_sol.time;
    state_est_post.nav_sol.C_b_e = (Matrix3::Identity() - 
                                skewSymmetric(Vector3(x_est_new.template block<3,1>(0,0)))) * 
                                state_est_prior.nav_sol.C_b_e;
    state_est_post.nav_sol.v_eb_e = state_est_prior.nav_sol.v_eb_e - x_est_new.template block<3,1>(3,0);
    state_est_post.nav_sol.r_eb_e = state_est_prior.nav_sol.r_eb_e - x_est_new.template block<3,1>(6,0);
    state_est_post.P_matrix = P_matrix_post;

    // Update IMU bias estimates
    state_est_post.acc_bias = state_est_prior.acc_bias + x_est_new.template block<3,1>(9,0);
    state_est_post.gyro_bias = state_est_prior.gyro_bias + x_est_new.template block<3,1>(12,0);

    // Update clock offset + drift estimates
    state_est_post.clock_offset = state_est_prior.clock_offset +  x_est_new(15,0);
    state_est_post.clock_drift = state_est_prior.clock_drift + x_est_new(16,0);

    // Update innovations and sigmas
    state_est_post.innovations_sigmas.clear();
    for(int i=0; i<delta_z.rows(); i++){
        state_est_post.innovations_sigmas.push_back(std::make_pair(delta_z(i), sqrt(S_matrix(i,i))));
    }

    return state_est_post;
}

template<typename T>
typename Types<T>::StateEstEcef 
Navigation<T>::lcUpdateKFPosRotEcef(const PosRotMeasEcef & pos_rot_meas, 
                                    const StateEstEcef & state_est_prior,
                                    const T & p_value) {
        
    // Set-up measurement matrix using (14.115)
    Eigen::Matrix<T,6,15> H_matrix = Eigen::Matrix<T,6,15>::Zero();
    H_matrix.template block<3,3>(0,6) = - Matrix3::Identity(); // Position
    H_matrix.template block<3,3>(3,0) = - Matrix3::Identity(); // Rotation

    // Set-up measurement noise covariance matrix
    Eigen::Matrix<T, 6,6> R_matrix = pos_rot_meas.cov_mat;

    // Formulate measurement innovations using (14.102), noting that zero
    // lever arm is assumed here. See (14.151) for attitude int
    Eigen::Matrix<T,6,1> delta_z;
    delta_z.template block<3,1>(0,0) = pos_rot_meas.r_eb_e - state_est_prior.nav_sol.r_eb_e; // pos
    delta_z.template block<3,1>(3,0) = deSkew(Matrix3(pos_rot_meas.C_b_e * state_est_prior.nav_sol.C_b_e.transpose() - Matrix3::Identity()));// rot

    // Do error-state Kalman filter update
    Eigen::Matrix<T, 6, 6> S_matrix;
    Eigen::Matrix<T, 15, 1> x_est_new;
    Eigen::Matrix<T, 15, 15> P_matrix_post;
    bool valid_update = updateKF<15, 6>(delta_z,
                                        state_est_prior.P_matrix.template block<15,15>(0,0),
                                        H_matrix,
                                        R_matrix,
                                        p_value,
                                        S_matrix,
                                        x_est_new,
                                        P_matrix_post);
    
    // CLOSED-LOOP CORRECTION

    // Correct attitude, velocity, and position using (14.7-9)
    StateEstEcef state_est_post;
    state_est_post.valid = valid_update;
    state_est_post.nav_sol.time = state_est_prior.nav_sol.time;
    state_est_post.nav_sol.C_b_e = (Matrix3::Identity() - 
                                skewSymmetric(Vector3(x_est_new.template block<3,1>(0,0)))) * 
                                state_est_prior.nav_sol.C_b_e;
    state_est_post.nav_sol.v_eb_e = state_est_prior.nav_sol.v_eb_e - x_est_new.template block<3,1>(3,0);
    state_est_post.nav_sol.r_eb_e = state_est_prior.nav_sol.r_eb_e - x_est_new.template block<3,1>(6,0);
    state_est_post.P_matrix.template block<15,15>(0,0) = P_matrix_post;
    
    // Update IMU bias estimates
    state_est_post.acc_bias = state_est_prior.acc_bias + x_est_new.template block<3,1>(9,0);
    state_est_post.gyro_bias = state_est_prior.gyro_bias + x_est_new.template block<3,1>(12,0);

    // Update innovations and sigmas
    state_est_post.innovations_sigmas.clear();
    for(int i=0; i<delta_z.rows(); i++){
        state_est_post.innovations_sigmas.push_back(std::make_pair(delta_z(i), sqrt(S_matrix(i,i))));
    }

    return state_est_post;

}


template<typename T>
typename Types<T>::StateEstEcef 
Navigation<T>::lcUpdateKFOptFlow(const cv::Mat & img_0,
                                const cv::Mat & img_1,
                                const Matrix3 & K,
                                const Matrix3 & C_c_b,
                                const StateEstEcef & state_est_prior,
                                const T & p_value) {

    cv::Mat K_cv;
    cv::eigen2cv(K, K_cv);
    
    StateEstEcef state_est_post = state_est_prior;

    // Validate input
    if(img_0.empty() || img_1.empty() || img_0.size() != img_1.size() || img_0.type() != img_1.type()) {
        state_est_post.valid = false;
        return state_est_post;
    }

    // Convert to grayscale
    cv::Mat gray_0, gray_1;
    if (img_0.channels() == 1) {
        gray_0 = img_0;
        gray_1 = img_1;
    } else {
        cv::cvtColor(img_0, gray_0, cv::COLOR_BGR2GRAY);
        cv::cvtColor(img_1, gray_1, cv::COLOR_BGR2GRAY);
    }

    // Detect keypoints in the first image
    std::vector<cv::Point2f> points0;
    cv::goodFeaturesToTrack(gray_0, points0, 500, 0.01, 10);

    // If not enough points, return
    if(points0.size() < 8) {
        LOG(WARNING) << "OPF update failed: not enough goodFeaturesToTrack";
        state_est_post.valid = false;
        return state_est_post;
    }

    // Track keypoints using optical flow
    std::vector<cv::Point2f> points1;
    std::vector<unsigned char> status;
    std::vector<float> err;
    cv::calcOpticalFlowPyrLK(gray_0, gray_1, points0, points1, status, err);

    // Filter valid tracked points
    std::vector<cv::Point2f> inlier_pts0, inlier_pts1;
    inlier_pts0.reserve(points0.size());
    inlier_pts1.reserve(points1.size());
    for (size_t i = 0; i < status.size(); ++i) {
        if (status[i]) {
            inlier_pts0.push_back(points0[i]);
            inlier_pts1.push_back(points1[i]);
        }
    }

    // If not enough inliers, return not valid
    if (inlier_pts0.size() < 8) {
        LOG(WARNING) << "OPF update failed: not enough calcOpticalFlowPyrLK inliers";
        state_est_post.valid = false;
        return state_est_post;
    }


    // Estimate Essential matrix: p1 * K^-T * E * K^-1 * p0 = 0
    cv::Mat mask;  // inliers
    cv::Mat E_cv = cv::findEssentialMat(inlier_pts0,
                                        inlier_pts1,
                                        K_cv,
                                        cv::USAC_MAGSAC,
                                        0.99,          // confidence
                                        1.0,            // epipolar threshold
                                        mask);

    if (E_cv.empty()) {
        LOG(WARNING) << "OPF update failed: not enough findEssentialMat inliers";
        state_est_post.valid = false;
        return state_est_post;
    }

    // Decompose into R_0_1 and t_10_1
    cv::Mat R_cv, t_cv;
    int inliers = cv::recoverPose(E_cv,
                                inlier_pts0,
                                inlier_pts1,
                                K_cv,
                                R_cv,
                                t_cv,
                                cv::USAC_MAGSAC,
                                mask);

    if (inliers < 8) {
        LOG(WARNING) << "OPF update failed: not enough recoverPose inliers";
        state_est_post.valid = false;
        return state_est_post;
    }

    // ========== Project and display the translation vector ==========
    // t_cv is the translation direction (up to scale) expressed in the first camera frame.
    // We project a point along this direction into the first image and draw an arrow from
    // the principal point to the projected point for debugging/visualization.
    try {
    // Ensure double precision for projection math
    cv::Mat t_cam, Kd;
    t_cv.convertTo(t_cam, CV_64F);
    K_cv.convertTo(Kd, CV_64F);

    // principal point (cx, cy)
    double cx = Kd.at<double>(0,2);
    double cy = Kd.at<double>(1,2);
    cv::Point2d pp(cx, cy);

    // Choose a scale along the translation direction so the projected point falls inside the image
    // A heuristic scale: some multiple of the focal length (so it's visible). If depth is negative,
    // flip sign so the arrow points into the scene.
    double focal = 0.5 * (Kd.at<double>(0,0) + Kd.at<double>(1,1));
    double scale = focal * 5.0; // 5 focal lengths out

    // If t points away (z < 0), flip to point forward
    double tz = t_cam.at<double>(2,0);
    if (tz < 0) {
    t_cam = -t_cam;
    tz = -tz;
    }

    cv::Mat P = t_cam * scale; // 3x1 point in camera frame
    cv::Mat ph = Kd * P; // projected homogeneous

    if (std::abs(ph.at<double>(2,0)) > 1e-8) {
    cv::Point2d proj_pt(ph.at<double>(0,0) / ph.at<double>(2,0), 
                        ph.at<double>(1,0) / ph.at<double>(2,0));

    // Draw on a copy of the image
    cv::Mat img_arrow = img_1.clone();
    // If grayscale, convert to color for visualization
    if (img_arrow.channels() == 1) cv::cvtColor(img_arrow, img_arrow, cv::COLOR_GRAY2BGR);

    // Draw principal point
    cv::circle(img_arrow, pp, 3, cv::Scalar(255,0,0), -1);
    // Draw arrow from principal point to projected point
    cv::arrowedLine(img_arrow, pp, proj_pt, cv::Scalar(0,0,255), 2, cv::LINE_AA, 0, 0.15);

    cv::imshow("TranslationVector", img_arrow);
    cv::waitKey(1);
    }
    } catch (const std::exception &e) {
    LOG(WARNING) << "Failed to draw translation vector: " << e.what();
    } catch (...) {
    LOG(WARNING) << "Failed to draw translation vector: unknown error";
    }

    // ========== Plot inlier matches for debugging ========== 

    cv::Mat img_matches;
    std::vector<cv::KeyPoint> kp0, kp1;
    cv::KeyPoint::convert(inlier_pts0, kp0);
    cv::KeyPoint::convert(inlier_pts1, kp1);
    std::vector<cv::DMatch> matches;
    for (int i = 0; i < inlier_pts0.size(); ++i) {
        if (mask.at<unsigned char>(i)) {
            matches.push_back(cv::DMatch(i, i, 0));
        }
    }
    cv::drawMatches(img_0, kp0,
                    img_1, kp1,
                    matches,
                    img_matches,
                    cv::Scalar(0, 255, 0), // matchColor
                    cv::Scalar::all(-1),   // singlePointColor
                    std::vector<char>(),   // matchesMask
                    cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS);
    cv::imshow("Matches", img_matches);
    cv::waitKey(1);

    // ====================================================

    // Convert to Eigen
    Matrix3 R;
    cv::cv2eigen(R_cv, R);
    Vector3 t;
    cv::cv2eigen(t_cv, t);
    Matrix3 E;
    cv::cv2eigen(E_cv, E);

    // Fundamental matrix
    Matrix3 F = K.transpose().inverse() * E * K.inverse();

    // Extract and normalize translation vector
    Vector3 vel_meas_opf = t.normalized();

    // It's in camera frame: rotate to body
    vel_meas_opf = C_c_b * vel_meas_opf;
    // Then rotate to ecef
    vel_meas_opf = state_est_prior.nav_sol.C_b_e * vel_meas_opf;

    // Set-up measurement matrix: jacobian of normalized velocity vector v / ||v|| w.r.t. state
    Eigen::Matrix<T,3,15> H_matrix = Eigen::Matrix<T,3,15>::Zero();
    Vector3 vel_est = state_est_prior.nav_sol.v_eb_e;
    T vel_est_norm = vel_est.norm();
    H_matrix.template block<3,3>(0,3) = - 1 / vel_est_norm * 
                                        (Eigen::Matrix<T,3,3>::Identity() + 
                                        1 / pow(vel_est_norm,2) * (vel_est * vel_est.transpose()));
    
    // Set-up measurement noise covariance matrix: compute from epipolar residuals
    // R = J * J' * cov(epipolar_residuals) --> sum for each residual

    // Compute epipolar residuals sigma
    const size_t m = inlier_pts0.size();
    Eigen::ArrayXd epipolar_residuals(m);
    for (size_t i = 0; i < m; ++i) {
        Vector3 x0 = Vector3(inlier_pts0[i].x, inlier_pts0[i].y, 1);
        Vector3 x1 = Vector3(inlier_pts1[i].x, inlier_pts1[i].y, 1);
        T epipolar_residual = x0.transpose() * F * x1;
        epipolar_residuals(i) = epipolar_residual;
    }

    // Suppose zero-mean
    // T epipolar_residuals_std = std::sqrt(epipolar_residuals.square().sum()/(epipolar_residuals.size()-1));
    T epipolar_residuals_std = 50;

    // Compute measurement noise covariance matrix
    Eigen::Matrix<T, 3, 3> R_matrix;
    for (size_t i = 0; i < m; ++i) {
        // Jacobian of i-th epipolar residual (1D) w.r.t. translation vector (3D): It's a 1x3 vec
        Vector3 x0 = Vector3(inlier_pts1[i].x, inlier_pts1[i].y, 1);
        Vector3 x1 = K.colPivHouseholderQr().solve(x0);
        Eigen::Matrix<T, 1, 3> J_r_t = (skewSymmetric(x1) * R * x0).transpose();
        R_matrix += (J_r_t.transpose() * J_r_t);
    }

    // Invert and mult with residual variance
    R_matrix = R_matrix.template selfadjointView<Eigen::Upper>().llt().solve(Eigen::Matrix<T, 3, 3>::Identity()) * pow(epipolar_residuals_std, 2);

    // Formulate measurement innovations
    Eigen::Matrix<T,3,1> delta_z = vel_meas_opf - vel_est.normalized();

    // Perform error-state Kalman filter update
    Eigen::Matrix<T, 3, 3> S_matrix;
    Eigen::Matrix<T, 15, 1> x_est_new;
    Eigen::Matrix<T, 15, 15> P_matrix_post;
    bool valid_update = updateKF<15, 3>(delta_z,
                                        state_est_prior.P_matrix.template block<15,15>(0,0),
                                        H_matrix,
                                        R_matrix,
                                        p_value,
                                        S_matrix,
                                        x_est_new,
                                        P_matrix_post);

    // Correct attitude, velocity, and position using (14.7-9)
    state_est_post.valid = valid_update;
    state_est_post.nav_sol.time = state_est_prior.nav_sol.time;
    state_est_post.nav_sol.C_b_e = (Matrix3::Identity() - 
                                skewSymmetric(Vector3(x_est_new.template block<3,1>(0,0)))) * 
                                state_est_prior.nav_sol.C_b_e;
    state_est_post.nav_sol.v_eb_e = state_est_prior.nav_sol.v_eb_e - x_est_new.template block<3,1>(3,0);
    state_est_post.nav_sol.r_eb_e = state_est_prior.nav_sol.r_eb_e - x_est_new.template block<3,1>(6,0);
    state_est_post.P_matrix.template block<15,15>(0,0) = P_matrix_post;
    
    // Update IMU bias estimates
    state_est_post.acc_bias = state_est_prior.acc_bias + x_est_new.template block<3,1>(9,0);
    state_est_post.gyro_bias = state_est_prior.gyro_bias + x_est_new.template block<3,1>(12,0);

    // Update innovations and sigmas
    state_est_post.innovations_sigmas.clear();
    for(int i=0; i<delta_z.rows(); i++){
        state_est_post.innovations_sigmas.push_back(std::make_pair(delta_z(i), sqrt(S_matrix(i,i))));
    }

    return state_est_post;
}

template<typename T>
typename Types<T>::GnssLsPosVelClock 
Navigation<T>::gnssLsPositionVelocityClock(const GnssMeasurements & gnss_measurements,
                                            const Vector3 & prior_r_ea_e,
                                            const Vector3 & prior_v_ea_e) {
    GnssLsPosVelClock est_pos_vel;
    // POSITION AND CLOCK OFFSET
    Vector4 x_pred, x_est;
    x_pred.template segment<3>(0) = prior_r_ea_e;
    x_pred(3) = 0;
    T test_convergence = 1.0;
    T epsilon_convergence = 0.0001;  
    int max_iters = 10000;     
    int count_iters = 0;                          
    while (test_convergence > epsilon_convergence && count_iters < max_iters) {
        count_iters++;
        Eigen::Matrix<T, Eigen::Dynamic, 4, 0, kMaxGnssSatellites> H_matrix(gnss_measurements.no_meas, 4);
        Eigen::Matrix<T, Eigen::Dynamic, 1, 0, kMaxGnssSatellites> pred_meas(gnss_measurements.no_meas);
        for (int j = 0; j < gnss_measurements.no_meas; ++j) {
            Vector3 delta_r = gnss_measurements.meas.template block<1,3>(j,2).transpose() - x_pred.template segment<3>(0);
            T approx_range = delta_r.norm();
            Matrix3 C_e_I;
            C_e_I << 1, kOmega_ie * approx_range / kC, 0,
                     -kOmega_ie * approx_range / kC, 1, 0,
                     0, 0, 1;
            delta_r = C_e_I * gnss_measurements.meas.template block<1,3>(j,2).transpose() - x_pred.template segment<3>(0);
            T range = delta_r.norm();
            pred_meas(j) = range + x_pred(3);
            H_matrix.template block<1,3>(j,0) = -delta_r.transpose() / range;
            H_matrix(j,3) = 1;
        }

        // Gauss-Newton: J'J * Dx = J' * cost -> Ax = b -> x = A.llt().solve(b)
        // x = x + Dx
        auto A = H_matrix.transpose() * H_matrix;
        auto b = H_matrix.transpose() * (gnss_measurements.meas.col(0).head(gnss_measurements.no_meas) - pred_meas.head(gnss_measurements.no_meas));
        auto delta_x = A.template selfadjointView<Eigen::Upper>().llt().solve(b);
        x_est = x_pred + delta_x;
        test_convergence = (x_est - x_pred).norm();
        x_pred = x_est;
    }

    est_pos_vel.r_ea_e = x_est.template segment<3>(0);
    est_pos_vel.clock(0) = x_est(3);
    // VELOCITY AND CLOCK DRIFT
    Matrix3 Omega_ie = skewSymmetric(Vector3(0, 0, kOmega_ie));
    x_pred.template segment<3>(0) = prior_v_ea_e;
    x_pred(3) = 0;
    test_convergence = 1.0;
    count_iters = 0;
    while (test_convergence > epsilon_convergence && count_iters < max_iters) {
        count_iters++;
        Eigen::Matrix<T, Eigen::Dynamic, 4, 0, kMaxGnssSatellites> H_matrix(gnss_measurements.no_meas, 4);
        Eigen::Matrix<T, Eigen::Dynamic, 1, 0, kMaxGnssSatellites> pred_meas(gnss_measurements.no_meas);
        for (int j = 0; j < gnss_measurements.no_meas; ++j) {
            Vector3 delta_r = gnss_measurements.meas.template block<1,3>(j,2).transpose() - est_pos_vel.r_ea_e;
            T approx_range = delta_r.norm();
            Matrix3 C_e_I;
            C_e_I << 1, kOmega_ie * approx_range / kC, 0,
                     -kOmega_ie * approx_range / kC, 1, 0,
                     0, 0, 1;
            delta_r = C_e_I * gnss_measurements.meas.template block<1,3>(j,2).transpose() - est_pos_vel.r_ea_e;
            T range = delta_r.norm();
            Vector3 u_as_e = delta_r / range;
            Vector3 sat_velocity = gnss_measurements.meas.template block<1,3>(j,5).transpose();
            Vector3 sat_position = gnss_measurements.meas.template block<1,3>(j,2).transpose();
            T range_rate = u_as_e.transpose() * (C_e_I * (sat_velocity + Omega_ie * sat_position) - (x_pred.template segment<3>(0) + Omega_ie * est_pos_vel.r_ea_e));
            pred_meas(j) = range_rate + x_pred(3);
            H_matrix.template block<1,3>(j,0) = -u_as_e.transpose();
            H_matrix(j,3) = 1;
        }

        // Gauss-Newton: J'J * Dx = J' * cost -> Ax = b -> x = A.llt().solve(b)
        // x = x + Dx
        auto A = H_matrix.transpose() * H_matrix;
        auto b = H_matrix.transpose() * (gnss_measurements.meas.col(1).head(gnss_measurements.no_meas) - pred_meas.head(gnss_measurements.no_meas));
        auto delta_x = A.template selfadjointView<Eigen::Upper>().llt().solve(b);
        x_est = x_pred + delta_x;

        test_convergence = (x_est - x_pred).norm();
        x_pred = x_est;
    }
    est_pos_vel.v_ea_e = x_est.template segment<3>(0);
    est_pos_vel.clock(1) = x_est(3);
    return est_pos_vel;
}

template<typename T>
typename Types<T>::GnssPosVelMeasEcef 
Navigation<T>::gnssLsPositionVelocity(const GnssMeasurements & gnss_measurements,
                                    const Vector3 & prior_r_ea_e,
                                    const Vector3 & prior_v_ea_e,
                                    const GnssConfig & gnss_config) {
    
    GnssLsPosVelClock est_pos_vel = Navigation<T>::gnssLsPositionVelocityClock(gnss_measurements, prior_r_ea_e, prior_v_ea_e);
    GnssPosVelMeasEcef pos_vel_gnss_meas_ecef;
    pos_vel_gnss_meas_ecef.r_ea_e = est_pos_vel.r_ea_e;
    pos_vel_gnss_meas_ecef.v_ea_e = est_pos_vel.v_ea_e;
    pos_vel_gnss_meas_ecef.cov_mat = Eigen::Matrix<T,6,6>::Identity();
    pos_vel_gnss_meas_ecef.cov_mat.template block<3,3>(0,0) = pow(gnss_config.lc_pos_sd,2.0) * Matrix3::Identity();
    pos_vel_gnss_meas_ecef.cov_mat.template block<3,3>(3,3) = pow(gnss_config.lc_vel_sd,2.0) * Matrix3::Identity();
    return pos_vel_gnss_meas_ecef;
}

template<typename T>
Eigen::Matrix<T,17,17> 
Navigation<T>::initializePMmatrix(const KfConfig & kf_config) {
    Eigen::Matrix<T,17,17> P_matrix;
    // Initialize error covariance matrix
    P_matrix = Eigen::Matrix<T,17,17>::Zero();
    P_matrix.template block<3,3>(0,0) = Matrix3::Identity() * pow(kf_config.init_att_unc,2); // attitude error
    P_matrix.template block<3,3>(3,3) = Matrix3::Identity() * pow(kf_config.init_vel_unc,2); // vel error
    P_matrix.template block<3,3>(6,6) = Matrix3::Identity() * pow(kf_config.init_pos_unc,2); // pos error
    P_matrix.template block<3,3>(9,9) = Matrix3::Identity() * pow(kf_config.init_b_a_unc,2); // acc bias error
    P_matrix.template block<3,3>(12,12) = Matrix3::Identity() * pow(kf_config.init_b_g_unc,2); // gyro bias error
    P_matrix(15,15) = pow(kf_config.init_clock_offset_unc,2); // clock offset error
    P_matrix(16,16) = pow(kf_config.init_clock_drift_unc,2); // clock drift error
    return P_matrix;
}

template<typename T>
typename Types<T>::StateEstEcef 
Navigation<T>::initStateFromGroundTruth(const NavSolutionEcef & true_nav_ecef, 
                                        const KfConfig & kf_config, 
                                        const GnssMeasurements & gnss_meas, 
                                        std::mt19937 & gen) {
    StateEstEcef state_est_ecef;
    state_est_ecef.valid = true;
    state_est_ecef.nav_sol = true_nav_ecef;
    std::normal_distribution att_d{T(0.0), kf_config.init_att_unc};
    std::normal_distribution vel_d{T(0.0), kf_config.init_vel_unc};
    std::normal_distribution pos_d{T(0.0), kf_config.init_pos_unc};
    state_est_ecef.nav_sol.C_b_e = true_nav_ecef.C_b_e * eulerToDcm(Vector3(att_d(gen), att_d(gen), att_d(gen)));
    state_est_ecef.nav_sol.r_eb_e += Vector3(pos_d(gen), pos_d(gen), pos_d(gen));
    state_est_ecef.nav_sol.v_eb_e += Vector3(vel_d(gen), vel_d(gen), vel_d(gen));
    state_est_ecef.acc_bias = Vector3::Zero();
    state_est_ecef.gyro_bias = Vector3::Zero();
    // Error covariance matrix
    state_est_ecef.P_matrix = Navigation<T>::initializePMmatrix(kf_config);
    // Init clock states using NLLS solver
    GnssLsPosVelClock gnss_pos_vel_clock_est = Navigation<T>::gnssLsPositionVelocityClock(gnss_meas, state_est_ecef.nav_sol.r_eb_e, state_est_ecef.nav_sol.v_eb_e);
    state_est_ecef.clock_offset = gnss_pos_vel_clock_est.clock(0);
    state_est_ecef.clock_drift = gnss_pos_vel_clock_est.clock(1);
    return state_est_ecef;
}

template struct Navigation<double>;
template struct Navigation<float>;

};