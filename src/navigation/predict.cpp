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

template struct Navigation<double>;
template struct Navigation<float>;

};