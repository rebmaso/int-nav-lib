#include "navigation.h"

namespace intnavlib {

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