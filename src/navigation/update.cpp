#include "navigation.h"

namespace intnavlib {

// IMPLEMENTATION OF MEMBER TEMPLATES

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
    // P_matrix_post = (Eigen::Matrix<T, n_x, n_x>::Identity() - K_matrix * H_matrix) * P_matrix;

    // Joseph form: so we ensure positive definiteness of P
    Eigen::Matrix<T, n_x, n_x> P_matrix_post_fac_1 = (Eigen::Matrix<T, n_x, n_x>::Identity() - K_matrix * H_matrix);
    P_matrix_post = P_matrix_post_fac_1 * P_matrix * P_matrix_post_fac_1.transpose() + K_matrix * R_matrix * K_matrix.transpose();

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
typename Types<T>::StateEstEcef 
Navigation<T>::lcUpdateKFPosEcef (const PosMeasEcef & pos_meas, 
                                const StateEstEcef & state_est_prior,
                                const T & p_value) {

    // A priori error state is always zero in closed loop filter
    Eigen::Matrix<T,15,1> x_est_propagated = Eigen::Matrix<T,15,1>::Zero();
        
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
Navigation<T>::lcUpdateKFPosRotEcef (const PosRotMeasEcef & pos_rot_meas, 
                                    const StateEstEcef & state_est_prior,
                                    const T & p_value) {

    // A priori error state is always zero in closed loop filter
    Eigen::Matrix<T,15,1> x_est_propagated = Eigen::Matrix<T,15,1>::Zero();
        
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

template struct Navigation<double>;
template struct Navigation<float>;

};