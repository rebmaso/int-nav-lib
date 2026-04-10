#include "helpers.h"

namespace intnavlib {

template <typename T>
typename Types<T>::EvalDataEcef 
Helpers<T>::getEvalDataEcef(const StateEstEcef & state_est_ecef, 
                            const NavSolutionEcef & true_nav_ecef) {

    EvalDataEcef eval_data;

    // Errors
    eval_data.time = state_est_ecef.nav_sol.time;
    eval_data.delta_r_eb_e = state_est_ecef.nav_sol.r_eb_e - true_nav_ecef.r_eb_e;
    eval_data.delta_v_eb_e = state_est_ecef.nav_sol.v_eb_e - true_nav_ecef.v_eb_e;
    eval_data.delta_rot_eb_e = deSkew(state_est_ecef.nav_sol.C_b_e * true_nav_ecef.C_b_e.transpose() - Matrix3::Identity());

    // Sigmas
    eval_data.sigma_delta_r_eb_e << sqrt(state_est_ecef.P_matrix(6,6)), sqrt(state_est_ecef.P_matrix(7,7)), sqrt(state_est_ecef.P_matrix(8,8));
    eval_data.sigma_delta_v_eb_e << sqrt(state_est_ecef.P_matrix(3,3)), sqrt(state_est_ecef.P_matrix(4,4)), sqrt(state_est_ecef.P_matrix(5,5));
    eval_data.sigma_delta_rot_eb_e << sqrt(state_est_ecef.P_matrix(0,0)), sqrt(state_est_ecef.P_matrix(1,1)), sqrt(state_est_ecef.P_matrix(2,2));
    
    // These are just error-state estimates and their estimated sigmas
    eval_data.delta_b_a = state_est_ecef.acc_bias;
    eval_data.delta_b_g = state_est_ecef.gyro_bias;
    eval_data.delta_clock_offset = state_est_ecef.clock_offset;
    eval_data.delta_clock_drift = state_est_ecef.clock_drift;
    eval_data.sigma_delta_b_a << sqrt(state_est_ecef.P_matrix(9,9)), sqrt(state_est_ecef.P_matrix(10,10)), sqrt(state_est_ecef.P_matrix(11,11));
    eval_data.sigma_delta_b_g << sqrt(state_est_ecef.P_matrix(12,12)), sqrt(state_est_ecef.P_matrix(13,13)), sqrt(state_est_ecef.P_matrix(14,14));
    eval_data.sigma_delta_clock_offset = sqrt(state_est_ecef.P_matrix(15,15));
    eval_data.sigma_delta_clock_drift = sqrt(state_est_ecef.P_matrix(16,16));

    // Innovations and sigmas
    eval_data.innovations_sigmas = state_est_ecef.innovations_sigmas;

    return eval_data;
}

template <typename T>
typename Helpers<T>::Vector3 
Helpers<T>::gravityEcef(const typename Helpers<T>::Vector3 & r_eb_e) {
    T mag_r = r_eb_e.norm();
    Vector3 g = Vector3::Zero();
    // If the input position is 0,0,0, produce a dummy output
    if (mag_r >= kEpsilon)
    // Calculate gravitational acceleration using (2.142)
    {
        T z_scale = 5.0 * pow(r_eb_e(2) / mag_r,2.0);
        Vector3 gamma_1;
        gamma_1 << (1.0 - z_scale) * r_eb_e(0), 
                    (1.0 - z_scale) * r_eb_e(1),
                    (3.0 - z_scale) * r_eb_e(2);
        Vector3 gamma;
        gamma = (-kGravConst / pow(mag_r,3.0)) *
                (r_eb_e + 1.5 * kJ2 * 
                pow(kR0 / mag_r,2.0) * gamma_1);
        // Add centripetal acceleration using (2.133)
        g(0) = gamma(0) + pow(kOmega_ie,2.0) * r_eb_e(0);
        g(1) = gamma(1) + pow(kOmega_ie,2.0) * r_eb_e(1);
        g(2) = gamma(2);
    }
    return g;
}

template <typename T>
typename Helpers<T>::Matrix3 
Helpers<T>::skewSymmetric(const typename Helpers<T>::Vector3 & a) {
    Matrix3 S;
    S << 0.0, -a(2),  a(1),
      a(2),     0.0, -a(0),
     -a(1),  a(0),     0.0;
    return S;
}

template <typename T>
typename Helpers<T>::Vector3 
Helpers<T>::deSkew(const typename Helpers<T>::Matrix3 & S){
    Vector3 a;
    a << -S(1,2), S(0,2), -S(0,1);
    return a; 
}

template <typename T>
std::string 
Helpers<T>::getCurrentDateTime() {
    auto now = std::time(nullptr);
    std::tm tm_now;
    localtime_r(&now, &tm_now); // Use localtime_s on Windows or localtime_r on Unix-like systems
    std::ostringstream oss;
    oss << std::put_time(&tm_now, "%Y%m%d_%H%M%S");
    return oss.str();
}

template <typename T>
typename Helpers<T>::Vector2
Helpers<T>::radiiOfCurvature(T L) {
    // Calculate meridian radius of curvature using (2.105)
    T temp = 1.0 - pow((kEccentricity * sin(L)),2.0); 
    T R_N = kR0 * (1.0 - pow(kEccentricity,2.0)) / pow(temp,1.5);
    // Calculate transverse radius of curvature using (2.105)
    T R_E = kR0 / sqrt(temp);
    Vector2 radii;
    radii << R_N, R_E;
    return radii;
}

template <typename T>
typename Types<T>::ErrorsNed 
Helpers<T>::calculateErrorsNed(const NavSolutionNed & true_nav_sol, 
                                const NavSolutionNed & est_nav_sol){
    // Position error calculation
    Vector2 radii = radiiOfCurvature(true_nav_sol.latitude);
    T R_N = radii(0);
    T R_E = radii(1);
    Vector3 delta_r_eb_n;
    delta_r_eb_n(0) = (est_nav_sol.latitude - true_nav_sol.latitude) * (R_N + true_nav_sol.height);
    delta_r_eb_n(1) = (est_nav_sol.longitude - true_nav_sol.longitude) * (R_E + true_nav_sol.height) * std::cos(true_nav_sol.latitude);
    delta_r_eb_n(2) = -(est_nav_sol.height - true_nav_sol.height);
    // Velocity error calculation
    Vector3 delta_v_eb_n = est_nav_sol.v_eb_n - true_nav_sol.v_eb_n;
    // Attitude error calculation
    Matrix3 delta_C_b_n = est_nav_sol.C_b_n * true_nav_sol.C_b_n.transpose();
    Vector3 delta_rot_nb_n = -dcmToEuler(delta_C_b_n);
    ErrorsNed errors_ned;
    errors_ned.time = true_nav_sol.time;
    errors_ned.delta_r_eb_n = delta_r_eb_n;
    errors_ned.delta_v_eb_n = delta_v_eb_n;
    errors_ned.delta_rot_nb_n = delta_rot_nb_n;
    return errors_ned;
}

template <typename T>
typename Types<T>::ImuErrors 
Helpers<T>::tacticalImuErrors(){
    ImuErrors imu_errors;
    imu_errors.b_a << 900.0,-1300.0,800.0;
    imu_errors.b_a = imu_errors.b_a * kMuGToMetersPerSecondSquared;
    imu_errors.b_g << -9.0, 13.0, -8.0;
    imu_errors.b_g = imu_errors.b_g * kDegToRad / 3600.0;
    imu_errors.M_a << 500.0, -300.0, 200.0,
            -150.0, -600.0, 250.0,
            -250.0,  100.0, 450.0;
    imu_errors.M_a = imu_errors.M_a * 1.0e-6;
    imu_errors.M_g << 400.0, -300.0,  250.0,
            0.0, -300.0, -150.0,
            0.0,    0.0, -350.0; 
    imu_errors.M_g = imu_errors.M_g * 1.0e-6;
    imu_errors.G_g << 0.9, -1.1, -0.6,
            -0.5,  1.9, -1.6,
            0.3,  1.1, -1.3;
    imu_errors.G_g = imu_errors.G_g * kDegToRad / (3600.0 * 9.80665);  
    imu_errors.accel_noise_root_psd = 100.0 * kMuGToMetersPerSecondSquared;
    imu_errors.gyro_noise_root_psd = 0.01 * kDegToRad / 60.0;
    imu_errors.accel_quant_level = 1.0e-2;
    imu_errors.gyro_quant_level = 2.0e-4;
    return imu_errors;
}

template <typename T>
typename Types<T>::GnssConfig 
Helpers<T>::defaultGnssConfig(){
    GnssConfig gnss_config;
    gnss_config.epoch_interval = 0.5;
    gnss_config.init_est_r_ea_e = Vector3::Zero();
    gnss_config.no_sat = 30.0;
    gnss_config.r_os = 2.656175E7;
    gnss_config.inclination = 55.0;
    gnss_config.const_delta_lambda = 0.0;
    gnss_config.const_delta_t = 0.0;
    gnss_config.mask_angle = 10.0;
    gnss_config.sis_err_sd = 1.0;
    gnss_config.zenith_iono_err_sd = 2.0;
    gnss_config.zenith_trop_err_sd = 0.2;
    gnss_config.code_track_err_sd = 1.0;
    gnss_config.rate_track_err_sd = 0.02;
    gnss_config.rx_clock_offset = 10000.0;
    gnss_config.rx_clock_drift = 100.0;
    gnss_config.lc_pos_sd = 2.5;
    gnss_config.lc_vel_sd = 0.1;
    gnss_config.pseudo_range_sd = 2.5;
    gnss_config.range_rate_sd = 0.1;
    return gnss_config;
}

template <typename T>
typename Types<T>::KfConfig 
Helpers<T>::tacticalImuKFConfig(){
    KfConfig kf_config;
    kf_config.init_att_unc = kDegToRad * 1.0;
    kf_config.init_vel_unc = 0.1;
    kf_config.init_pos_unc = 10.0;
    kf_config.init_b_a_unc = 1000.0 * kMuGToMetersPerSecondSquared;
    kf_config.init_b_g_unc = 10.0 * kDegToRad / 3600.0;
    kf_config.init_clock_offset_unc = 10.0;
    kf_config.init_clock_drift_unc = 0.1;
    kf_config.gyro_noise_psd = pow(0.02 * kDegToRad / 60.0, 2.0);
    kf_config.accel_noise_psd = pow(200.0 * kMuGToMetersPerSecondSquared, 2.0);
    kf_config.accel_bias_psd = 1.0E-7;
    kf_config.gyro_bias_psd = 2.0E-12;
    kf_config.clock_freq_psd = 1;
    kf_config.clock_phase_psd = 1;
    return kf_config;
}

template struct Helpers<float>;
template struct Helpers<double>;

};