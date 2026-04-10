#include "simulation.h"

namespace intnavlib {

template<typename T>
typename Types<T>::SatPosVel
Simulation<T>::satellitePositionsAndVelocities(const T & time, const GnssConfig & gnss_config) {

    SatPosVel gnssPosVel;

    // Convert inclination angle to radians
    T inclination = gnss_config.inclination * kDegToRad; 

    // Determine orbital angular rate
    T omega_is = std::sqrt(kGravConst / std::pow(gnss_config.r_os, 3));

    // Determine constellation time
    T const_time = time + gnss_config.const_delta_t;

    // Resize output matrices
    int no_sat = static_cast<int>(gnss_config.no_sat);

    gnssPosVel.sat_r_es_e = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(no_sat, 3);
    gnssPosVel.sat_v_es_e = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(no_sat, 3);

    // Loop over satellites
    for (int j = 0; j < no_sat; ++j) {
        // Argument of latitude
        T u_os_o = 2 * M_PI * j / no_sat + omega_is * const_time;
        
        // Satellite position in the orbital frame
        Vector3 r_os_o;
        r_os_o << gnss_config.r_os * std::cos(u_os_o),
                    gnss_config.r_os * std::sin(u_os_o), 
                    0;

        // Longitude of the ascending node
        T Omega = (M_PI * (j % 6) / 3 + gnss_config.const_delta_lambda * kDegToRad) - kOmega_ie * const_time;

        // ECEF Satellite Position
        gnssPosVel.sat_r_es_e(j, 0) = r_os_o(0) * std::cos(Omega) - r_os_o(1) * std::cos(inclination) * std::sin(Omega);
        gnssPosVel.sat_r_es_e(j, 1) = r_os_o(0) * std::sin(Omega) + r_os_o(1) * std::cos(inclination) * std::cos(Omega);
        gnssPosVel.sat_r_es_e(j, 2) = r_os_o(1) * std::sin(inclination);

        // Satellite velocity in the orbital frame
        Vector3 v_os_o;
        v_os_o << -gnss_config.r_os * omega_is * std::sin(u_os_o),
                    gnss_config.r_os * omega_is * std::cos(u_os_o), 0;

        // ECEF Satellite velocity
        gnssPosVel.sat_v_es_e(j, 0) = v_os_o(0) * std::cos(Omega) - v_os_o(1) * std::cos(inclination) * std::sin(Omega) + kOmega_ie * gnssPosVel.sat_r_es_e(j, 1);
        gnssPosVel.sat_v_es_e(j, 1) = v_os_o(0) * std::sin(Omega) + v_os_o(1) * std::cos(inclination) * std::cos(Omega) - kOmega_ie * gnssPosVel.sat_r_es_e(j, 0);
        gnssPosVel.sat_v_es_e(j, 2) = v_os_o(1) * std::sin(inclination);
    }

    return gnssPosVel;
}

template<typename T>
typename Types<T>::GnssMeasurements 
Simulation<T>::generateGnssMeasurements(const T & time,
                                        const SatPosVel & gnss_pos_vel,
                                        const NavSolutionNed& true_nav_ned,
                                        const NavSolutionEcef& true_nav_ecef,
                                        const Eigen::Matrix<T, Eigen::Dynamic, 1, 0, kMaxGnssSatellites, 1> & gnss_biases, 
                                        const GnssConfig& gnss_config,
                                        std::mt19937 & gen) {

    // nomally distributed error
    std::normal_distribution<T> randn(0.0, 1.0);

    GnssMeasurements gnss_measurements;

    // Initialize number of GNSS measurements
    gnss_measurements.no_meas = 0;

    // Calculate ECEF to NED coordinate transformation matrix
    T cos_lat = std::cos(true_nav_ned.latitude);
    T sin_lat = std::sin(true_nav_ned.latitude);
    T cos_long = std::cos(true_nav_ned.longitude);
    T sin_long = std::sin(true_nav_ned.longitude);
    Matrix3 C_e_n;
    C_e_n << -sin_lat * cos_long, -sin_lat * sin_long,  cos_lat,
             -sin_long,            cos_long,        0,
             -cos_lat * cos_long, -cos_lat * sin_long, -sin_lat;
     
    // Skew symmetric matrix of Earth rate
    Matrix3 Omega_ie = skewSymmetric(Vector3(0, 0, kOmega_ie));
       
    // Resize the GNSS measurements matrix to accommodate the maximum possible measurements
    gnss_measurements.meas = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(gnss_config.no_sat, 8);

    // Loop over satellites
    for (int j = 0; j < gnss_config.no_sat; ++j) {
        // Determine ECEF line-of-sight vector
        Vector3 delta_r = gnss_pos_vel.sat_r_es_e.row(j).transpose() - true_nav_ecef.r_eb_e;
        T approx_range = delta_r.norm();
        Vector3 u_as_e = delta_r / approx_range;
    
        // Convert line-of-sight vector to NED and determine elevation
        T elevation = -std::asin(C_e_n.row(2).dot(u_as_e));
    
        // Determine if satellite is above the masking angle
        if (elevation >= gnss_config.mask_angle * kDegToRad) {
            // Increment number of measurements
            gnss_measurements.no_meas++;
    
            // Calculate frame rotation during signal transit time
            Matrix3 C_e_I;
            C_e_I << 1, kOmega_ie * approx_range / kC, 0,
                     -kOmega_ie * approx_range / kC, 1, 0,
                     0, 0, 1;

            // Calculate range
            delta_r = C_e_I * gnss_pos_vel.sat_r_es_e.row(j).transpose() - true_nav_ecef.r_eb_e;
            T range = delta_r.norm();
        
            // Calculate range rate
            T range_rate = u_as_e.dot(C_e_I * (gnss_pos_vel.sat_v_es_e.row(j).transpose() + Omega_ie * gnss_pos_vel.sat_r_es_e.row(j).transpose()) - (true_nav_ecef.v_eb_e + Omega_ie * true_nav_ecef.r_eb_e));
    
            // Calculate pseudo-range measurement
            gnss_measurements.meas(gnss_measurements.no_meas-1, 0) = range + gnss_biases(j) + gnss_config.rx_clock_offset + gnss_config.rx_clock_drift * time + gnss_config.code_track_err_sd * randn(gen);
    
            // Calculate pseudo-range rate measurement
            gnss_measurements.meas(gnss_measurements.no_meas-1, 1) = range_rate + gnss_config.rx_clock_drift + gnss_config.rate_track_err_sd * randn(gen);
    
            // Append satellite position and velocity to output data
            gnss_measurements.meas.template block<1, 3>(gnss_measurements.no_meas-1, 2) = gnss_pos_vel.sat_r_es_e.row(j);
            gnss_measurements.meas.template block<1, 3>(gnss_measurements.no_meas-1, 5) = gnss_pos_vel.sat_v_es_e.row(j);
        }
    }

    // Resize the GNSS measurements matrix to the actual number of measurements
    // Which is <= n of satellites
    gnss_measurements.meas.conservativeResize(gnss_measurements.no_meas, 8);

    // 6. Set-up measurement noise covariance matrix assuming all measurements
    // are independent and have equal variance for a given measurement type.
    gnss_measurements.cov_mat =
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, 0, 2* kMaxGnssSatellites, 2* kMaxGnssSatellites>::Identity(2*gnss_measurements.no_meas, 2*gnss_measurements.no_meas);

    // Ranges
    gnss_measurements.cov_mat.block(0,0,gnss_measurements.no_meas,gnss_measurements.no_meas) 
                *= pow(gnss_config.pseudo_range_sd,2.0);
    // Range rates
    gnss_measurements.cov_mat.block(gnss_measurements.no_meas,gnss_measurements.no_meas,gnss_measurements.no_meas,gnss_measurements.no_meas)
                *= pow(gnss_config.range_rate_sd,2.0);

    return gnss_measurements;
}

template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1, 0, Simulation<T>::kMaxGnssSatellites, 1>
Simulation<T>::initializeGnssBiases(const NavSolutionEcef & true_nav_ecef,
                                    const NavSolutionNed & true_nav_ned,
                                    const SatPosVel & gnss_pos_vel,
                                    const GnssConfig& gnss_config,
                                    std::mt19937 & gen) { // gen not templated

    // nomally distributed error
    std::normal_distribution<T> randn(0.0, 1.0);

    Eigen::Matrix<T, Eigen::Dynamic, 1, 0, kMaxGnssSatellites, 1> gnss_biases(gnss_config.no_sat);

    // Calculate ECEF to NED coordinate transformation matrix
    T cos_lat = std::cos(true_nav_ned.latitude);
    T sin_lat = std::sin(true_nav_ned.latitude);
    T cos_long = std::cos(true_nav_ned.longitude);
    T sin_long = std::sin(true_nav_ned.longitude);
    Matrix3 C_e_n;
    C_e_n << -sin_lat * cos_long, -sin_lat * sin_long,  cos_lat,
             -sin_long,            cos_long,        0,
             -cos_lat * cos_long, -cos_lat * sin_long, -sin_lat;

    // Loop over satellites
    for (int j = 0; j < gnss_config.no_sat; ++j) {
        // Determine ECEF line-of-sight vector
        Vector3 delta_r = gnss_pos_vel.sat_r_es_e.row(j).transpose() - true_nav_ecef.r_eb_e;
        Vector3 u_as_e = delta_r / delta_r.norm();
    
        // Convert line-of-sight vector to NED and determine elevation
        T elevation = -std::asin(C_e_n.row(2).dot(u_as_e));
    
        // Limit the minimum elevation angle to the masking angle
        elevation = std::max(elevation, gnss_config.mask_angle * kDegToRad);
    
        // Calculate ionosphere and troposphere error standard deviations
        T iono_SD = gnss_config.zenith_iono_err_sd / std::sqrt(1 - 0.899 * std::cos(elevation) * std::cos(elevation));
        T trop_SD = gnss_config.zenith_trop_err_sd / std::sqrt(1 - 0.998 * std::cos(elevation) * std::cos(elevation));
    
        // Determine range bias
        gnss_biases(j) = gnss_config.sis_err_sd * randn(gen) + iono_SD * randn(gen) + trop_SD * randn(gen);
    }

    return gnss_biases;
}

template struct Simulation<double>;
template struct Simulation<float>;

};