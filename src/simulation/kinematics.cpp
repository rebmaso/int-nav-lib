#include "simulation.h"

namespace intnavlib {

template<typename T>
typename Types<T>::ImuMeasurements 
Simulation<T>::kinematicsEcef(const NavSolutionEcef & new_nav,
                                const NavSolutionEcef & old_nav) {
    
    T tor_i = new_nav.time - old_nav.time;

    // Init measurements to 0
    ImuMeasurements true_imu_meas;
    true_imu_meas.time = new_nav.time;
    true_imu_meas.f = Vector3::Zero();
    true_imu_meas.omega = Vector3::Zero();

    if (tor_i > 0.0) {

    // From (2.145) determine the Earth rotation over the update interval
    // C_Earth = C_e_i' * old_C_e_i

    T alpha_ie = kOmega_ie * tor_i;

    T cos_alpha_ie = cos(alpha_ie);
    T sin_alpha_ie = sin(alpha_ie);

    Matrix3 C_Earth;
    C_Earth << cos_alpha_ie, sin_alpha_ie, 0.0,
                -sin_alpha_ie, cos_alpha_ie, 0.0,
                    0.0,             0.0,  1.0;

    // Obtain coordinate transformation matrix from the old attitude (w.r.t.
    // an inertial frame) to the new (compensate for earth rotation)
    Matrix3 C_old_new = new_nav.C_b_e.transpose() * C_Earth * old_nav.C_b_e;

    // Calculate the approximate angular rate w.r.t. an inertial frame
    Vector3 alpha_ib_b;
    alpha_ib_b(0) = 0.5 * (C_old_new(1,2) - C_old_new(2,1));
    alpha_ib_b(1) = 0.5 * (C_old_new(2,0) - C_old_new(0,2));
    alpha_ib_b(2) = 0.5 * (C_old_new(0,1) - C_old_new(1,0));

    // Calculate and apply the scaling factor
    T temp = acos(0.5 * (C_old_new(0,0) + C_old_new(1,1) + C_old_new(2,2) - 1.0));
    if (temp > 2.0e-5) // scaling is 1 if temp is less than this
        alpha_ib_b = alpha_ib_b * temp/sin(temp);
    
    // Calculate the angular rate
    true_imu_meas.omega = alpha_ib_b / tor_i;

    // Calculate the specific force resolved about ECEF-frame axes
    // From (5.36)
    Vector3 kOmega_ie_vec;
    kOmega_ie_vec << 0.0 , 0.0 , kOmega_ie;
    Vector3 f_ib_e = ((new_nav.v_eb_e - old_nav.v_eb_e) / tor_i) - gravityEcef(old_nav.r_eb_e)
        + 2.0 * skewSymmetric(kOmega_ie_vec) * old_nav.v_eb_e;

    // Calculate the average body-to-ECEF-frame coordinate transformation
    // matrix over the update interval using (5.84) and (5.85)

    T mag_alpha = alpha_ib_b.norm();
    Matrix3 Alpha_ib_b = skewSymmetric(alpha_ib_b);

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
    
    // Transform specific force to body-frame resolving axes using (5.81)
    // plain inverse as matrix is 3x3
    // true_imu_meas.f = ave_C_b_e.inverse() * f_ib_e;
    true_imu_meas.f = ave_C_b_e.colPivHouseholderQr().solve(f_ib_e);

    }

    return true_imu_meas;

}


template<typename T>
typename Types<T>::ImuMeasurements 
Simulation<T>::imuModel(const ImuMeasurements & true_imu_meas,
                        const ImuMeasurements & old_imu_meas,
                        const ImuErrors & imu_errors,
                        const T & tor_i,
                        std::mt19937 & gen) {
    
    ImuMeasurements imu_measurements;

    imu_measurements.time = true_imu_meas.time;

    // Init noises
    Vector3 accel_noise = Vector3::Zero();
    Vector3 gyro_noise = Vector3::Zero();

    // Normal distribution init
    std::normal_distribution<T> randn(0.0, 1.0); 

    if(tor_i > 0.0) {
        accel_noise << randn(gen) * imu_errors.accel_noise_root_psd / sqrt(tor_i), 
                        randn(gen) * imu_errors.accel_noise_root_psd / sqrt(tor_i),
                        randn(gen) * imu_errors.accel_noise_root_psd / sqrt(tor_i);  
        gyro_noise << randn(gen) * imu_errors.gyro_noise_root_psd / sqrt(tor_i),
                        randn(gen) * imu_errors.gyro_noise_root_psd / sqrt(tor_i),
                        randn(gen) * imu_errors.gyro_noise_root_psd / sqrt(tor_i);
    }

    // Calculate accelerometer and gyro outputs using (4.16) and (4.17)

    // Specific force
    Vector3 uq_f_ib_b = imu_errors.b_a + (Matrix3::Identity() + imu_errors.M_a) * true_imu_meas.f 
                                + accel_noise;

    // Angular velocity
    Vector3 uq_omega_ib_b = imu_errors.b_g + (Matrix3::Identity() + imu_errors.M_g) * true_imu_meas.omega 
                                    + imu_errors.G_g * true_imu_meas.f 
                                    + gyro_noise;

    // Quantize accelerometer outputs
    if (imu_errors.accel_quant_level > kEpsilon) {
        Vector3 temp = (uq_f_ib_b + old_imu_meas.quant_residuals_f) / imu_errors.accel_quant_level;
        for(size_t i = 0; i < 3; i++) temp(i) = round(temp(i));
        imu_measurements.f = imu_errors.accel_quant_level * temp;
        imu_measurements.quant_residuals_f = uq_f_ib_b + old_imu_meas.quant_residuals_f -
            imu_measurements.f;
    }
    else{
        imu_measurements.f = uq_f_ib_b;
        imu_measurements.quant_residuals_f = Vector3::Zero();
    }

    // Quantize gyro outputs
    if (imu_errors.gyro_quant_level > kEpsilon) {
        Vector3 temp = (uq_omega_ib_b + old_imu_meas.quant_residuals_omega) / imu_errors.gyro_quant_level;
        for(size_t i = 0; i < 3; i++) temp(i) = round(temp(i)); 
        imu_measurements.omega = imu_errors.gyro_quant_level * temp;
        imu_measurements.quant_residuals_omega = uq_omega_ib_b + old_imu_meas.quant_residuals_omega -
            imu_measurements.omega;
    }
    else {
        imu_measurements.omega = uq_omega_ib_b;
        imu_measurements.quant_residuals_omega = Vector3::Zero();
    }

    // // No quantization
    // imu_measurements.f = uq_f_ib_b;
    // imu_measurements.omega = uq_omega_ib_b;

    return imu_measurements;
}

template<typename T>
typename Types<T>::PosMeasEcef 
Simulation<T>::genericPosSensModel(const NavSolutionEcef & true_nav, 
                                const T & pos_sigma,
                                std::mt19937 & gen){

    PosMeasEcef pos_meas;
    pos_meas.time = true_nav.time;

    // nomally distributed error
    std::normal_distribution<T> randn(0.0, pos_sigma);

    Vector3 pos_error;
    pos_error << randn(gen), randn(gen), randn(gen);

    Matrix3 cov_mat = Matrix3::Identity() * pos_sigma * pos_sigma;

    // add error
    pos_meas.r_eb_e = true_nav.r_eb_e + pos_error;
    pos_meas.cov_mat = cov_mat;

    return pos_meas;
}

template<typename T>
typename Types<T>::PosRotMeasEcef 
Simulation<T>::genericPosRotSensModel(const NavSolutionEcef & true_nav, 
                                        const T & pos_sigma, // This was correct
                                        const T & rot_sigma,
                                        std::mt19937 & gen){ // gen not templated
    PosRotMeasEcef pos_rot_meas;
    pos_rot_meas.time = true_nav.time;

    // nomally distributed error
    std::normal_distribution<T> randn_pos(0.0, pos_sigma);
    std::normal_distribution<T> randn_rot(0.0, rot_sigma);

    Vector3 pos_error;
    pos_error << randn_pos(gen), randn_pos(gen), randn_pos(gen);

    Vector3 rot_error;
    rot_error << randn_rot(gen), randn_rot(gen), randn_rot(gen);
    Matrix3 C_b_b = eulerToDcm(rot_error); // perturbation Rot mat

    Eigen::Matrix<T,6,6> cov_mat = Eigen::Matrix<T,6,6>::Identity();
    cov_mat.template block<3,3>(0,0) *= pos_sigma * pos_sigma;
    cov_mat.template block<3,3>(3,3) *= rot_sigma * rot_sigma;

    // add error
    pos_rot_meas.r_eb_e = true_nav.r_eb_e + pos_error;
    pos_rot_meas.C_b_e = true_nav.C_b_e * C_b_b;
    pos_rot_meas.cov_mat = cov_mat;

    return pos_rot_meas;
}

template struct Simulation<double>;
template struct Simulation<float>;

};