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

#ifdef BUILD_VISION

template<typename T>
Simulation<T>::ViewerInterface::ViewerInterface(const std::string& earth_file_path,
                                                double fx, 
                                                double fy, 
                                                double cx, 
                                                double cy, 
                                                unsigned int width, 
                                                unsigned int height, 
                                                double near, 
                                                double far, 
                                                unsigned int delay_frames,
                                                unsigned int init_delay_frames) :

                                                fx(fx),
                                                fy(fy),
                                                cx(cx),
                                                cy(cy),
                                                width(width),
                                                height(height),
                                                near(near),
                                                far(far),
                                                delay_frames(delay_frames),
                                                init_delay_frames(init_delay_frames) {

        osgEarth::initialize();

        std::vector<std::string> fake_argv_str;
        fake_argv_str.push_back("nav_sim_ol"); // dummy program name
        fake_argv_str.push_back("--nologdepth");
        fake_argv_str.push_back(earth_file_path);

        std::vector<char*> fake_argv;
        for (auto& s : fake_argv_str) {
            fake_argv.push_back(&s[0]);
        }
        int fake_argc = fake_argv.size();
        arguments = std::make_unique<osg::ArgumentParser>(&fake_argc, fake_argv.data());
        
        viewer = new osgViewer::Viewer(*arguments);
        
        // // EGL graphics context for headless rendering
        // // When headless, only works on NVIDIA GPUs it seems
        // // If i call eglinitialize (eg with eglinfo) it fails when headless
        // // TODO: try with the exact same impl from the nvidia example
        // // https://community.khronos.org/t/offscreen-rendering-using-opengl-and-egl/109487/19
        // // https://developer.nvidia.com/blog/egl-eye-opengl-visualization-without-x-server/
        // // https://github.com/omaralvarez/osgEGL/tree/master
        gc = new EGLGraphicsWindowEmbedded(static_cast<int>(width), static_cast<int>(height));
        viewer->getCamera()->setGraphicsContext(gc);
        viewer->getCamera()->setViewport(new osg::Viewport(0, 0, width, height));
        viewer->realize();

        // Setup hidden viewer using pbuffer
        // Should work without a display on intel GPUs
        // https://osg-users.openscenegraph.narkive.com/pIDGiyl0/hidden-viewer
        // traits = new osg::GraphicsContext::Traits;
        // traits->x = 0;
        // traits->y = 0;
        // traits->width = width;
        // traits->height = height;
        // traits->windowDecoration = false;
        // traits->doubleBuffer = false;
        // traits->sharedContext = 0;
        // traits->pbuffer = true; // Set pixel buffer
        // traits->readDISPLAY();
        // gc = osg::GraphicsContext::createGraphicsContext(traits.get());
        // gc->realize();
        // gc->makeCurrent();
        // viewer->getCamera()->setGraphicsContext(gc);
        // viewer->getCamera()->setViewport(new osg::Viewport(0, 0, width, height));

        manip = new SimpleManipulator();
        viewer->setCameraManipulator(manip.get());

        viewer->getCamera()->setSmallFeatureCullingPixelSize(-1.0f);
        viewer->getCamera()->setProjectionResizePolicy(osg::Camera::FIXED); 
        viewer->getCamera()->setComputeNearFarMode(osg::CullSettings::DO_NOT_COMPUTE_NEAR_FAR);

        auto map_node = MapNodeHelper().load(*arguments, viewer.get());
        if (map_node.valid() && MapNode::get(map_node))
            viewer->setSceneData(map_node);

        setProj(fx, fy, cx, cy, width, height, near, far);

        // LOG(INFO) << "" << std::endl;
        // LOG(INFO) << "" << "OpenGL Vendor: " << glGetString(GL_VENDOR) << std::endl;
        // LOG(INFO) << "" << "OpenGL Renderer: " << glGetString(GL_RENDERER) << std::endl;
        // LOG(INFO) << "" << "OpenGL Version: " << glGetString(GL_VERSION) << std::endl;
        // LOG(INFO) << "" << "OpenGL Shading Language Version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << std::endl;
        // LOG(INFO) << "" << std::endl;

        // Allocate and attach a persistent color image.
        color_image = new osg::Image();
        color_image->allocateImage(static_cast<int>(width), static_cast<int>(height), 1, GL_RGB, GL_UNSIGNED_BYTE);
        viewer->getCamera()->attach(osg::Camera::COLOR_BUFFER, color_image.get());

        // Allocate and attach a persistent depth image.
        depth_image = new osg::Image();
        depth_image->allocateImage(static_cast<int>(width), static_cast<int>(height), 1, GL_DEPTH_COMPONENT, GL_FLOAT);
        viewer->getCamera()->attach(osg::Camera::DEPTH_BUFFER, depth_image.get());

}

template<typename T>
void Simulation<T>::ViewerInterface::setPose(const Eigen::Matrix3d &rotation, const Eigen::Vector3d &translation) {

	// Construct 4x4 Eigen transformation matrix
	Eigen::Matrix4d eigenViewMatrix = Eigen::Matrix4d::Identity();

	// Invert here
	eigenViewMatrix.block<3, 3>(0, 0) = rotation.transpose();
	eigenViewMatrix.block<3, 1>(0, 3) = - rotation.transpose() * translation;


	// Flip Y and Z axes
	eigenViewMatrix.row(1) *= -1;
	eigenViewMatrix.row(2) *= -1;

	// Transpose the matrix
	eigenViewMatrix.transposeInPlace();

	// Construct osg::Matrixd 
	osg::Matrixd viewMatrix(
	    eigenViewMatrix(0, 0), eigenViewMatrix(0, 1), eigenViewMatrix(0, 2), eigenViewMatrix(0, 3),
	    eigenViewMatrix(1, 0), eigenViewMatrix(1, 1), eigenViewMatrix(1, 2), eigenViewMatrix(1, 3),
	    eigenViewMatrix(2, 0), eigenViewMatrix(2, 1), eigenViewMatrix(2, 2), eigenViewMatrix(2, 3),
	    eigenViewMatrix(3, 0), eigenViewMatrix(3, 1), eigenViewMatrix(3, 2), eigenViewMatrix(3, 3));
	    
	// Apply the matrix
    	manip->setByInverseMatrix(viewMatrix);
};

template<typename T>
void Simulation<T>::ViewerInterface::setProj(const double & fx, 
	const double & fy, 
	const double & cx, 
	const double & cy, 
	const double & w, 
	const double & h, 
	const double & near, 
	const double & far) {
	osg::Matrixd projMatrix(2 * fx / w, 0.0, (w - 2 * cx) / w, 0.0,
	    0.0, 2 * fy / h, (h - 2 * cy) / h, 0.0,
	    0.0, 0.0, (-far - near) / (far - near), -2.0 * far * near / (far - near),
	    0.0, 0.0, -1.0, 0.0);
	projMatrix.transpose(projMatrix);

	viewer->getCamera()->setProjectionMatrix(projMatrix);
}

template<typename T>
void Simulation<T>::ViewerInterface::startWithInit(const Eigen::Matrix3d &rotation, 
                                        const Eigen::Vector3d &translation) {

	setPose(rotation, translation);
	// LOG(INFO) << "Scene Loading...";
	for(int i = 0; i < init_delay_frames; i++) viewer->frame();
	// LOG(INFO) << "Ready!";                      
	renderThread = std::thread(&ViewerInterface::renderLoop, this);
}

template<typename T>
cv::Mat Simulation<T>::ViewerInterface::getSnapShot() {

	// Wait some frames - producer-consumer model
	std::unique_lock<std::mutex> lk(viewer_mutex);
	viewer_flag = true;
	viewer_cv.wait(lk,[&]{return !viewer_flag;});

	// Read the attached color image, flip it vertically, and encode.
	cv::Mat raw_image;
	raw_image = cv::Mat(color_image->t(), color_image->s(), CV_8UC3, color_image->data()).clone();

	cv::Mat final_image;
	cv::flip(raw_image, final_image, 0);
	cv::cvtColor(final_image, final_image, cv::COLOR_BGR2GRAY);

	return final_image.clone();
}

template<typename T>
void Simulation<T>::ViewerInterface::renderLoop () {
	// Wait for first pose to be set, before looping
	while(!viewer->done()) {
	    viewer->frame();
	    if(viewer_flag) {
		for(size_t i = 0; i<delay_frames; i++) viewer->frame();
		viewer_flag = false;
		viewer_cv.notify_one();
	    }
	}
}

#endif

template struct Simulation<double>;
template struct Simulation<float>;

};