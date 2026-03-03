#include "helpers.h"

#ifdef BUILD_VISION

#include <mutex>
#include <opencv2/calib3d.hpp>
#include <gdal_priv.h>
#include <cpl_conv.h>
#include <ogr_spatialref.h>
#include <ogrsf_frmts.h>

#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#endif


namespace intnavlib {

template <typename T>
typename Types<T>::NavSolutionEcef 
Helpers<T>::nedToEcef(const NavSolutionNed& nav_sol_ned) {
    // Extract the inputs
    T L_b = nav_sol_ned.latitude;
    T lambda_b = nav_sol_ned.longitude;
    T h_b = nav_sol_ned.height;
    Vector3 v_eb_n = nav_sol_ned.v_eb_n;
    Matrix3 C_b_n = nav_sol_ned.C_b_n;
    T sin_lat = std::sin(L_b);
    T cos_lat = std::cos(L_b);
    T sin_long = std::sin(lambda_b);
    T cos_long = std::cos(lambda_b);
    // Compute transverse radius of curvature
    T R_E = kR0 / std::sqrt(1.0 - pow((kEccentricity * sin_lat),2.0));
    // Convert position from curvilinear to Cartesian ECEF coordinates
    Vector3 r_eb_e;
    r_eb_e(0) = (R_E + h_b) * cos_lat * cos_long;
    r_eb_e(1) = (R_E + h_b) * cos_lat * sin_long;
    r_eb_e(2) = ((1.0 - kEccentricity * kEccentricity) * R_E + h_b) * sin_lat;
    // Compute the ECEF to NED coordinate transformation matrix
    Matrix3 C_e_n;
    C_e_n << -sin_lat * cos_long, -sin_lat * sin_long, cos_lat,
             -sin_long,            cos_long,          0.0,
             -cos_lat * cos_long, -cos_lat * sin_long, -sin_lat;
    // Transform velocity from NED to ECEF frame
    Vector3 v_eb_e = C_e_n.transpose() * v_eb_n;
    // Transform attitude from NED to ECEF frame
    Matrix3 C_b_e = C_e_n.transpose() * C_b_n;
    // Construct the output structure
    NavSolutionEcef nav_sol_ecef;
    nav_sol_ecef.time = nav_sol_ned.time;
    nav_sol_ecef.r_eb_e = r_eb_e;
    nav_sol_ecef.v_eb_e = v_eb_e;
    nav_sol_ecef.C_b_e = C_b_e;
    return nav_sol_ecef;
}

template <typename T>
typename Types<T>::NavSolutionNed 
Helpers<T>::ecefToNed(const NavSolutionEcef & nav_sol_ecef){
    // Convert position using Borkowski closed-form exact solution
    // From (2.113)
    T lambda_b = atan2(nav_sol_ecef.r_eb_e(1), nav_sol_ecef.r_eb_e(0));
    // From (C.29) and (C.30)
    T k1 = sqrt(1.0 - pow(kEccentricity,2.0)) * abs(nav_sol_ecef.r_eb_e(2));
    T k2 = pow(kEccentricity,2.0) * kR0;
    T beta = sqrt(pow(nav_sol_ecef.r_eb_e(0),2.0) + pow(nav_sol_ecef.r_eb_e(1),2.0));
    T E = (k1 - k2) / beta;
    T F = (k1 + k2) / beta;
    // From (C.31)
    T P = 4.0/3.0 * (E*F + 1.0);
    // From (C.32)
    T Q = 2.0 * (pow(E,2.0) - pow(F,2.0));
    // From (C.33)
    T D = pow(P,3.0) + pow(Q,2.0);
    // From (C.34)
    T V = pow(sqrt(D) - Q, 1.0/3.0) - pow(sqrt(D) + Q,1.0/3.0);
    // From (C.35)
    T G = 0.5 * (sqrt(pow(E,2.0) + V) + E);
    // From (C.36)
    T T_ = sqrt(pow(G,2.0) + (F - V * G) / (2.0 * G - E)) - G;
    // From (C.37)
    T L_b = sgn(nav_sol_ecef.r_eb_e(2)) * atan((1.0 - pow(T_,2.0)) / (2.0 * T_ * sqrt (1.0 - pow(kEccentricity,2.0))));
    // From (C.38)
    T h_b = (beta - kR0 * T_) * cos(L_b) +
        (nav_sol_ecef.r_eb_e(2) - sgn(nav_sol_ecef.r_eb_e(2)) * kR0 * sqrt(1.0 - pow(kEccentricity,2.0))) * sin(L_b);  
    // Calculate ECEF to NED coordinate transformation matrix using (2.150)
    T cos_lat = cos(L_b);
    T sin_lat = sin(L_b);
    T cos_long = cos(lambda_b);
    T sin_long = sin(lambda_b);
    Matrix3 C_e_n;
    C_e_n << -sin_lat * cos_long, -sin_lat * sin_long,  cos_lat,
                    -sin_long,            cos_long,        0.0,
            -cos_lat * cos_long, -cos_lat * sin_long, -sin_lat;
    // Transform velocity using (2.73)
    Vector3 v_eb_n = C_e_n * nav_sol_ecef.v_eb_e;
    // Transform attitude using (2.15)
    Matrix3 C_b_n = C_e_n * nav_sol_ecef.C_b_e;
    NavSolutionNed nav_sol_ned;
    nav_sol_ned.time = nav_sol_ecef.time;
    nav_sol_ned.latitude = L_b;
    nav_sol_ned.longitude = lambda_b;
    nav_sol_ned.height = h_b;
    nav_sol_ned.v_eb_n = v_eb_n;
    nav_sol_ned.C_b_n = C_b_n;
    return nav_sol_ned;
}

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
typename Helpers<T>::Matrix3 
Helpers<T>::eulerToDcm(const typename Helpers<T>::Vector3 & rpy) {
    T roll = rpy(0);
    T pitch = rpy(1);
    T yaw = rpy(2);
    // Precompute sines and cosines of Euler angles
    T sin_roll = std::sin(roll);
    T cos_roll = std::cos(roll);
    T sin_pitch = std::sin(pitch);
    T cos_pitch = std::cos(pitch);
    T sin_yaw = std::sin(yaw);
    T cos_yaw = std::cos(yaw);
    // Calculate the coordinate transformation matrix R = Rz(yaw) * Ry(pitch) * Rx(roll)
    Matrix3 C;
    C(0,0) = cos_pitch * cos_yaw;
    C(1,0) = cos_pitch * sin_yaw;
    C(2,0) = -sin_pitch;
    C(0,1) = -cos_roll * sin_yaw + sin_roll * sin_pitch * cos_yaw;
    C(1,1) = cos_roll * cos_yaw + sin_roll * sin_pitch * sin_yaw;
    C(2,1) = sin_roll * cos_pitch;
    C(0,2) = sin_roll * sin_yaw + cos_roll * sin_pitch * cos_yaw;
    C(1,2) = -sin_roll * cos_yaw + cos_roll * sin_pitch * sin_yaw;
    C(2,2) = cos_roll * cos_pitch;
    return C;
}

template <typename T>
typename Helpers<T>::Vector3 
Helpers<T>::dcmToEuler(const typename Helpers<T>::Matrix3 & C) {
    Vector3 rpy;
    rpy(0) = atan2(C(2,1),C(2,2));
    rpy(1) = - asin(C(2,0));      
    rpy(2) = atan2(C(1,0),C(0,0));
    return rpy;
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

#ifdef BUILD_VISION
template <typename T>
cv::Mat Helpers<T>::getGdalReferenceImage(const typename Types<T>::NavSolutionEcef& nav_sol,
                                          T fx, T fy, T cx, T cy,
                                          unsigned int width, unsigned int height,
                                          const typename Helpers<T>::Matrix3& C_b_c,
                                          const std::string& dataset_path,
                                          const std::string& dem_dataset_path) {
    
    // 1. GDAL setup
    static std::once_flag gdal_initialized;
    std::call_once(gdal_initialized, [](){ GDALAllRegister(); });

    // 2. Get camera pose in NED
    NavSolutionNed nav_sol_ned = ecefToNed(nav_sol);
    Matrix3 C_c_n = nav_sol_ned.C_b_n * C_b_c;

    // 3. Get ground height from DEM
    double ground_height = 0.0;
    GDALDataset* dem_dataset = (GDALDataset*) GDALOpen(dem_dataset_path.c_str(), GA_ReadOnly);
    if (dem_dataset) {
        OGRSpatialReference oSRS, iSRS;
        const char* wkt_projection = dem_dataset->GetProjectionRef();
        oSRS.importFromWkt(wkt_projection);
        iSRS.SetWellKnownGeogCS("WGS84");
        iSRS.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
        oSRS.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

        OGRCoordinateTransformation *poCT = OGRCreateCoordinateTransformation(&iSRS, &oSRS);
        if (poCT) {
            double x = nav_sol_ned.longitude * kRadToDeg;
            double y = nav_sol_ned.latitude * kRadToDeg;
            if(poCT->Transform(1, &x, &y)) {
                double adfGeoTransform[6];
                dem_dataset->GetGeoTransform(adfGeoTransform);
                double invGeoTransform[6];
                if (GDALInvGeoTransform(adfGeoTransform, invGeoTransform)) {
                    double pixel, line;
                    GDALApplyGeoTransform(invGeoTransform, x, y, &pixel, &line);
                    int px = static_cast<int>(std::floor(pixel));
                    int py = static_cast<int>(std::floor(line));
                    
                    if (px >= 0 && px < dem_dataset->GetRasterXSize() && py >= 0 && py < dem_dataset->GetRasterYSize()) {
                        float val;
                        if(dem_dataset->GetRasterBand(1)->RasterIO(GF_Read, px, py, 1, 1, &val, 1, 1, GDT_Float32, 0, 0) == CE_None) {
                            ground_height = static_cast<double>(val);
                        }
                    }
                }
            }
            OCTDestroyCoordinateTransformation(poCT);
        }
        GDALClose(dem_dataset);
    } else {
        std::cerr << "Failed to open DEM dataset: " << dem_dataset_path << std::endl;
    }

    GDALDataset* dataset = (GDALDataset*) GDALOpen(dataset_path.c_str(), GA_ReadOnly);
    if (!dataset) {
        std::cerr << "Failed to open GDAL dataset: " << dataset_path << std::endl;
        return cv::Mat();
    }

    // 4. Define FOV corner rays in camera frame
    std::vector<Vector3> corner_rays_c;
    corner_rays_c.push_back(Vector3((0 - cx) / fx, (0 - cy) / fy, 1.0).normalized());
    corner_rays_c.push_back(Vector3((width - cx) / fx, (0 - cy) / fy, 1.0).normalized());
    corner_rays_c.push_back(Vector3((width - cx) / fx, (height - cy) / fy, 1.0).normalized());
    corner_rays_c.push_back(Vector3((0 - cx) / fx, (height - cy) / fy, 1.0).normalized());

    // 5. Intersect rays with ground plane in NED
    std::vector<cv::Point2d> ground_corners_lon_lat;
    Vector2 radii = radiiOfCurvature(nav_sol_ned.latitude);
    T R_N = radii(0);
    T R_E = radii(1);

    for (const auto& ray_c : corner_rays_c) {
        Vector3 ray_n = C_c_n * ray_c;
        if (ray_n(2) <= 0) { // Ray is not pointing down
            GDALClose(dataset);
            // std::cerr << "Camera is not pointing towards the ground." << std::endl;
            return cv::Mat();
        }
        T t = (nav_sol_ned.height - ground_height) / ray_n(2);
        Vector3 p_intersect_n = t * ray_n;

        // 6. Convert intersection points to Lat/Lon
        T corner_lat = nav_sol_ned.latitude + p_intersect_n(0) / (R_N + nav_sol_ned.height);
        T corner_lon = nav_sol_ned.longitude + p_intersect_n(1) / ((R_E + nav_sol_ned.height) * cos(nav_sol_ned.latitude));
        ground_corners_lon_lat.push_back(cv::Point2d(corner_lon * kRadToDeg, corner_lat * kRadToDeg));
    }

    // 7. Transform corner coordinates to dataset's pixel coordinates
    OGRSpatialReference oSRS, iSRS;
    const char* wkt_projection = dataset->GetProjectionRef();
    oSRS.importFromWkt(wkt_projection);
    iSRS.SetWellKnownGeogCS("WGS84");
    iSRS.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    oSRS.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

    OGRCoordinateTransformation *poCT = OGRCreateCoordinateTransformation(&iSRS, &oSRS);
    if (!poCT) {
        GDALClose(dataset);
        std::cerr << "Failed to create coordinate transformation." << std::endl;
        return cv::Mat();
    }

    double adfGeoTransform[6];
    dataset->GetGeoTransform(adfGeoTransform);
    double invGeoTransform[6];
    if (GDALInvGeoTransform(adfGeoTransform, invGeoTransform) == 0) {
        GDALClose(dataset);
        OCTDestroyCoordinateTransformation(poCT);
        std::cerr << "Failed to compute inverse geotransform." << std::endl;
        return cv::Mat();
    }

    std::vector<cv::Point2f> src_pts;
    double min_x = std::numeric_limits<double>::max(), max_x = std::numeric_limits<double>::min();
    double min_y = std::numeric_limits<double>::max(), max_y = std::numeric_limits<double>::min();

    for (const auto& pt_lon_lat : ground_corners_lon_lat) {
        double x = pt_lon_lat.x;
        double y = pt_lon_lat.y;
        poCT->Transform(1, &x, &y);

        double pixel, line;
        GDALApplyGeoTransform(invGeoTransform, x, y, &pixel, &line);
        src_pts.push_back(cv::Point2f(pixel, line));

        min_x = std::min(min_x, pixel);
        max_x = std::max(max_x, pixel);
        min_y = std::min(min_y, line);
        max_y = std::max(max_y, line);
    }
    OCTDestroyCoordinateTransformation(poCT);

    // 8. Read data from dataset (Crop)
    int xoff = std::max(0, (int)floor(min_x));
    int yoff = std::max(0, (int)floor(min_y));
    int req_width = (int)ceil(max_x) - xoff;
    int req_height = (int)ceil(max_y) - yoff;

    if (xoff + req_width > dataset->GetRasterXSize()) req_width = dataset->GetRasterXSize() - xoff;
    if (yoff + req_height > dataset->GetRasterYSize()) req_height = dataset->GetRasterYSize() - yoff;

    if (req_width <= 0 || req_height <= 0) {
        GDALClose(dataset);
        return cv::Mat();
    }

    cv::Mat gdal_chip(req_height, req_width, CV_8UC1);
    CPLErr err = dataset->GetRasterBand(1)->RasterIO(GF_Read, xoff, yoff, req_width, req_height, gdal_chip.data, req_width, req_height, GDT_Byte, 0, 0);
    GDALClose(dataset);
    if (err != CE_None) { std::cerr << "RasterIO failed." << std::endl; return cv::Mat(); }

    return gdal_chip;
}
#endif

template struct Helpers<float>;
template struct Helpers<double>;

};