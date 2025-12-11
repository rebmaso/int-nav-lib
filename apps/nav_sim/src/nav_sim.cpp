
#include <iostream>
#include <random>
#include <filesystem>

#include <glog/logging.h>
#include <opencv2/opencv.hpp>

#include "intnavlib/intnavlib.h"

#include <FGFDMExec.h>

/// @example nav_sim_cl.cpp
/// Integrated navigation demo script - both open and closed loop with JSBSim

using namespace intnavlib;

// Typedefs for convenience

using ScalarType = double;

using Vector3 = Eigen::Matrix<ScalarType,3,1>;
using Vector2 = Eigen::Matrix<ScalarType,2,1>;
using Matrix3 = Eigen::Matrix<ScalarType,3,3>;

using NavSolutionNed = Types<ScalarType>::NavSolutionNed;
using ImuErrors = Types<ScalarType>::ImuErrors;
using GnssConfig = Types<ScalarType>::GnssConfig;
using KfConfig = Types<ScalarType>::KfConfig;
using FileWriter = Helpers<ScalarType>::FileWriter;
using NavSolutionEcef = Types<ScalarType>::NavSolutionEcef;
using ImuMeasurements = Types<ScalarType>::ImuMeasurements;
using SatPosVel = Types<ScalarType>::SatPosVel;
using GnssMeasurements = Types<ScalarType>::GnssMeasurements;
using NavKF = Navigation<ScalarType>::NavKF;
using PosMeasEcef = Types<ScalarType>::PosMeasEcef;
using StateEstEcef = Types<ScalarType>::StateEstEcef;
using EvalDataEcef = Types<ScalarType>::EvalDataEcef;
using MotionProfileReader = Helpers<ScalarType>::MotionProfileReader;
using PosRotMeasEcef = Types<ScalarType>::PosRotMeasEcef;

static constexpr ScalarType kFeetToMeters = 0.3048;

enum SimType {
    INS,
    INS_POS,
    INS_POS_ROT,
    INS_GNSS_LC,
    INS_GNSS_TC,
    INS_OF,
    UNKNOWN
};

SimType parseSimType(const std::string& sim_type) {
    if (sim_type == "ins")                  return SimType::INS;
    if (sim_type == "ins_pos")              return SimType::INS_POS;
    if (sim_type == "ins_pos_rot")          return SimType::INS_POS_ROT;
    if (sim_type == "ins_gnss_lc")          return SimType::INS_GNSS_LC;
    if (sim_type == "ins_gnss_tc")          return SimType::INS_GNSS_TC;
    if (sim_type == "ins_of")               return SimType::INS_OF;
    return SimType::UNKNOWN;
}

// Helper to wrap angle between -pi and pi
inline double wrapPi(double angle) {
    return std::atan2(std::sin(angle), std::cos(angle));
}

// Helper to get na sol from FDM
NavSolutionNed getNavSolFromFdm(const JSBSim::FGFDMExec & fdm) {

    std::shared_ptr<JSBSim::FGPropagate> prop = fdm.GetPropagate();

    NavSolutionNed true_nav_ned;

    true_nav_ned.time = fdm.GetSimTime();
    true_nav_ned.latitude = prop->GetGeodLatitudeRad();
    true_nav_ned.longitude = prop->GetLongitude();

    // Height: convert to m
    true_nav_ned.height = prop->GetGeodeticAltitude() * kFeetToMeters;

    // Velocity: convert to m/s
    true_nav_ned.v_eb_n = Vector3(prop->GetVel(1), prop->GetVel(2), prop->GetVel(3)) * kFeetToMeters;

    auto C = prop->GetTb2l();
    true_nav_ned.C_b_n << C(1,1), C(1,2), C(1,3),
                            C(2,1), C(2,2), C(2,3),
                            C(3,1), C(3,2), C(3,3);
    
    return true_nav_ned;
}

int main(int argc, char* argv[]) {

    google::InitGoogleLogging(argv[0]);
    FLAGS_stderrthreshold = 0; // INFO: 0, WARNING: 1, ERROR: 2, FATAL: 3
    FLAGS_colorlogtostderr = 1;

    if (argc != 2) {
        LOG(ERROR) << "Usage: " << argv[0] << " <config.json>";
        LOG(ERROR) << "Example: ./build/" << argv[0] << " config/closed_loop.json";
        return 1;
    }

    // ============ Parse JSON config file =============

    cv::FileStorage fs(argv[1], cv::FileStorage::READ | cv::FileStorage::FORMAT_JSON);

    if (!fs.isOpened()) {
        LOG(ERROR) << "Failed to open file.";
        return 1;
    }

    const std::string sim_type_str = fs["sim_type"];
    const bool closed_loop = bool(int(fs["closed_loop"]));
    double fdm_end_time = fs["fdm_end_time"]; // TODO Get from xml
    std::string fdm_script_path = fs["fdm_script_path"];
    std::string motion_profile_path = fs["motion_profile_path"];
    ScalarType epoch_interval = fs["epoch_interval"];

    // Print config
    LOG(INFO) << "Simulation type: " << sim_type_str;
    LOG(INFO) << "Closed loop: " << closed_loop;    
    LOG(INFO) << "FDM script path: " << fdm_script_path;
    LOG(INFO) << "FDM end time: " << fdm_end_time;
    LOG(INFO) << "Motion profile path: " << motion_profile_path;

    // ==============================================

    // Get simulation type
    SimType sim_type = parseSimType(sim_type_str);
    if(sim_type == SimType::UNKNOWN) {
        LOG(ERROR) << "Unknown simulation type";
        return 1;
    }

    // Load JSB Script / Init motion profile reader
    JSBSim::FGFDMExec fdm;
    MotionProfileReader reader;
    if(closed_loop) {
        // Load fdm script
        if (!fdm.LoadScript(SGPath(fdm_script_path))) {
            LOG(ERROR) << "Failed to load JSBSim script: " << fdm_script_path;
            return 1;
        }
        // Initialize the fdm
        fdm.RunIC();
    }
    else {
        // Read motion profile
        reader = MotionProfileReader(motion_profile_path);
    }

    // ============== Init sim ==============

    // Camera parameters
    Matrix3 K;
    K << 400.0, 0.0, 200.0,
         0.0, 400.0, 200.0,
         0.0, 0.0, 1.0;
    Matrix3 C_b_c = Matrix3::Identity();

    // OSG viewer
    Simulationd::ViewerInterface viewer_interface("config/faster_world.earth", 
                                                K(0,0), 
                                                K(1,1), 
                                                K(0,2), 
                                                K(1,2),
                                                400, 
                                                400, 
                                                1.0, 
                                                20000.0, 
                                                10, 
                                                1);
    
    // Start render loop
    viewer_interface.start();

    std::random_device rd;
    // Mersenne twister PRNG, initialized with random seed
    std::mt19937 gen(rd());

    // Tactcal grade IMU errors
    ImuErrors imu_errors = Helpers<ScalarType>::tacticalImuErrors();

    // Default GNSS config
    GnssConfig gnss_config = Helpers<ScalarType>::defaultGnssConfig();

    // Tactcal grade IMU - KF config
    KfConfig kf_config = Helpers<ScalarType>::tacticalImuKFConfig();

    // Init profile writer
    std::string new_directory = "./results";
    std::string base_filename;
    if(closed_loop)
        base_filename = std::filesystem::path(fdm_script_path).filename().string();
    else
        base_filename = std::filesystem::path(motion_profile_path).filename().string();
    std::string filename_without_extension = base_filename.substr(0, base_filename.find_last_of('.'));
    std::string eval_filename_out = new_directory + "/" + filename_without_extension + "_eval.csv";
    std::string profile_filename_out = new_directory + "/" + filename_without_extension + "_profile.csv";
    if (!std::filesystem::exists(new_directory)) {
        std::filesystem::create_directory(new_directory);
    }
    FileWriter eval_data_writer(eval_filename_out);
    FileWriter profile_writer(profile_filename_out);

    // True nav solution
    NavSolutionNed true_nav_ned;
    if(closed_loop) {
        fdm.Run();
        true_nav_ned = getNavSolFromFdm(fdm); 
    }
    else {
        reader.readNextRow(true_nav_ned);
    }
    NavSolutionEcef true_nav_ecef = nedToEcef(true_nav_ned);   
    NavSolutionEcef true_nav_ecef_old = true_nav_ecef;

    // Old IMU measurements
    ImuMeasurements imu_meas_old;
    imu_meas_old.quant_residuals_f = Vector3::Zero();
    imu_meas_old.quant_residuals_omega = Vector3::Zero();

    // Images
    cv::Mat img_old;
    cv::Mat img;

    // Time of last KF update
    ScalarType time_last_update = -1.0;

    // Init GNSS range biases
    SatPosVel sat_pos_vel_0 = satellitePositionsAndVelocities(true_nav_ned.time,  gnss_config);
    Eigen::Matrix<ScalarType, Eigen::Dynamic, 1, 0, Constants<ScalarType>::kMaxGnssSatellites> gnss_biases = initializeGnssBiases(true_nav_ecef,
                                                                                                        true_nav_ned,
                                                                                                        sat_pos_vel_0,
                                                                                                        gnss_config,
                                                                                                        gen);
    // GNSS measurements at t0
    GnssMeasurements gnss_meas_t0 = generateGnssMeasurements(true_nav_ned.time,
                                                            sat_pos_vel_0,
                                                            true_nav_ned,
                                                            true_nav_ecef,
                                                            gnss_biases, 
                                                            gnss_config,
                                                            gen);

    // Init navigation filter
    StateEstEcef state_est_ecef_init = initStateFromGroundTruth(true_nav_ecef, kf_config, gnss_meas_t0, gen);
    // state_est_ecef_init.nav_sol = true_nav_ecef; // Cheating: ideal init - use to debug INS
    NavKF nav_filter(state_est_ecef_init, kf_config);

    while (true) {

        // ========= Get ground truth ============

        if(closed_loop) {
            if(fdm.GetSimTime() > fdm_end_time) break;
            // Run sim step
            fdm.Run();
            true_nav_ned = getNavSolFromFdm(fdm);
        }
        else {
            if(!reader.readNextRow(true_nav_ned)) break;
        }
        
        true_nav_ecef = nedToEcef(true_nav_ned);

        ScalarType tor_i = true_nav_ned.time - nav_filter.getTime();

        // ========== IMU Simulation ==========

        // Get true specific force and angular rates
        ImuMeasurements ideal_imu_meas = kinematicsEcef(true_nav_ecef, true_nav_ecef_old);

        // Get imu measurements by applying IMU model
        ImuMeasurements imu_meas = imuModel(ideal_imu_meas, imu_meas_old, imu_errors, tor_i, gen);
        // imu_meas = ideal_imu_meas; // Debug: ideal IMU
        
        // ========== Predict ==========

        if (sim_type == SimType::INS_GNSS_TC) {
            nav_filter.tcPredict(imu_meas, tor_i);
        }
        else {
            nav_filter.lcPredict(imu_meas, tor_i);
        }

        // Set camera pose
        viewer_interface.setPose(true_nav_ecef.C_b_e, true_nav_ecef.r_eb_e);

        // ========== Update =========

        ScalarType tor_s = true_nav_ned.time - time_last_update;
        if(tor_s >= epoch_interval && sim_type != SimType::INS) {

            // Simulate GNSS measurements
            SatPosVel sat_pos_vel = satellitePositionsAndVelocities(true_nav_ned.time, gnss_config);
            GnssMeasurements gnss_meas = generateGnssMeasurements(true_nav_ned.time,
                                                                    sat_pos_vel,
                                                                    true_nav_ned,
                                                                    true_nav_ecef,
                                                                    gnss_biases, 
                                                                    gnss_config,
                                                                    gen);
            // Get camera image
            img = viewer_interface.getSnapShot();
            cv::imshow("Synthetic Camera Image", img);
            cv::waitKey(1);
            
            // Loose GNSS Update
            if(sim_type == SimType::INS_GNSS_LC) {
                nav_filter.lcUpdateGnssEcef(gnss_meas, gnss_config);
            }
            // Tight GNSS Update
            else if(sim_type == SimType::INS_GNSS_TC) {
                nav_filter.tcUpdateGnssEcef(gnss_meas, tor_s);
            }
            // Loose position update
            else if(sim_type == SimType::INS_POS) {
                // Simulate position + attitude sensor measurement
                PosMeasEcef pos_meas_ecef = genericPosSensModel(true_nav_ecef, 10.0, gen);
                nav_filter.lcUpdatePosEcef(pos_meas_ecef);
            }
            // Loose position + attitude update
            else if(sim_type == SimType::INS_POS_ROT) {
                // Simulate position + attitude sensor measurement
                PosRotMeasEcef pos_rot_meas_ecef = genericPosRotSensModel(true_nav_ecef, 10.0, 0.01, gen);
                nav_filter.lcUpdatePosRotEcef(pos_rot_meas_ecef);
            }
            else if (sim_type == SimType::INS_OF)
            {   

                if(true_nav_ned.time < 100.0) {
                    // Simulate position + attitude sensor measurement
                    PosRotMeasEcef pos_rot_meas_ecef = genericPosRotSensModel(true_nav_ecef, 10.0, 0.01, gen);
                    nav_filter.lcUpdatePosRotEcef(pos_rot_meas_ecef);
                }
                else{

                    // Fuse Optical Flow
                    if (!img_old.empty()) {
                        nav_filter.lcUpdateOptFlow(img_old, img, K, C_b_c);
                    }
                }
                
            }
            else {
                LOG(ERROR) << "Unknown simulation type";
                return 1;
            }
            
            time_last_update = true_nav_ned.time;
        }

        // ========== Write Results ==========

        StateEstEcef state_est_ecef = nav_filter.getStateEst();
        EvalDataEcef eval_data_ecef = getEvalDataEcef(state_est_ecef, true_nav_ecef);
        eval_data_writer.writeEvalDataRow(eval_data_ecef);
        profile_writer.writeProfileRow(true_nav_ned);
        
        // ============= Update simulation state ==============

        true_nav_ecef_old = true_nav_ecef;
        imu_meas_old = imu_meas;
        img_old = img;

        if(closed_loop) {

            // ============== Close the guidance loop ==============

            NavSolutionNed est_nav_ned = ecefToNed(state_est_ecef.nav_sol);
            Vector3 rpy = dcmToEuler(est_nav_ned.C_b_n);

            // Autopilot properties
            fdm.SetPropertyValue("estimation/attitude/phi-rad", wrapPi(rpy(0))); // Roll
            fdm.SetPropertyValue("estimation/attitude/heading-true-rad", wrapPi(rpy(2))); // Heading
            fdm.SetPropertyValue("estimation/position/h-sl-ft", est_nav_ned.height / kFeetToMeters); // Altitude (feet) 
            fdm.SetPropertyValue("estimation/velocities/h-dot-fps", - est_nav_ned.v_eb_n(2) / kFeetToMeters); // Vertical speed (feet/s)

            // Guidance properties
            fdm.SetPropertyValue("estimation/position/lat-geod-rad", wrapPi(est_nav_ned.latitude));
            fdm.SetPropertyValue("estimation/position/long-geod-rad", wrapPi(est_nav_ned.longitude));
        }

    }

    LOG(INFO) << "Done!";

    return 0;
}
