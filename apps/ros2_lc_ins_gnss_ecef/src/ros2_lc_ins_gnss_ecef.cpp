#include <chrono>
#include <memory>
#include <queue>
#include <thread>
#include <mutex>

#include "rclcpp/rclcpp.hpp"
#include "sensor_msgs/msg/imu.hpp"
#include "sensor_msgs/msg/nav_sat_fix.hpp"
#include "geometry_msgs/msg/pose_with_covariance_stamped.hpp"
#include <tf2_ros/transform_broadcaster.h>
#include <geometry_msgs/msg/transform_stamped.hpp>
#include <nav_msgs/msg/path.hpp>
#include <geometry_msgs/msg/pose_stamped.hpp>
#include "intnavlib/intnavlib.h"

#define MAX_BUFFER_SIZE_IMU 200
#define MAX_BUFFER_SIZE_GNSS 50
#define TIME_EPSILON_S 0.01

/// @example ros2_lc_ins_gnss_ecef.cpp
/// ROS2 Loosely-coupled INS/GNSS example in ECEF frame

using namespace std::chrono_literals;
using namespace intnavlib;

// Typedefs for convenience

using ScalarType = double;

using Vector3 = Eigen::Matrix<ScalarType,3,1>;
using Vector2 = Eigen::Matrix<ScalarType,2,1>;
using Matrix3 = Eigen::Matrix<ScalarType,3,3>;

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

class InsGnssNode : public rclcpp::Node {

public:

    InsGnssNode() : Node("ins_gnss_node") {

        // Initialize variables from ROS2 parameters
        init();

        done = false;

        // Initialize subscribers and publisher
        imu_sub_ = this->create_subscription<sensor_msgs::msg::Imu>(
            imu_topic_, 100, std::bind(&InsGnssNode::imuCallback, this, std::placeholders::_1));
        gnss_sub_ = this->create_subscription<sensor_msgs::msg::NavSatFix>(
            gnss_topic_, 10, std::bind(&InsGnssNode::gnssCallback, this, std::placeholders::_1));
        pose_pub_ = this->create_publisher<geometry_msgs::msg::PoseWithCovarianceStamped>(out_pose_topic_, 100);
        path_pub_ = this->create_publisher<nav_msgs::msg::Path>(out_pose_topic_ + "_path", 10);
        path_msg_.header.frame_id = "ecef";
        path_counter_ = 0;
        tf_broadcaster_ = std::make_shared<tf2_ros::TransformBroadcaster>(this);

        // Start processing thread
        processing_thread_ = std::thread(&InsGnssNode::processData, this);
    }

    ~InsGnssNode() {
        done = true;
        imu_cv_.notify_all();
        processing_thread_.join();
        if (processing_thread_.joinable()) {
            processing_thread_.join();
        }
    }

private:

    bool done;

    // ROS2 Topics
    std::string imu_topic_;
    std::string gnss_topic_;
    std::string out_pose_topic_;

    // Pubs & Subs
    rclcpp::Subscription<sensor_msgs::msg::Imu>::SharedPtr imu_sub_;
    rclcpp::Subscription<sensor_msgs::msg::NavSatFix>::SharedPtr gnss_sub_;
    rclcpp::Publisher<geometry_msgs::msg::PoseWithCovarianceStamped>::SharedPtr pose_pub_;
    rclcpp::Publisher<nav_msgs::msg::Path>::SharedPtr path_pub_;
    unsigned int path_counter_;
    nav_msgs::msg::Path path_msg_;
    std::shared_ptr<tf2_ros::TransformBroadcaster> tf_broadcaster_;

    // Buffers
    std::queue<sensor_msgs::msg::Imu::SharedPtr> imu_buffer_;
    std::queue<sensor_msgs::msg::NavSatFix::SharedPtr> gnss_buffer_;
    std::mutex imu_buffer_mutex_;
    std::condition_variable imu_cv_;
    std::mutex gnss_buffer_mutex_;

    // Init params
    std::vector<ScalarType> init_lla_;
    std::vector<ScalarType> init_rpy_b_n_;
    std::vector<ScalarType> init_v_eb_n_;

    // Navigation state
    StateEstEcef state_est_ecef_;

    // EKF config
    KfConfig kf_config_;

    std::thread processing_thread_;

    // If log, write results to file
    bool write_log_;
    std::unique_ptr<FileWriter> log_writer_;

    void init() {

        // Declare parameters

        // Topics
        declare_parameter<std::string>("imu_topic");
        declare_parameter<std::string>("gnss_topic");
        declare_parameter<std::string>("out_pose_topic");

        // Logging
        declare_parameter<bool>("write_log");
        declare_parameter<std::string>("log_path");

        // EKF config
        declare_parameter<ScalarType>("kf_config.init_att_unc");
        declare_parameter<ScalarType>("kf_config.init_vel_unc");
        declare_parameter<ScalarType>("kf_config.init_pos_unc");
        declare_parameter<ScalarType>("kf_config.init_b_a_unc");
        declare_parameter<ScalarType>("kf_config.init_b_g_unc");
        declare_parameter<ScalarType>("kf_config.gyro_noise_psd");
        declare_parameter<ScalarType>("kf_config.accel_noise_psd");
        declare_parameter<ScalarType>("kf_config.accel_bias_psd");
        declare_parameter<ScalarType>("kf_config.gyro_bias_psd");

        // Init 
        declare_parameter<std::vector<ScalarType>>("init_lla");
        declare_parameter<std::vector<ScalarType>>("init_v_eb_n");
        declare_parameter<std::vector<ScalarType>>("init_rpy_b_n");

        // Get parameters

        // Topics 
        imu_topic_ = get_parameter("imu_topic").as_string();
        gnss_topic_ = get_parameter("gnss_topic").as_string();
        out_pose_topic_ = get_parameter("out_pose_topic").as_string();

        // Logging
        write_log_ = get_parameter("write_log").as_bool();
        std::string log_path = get_parameter("log_path").as_string();
        if(write_log_){
            log_writer_ = std::make_unique<FileWriter>(log_path);
        }

        // EKF config
        kf_config_.init_att_unc = get_parameter("kf_config.init_att_unc").as_double();
        kf_config_.init_vel_unc = get_parameter("kf_config.init_vel_unc").as_double();
        kf_config_.init_pos_unc = get_parameter("kf_config.init_pos_unc").as_double();
        kf_config_.init_b_a_unc = get_parameter("kf_config.init_b_a_unc").as_double();
        kf_config_.init_b_g_unc = get_parameter("kf_config.init_b_g_unc").as_double();
        kf_config_.gyro_noise_psd = get_parameter("kf_config.gyro_noise_psd").as_double();
        kf_config_.accel_noise_psd = get_parameter("kf_config.accel_noise_psd").as_double();
        kf_config_.accel_bias_psd = get_parameter("kf_config.accel_bias_psd").as_double();
        kf_config_.gyro_bias_psd = get_parameter("kf_config.gyro_bias_psd").as_double();

        // Init nav solution
        init_lla_ = get_parameter("init_lla").as_double_array();
        init_v_eb_n_ = get_parameter("init_v_eb_n").as_double_array();
        init_rpy_b_n_ = get_parameter("init_rpy_b_n").as_double_array();

        NavSolutionNed est_nav_ned = NavSolutionNed{0.0,
                                                    Constants<ScalarType>::kDegToRad * init_lla_[0], 
                                                    Constants<ScalarType>::kDegToRad * init_lla_[1], 
                                                    init_lla_[2], 
                                                    Vector3(init_v_eb_n_[0], init_v_eb_n_[1], init_v_eb_n_[2]), 
                                                    Constants<ScalarType>::kDegToRad * eulerToDcm(Vector3(init_rpy_b_n_[0], init_rpy_b_n_[1], init_rpy_b_n_[2]))};
        
        state_est_ecef_.valid = true;
        state_est_ecef_.nav_sol = nedToEcef(est_nav_ned);
        state_est_ecef_.acc_bias = Vector3::Zero();
        state_est_ecef_.gyro_bias = Vector3::Zero();
        state_est_ecef_.P_matrix = initializePMmatrix(kf_config_);

        // Log loaded parameters for debugging
        
        RCLCPP_INFO(this->get_logger(), "IMU topic: %s", imu_topic_.c_str());
        RCLCPP_INFO(this->get_logger(), "GNSS topic: %s", gnss_topic_.c_str());

        RCLCPP_INFO(this->get_logger(), "Init_lla = %f %f %f", init_lla_[0], init_lla_[1], init_lla_[2]);
        RCLCPP_INFO(this->get_logger(), "Init_v_eb_n = %f %f %f", init_v_eb_n_[0], init_v_eb_n_[1], init_v_eb_n_[2]);
        RCLCPP_INFO(this->get_logger(), "Init_rpy_b_n = %f %f %f", init_rpy_b_n_[0], init_rpy_b_n_[1], init_rpy_b_n_[2]);

        RCLCPP_INFO(this->get_logger(), "Init_att_unc = %f", kf_config_.init_att_unc);
        RCLCPP_INFO(this->get_logger(), "Init_vel_unc = %f", kf_config_.init_vel_unc);
        RCLCPP_INFO(this->get_logger(), "Init_pos_unc = %f", kf_config_.init_pos_unc);

        RCLCPP_INFO(this->get_logger(), "Nav Filter Initialized! Waiting on measurements...");
    }

    void imuCallback(const sensor_msgs::msg::Imu::SharedPtr msg) {
        //  Threadsafe push to buffer, notify condition variable
        std::lock_guard<std::mutex> lock(imu_buffer_mutex_);
        if(imu_buffer_.size() < MAX_BUFFER_SIZE_IMU) {
            imu_buffer_.push(msg);
            imu_cv_.notify_one();
        }
        else RCLCPP_ERROR(this->get_logger(), "IMU buffer full: dropping!");
    }

    void gnssCallback(const sensor_msgs::msg::NavSatFix::SharedPtr msg) {
        //  Threadsafe push to buffer
        std::lock_guard<std::mutex> lock(gnss_buffer_mutex_);
        if(gnss_buffer_.size() < MAX_BUFFER_SIZE_GNSS)
            gnss_buffer_.push(msg);
        else RCLCPP_ERROR(this->get_logger(), "GNSS buffer full: dropping!");
    }

    void processData() {
        while (rclcpp::ok() && !done) {

            sensor_msgs::msg::Imu::SharedPtr imu_msg;
            sensor_msgs::msg::NavSatFix::SharedPtr gnss_msg;

            // Threadsafe read from message buffers
            // Want to get oldest IMU measurement (front)
            ScalarType imu_time = 0;
            {   
                // Wait on IMU buffer
                std::unique_lock<std::mutex> lock(imu_buffer_mutex_);
                // If exceed wait time, retry loop
                if (!imu_cv_.wait_for(lock, 5s, [this] { return !imu_buffer_.empty();})) {
                    continue;  // No data received, retry loop
                }

                // Get measurement timestamp
                rclcpp::Time imu_stamp(imu_buffer_.front()->header.stamp);
                imu_time = imu_stamp.seconds();
                
                // If new, use measurement, else just pop
                if (imu_time > state_est_ecef_.nav_sol.time)
                    imu_msg = imu_buffer_.front();
                imu_buffer_.pop();
                
            }

            // If msg received but not good (old), continue
            if(!imu_msg) continue;

            // We received a new IMU measurement: do propagation

            RCLCPP_INFO(this->get_logger(), "Predict: %f", imu_time);

            ImuMeasurements imu_meas;
            imu_meas.f = Vector3(
                imu_msg->linear_acceleration.x,
                imu_msg->linear_acceleration.y,
                imu_msg->linear_acceleration.z);
            imu_meas.omega = Vector3(
                imu_msg->angular_velocity.x,
                imu_msg->angular_velocity.y,
                imu_msg->angular_velocity.z);
            imu_meas.time = imu_time;

            ScalarType tor_i = imu_meas.time - state_est_ecef_.nav_sol.time;

            // Predict
            state_est_ecef_ = lcPredictKF(state_est_ecef_, imu_meas, kf_config_, tor_i);

            // Now, see if we also have a new GNSS measurement

            {
                std::lock_guard<std::mutex> lock(gnss_buffer_mutex_);

                if (!gnss_buffer_.empty()) {
                    rclcpp::Time gnss_stamp(gnss_buffer_.front()->header.stamp);
                    // RCLCPP_INFO(this->get_logger(), "GNSS time %f", gnss_stamp.seconds());
                    // RCLCPP_INFO(this->get_logger(), "IMU time %f", imu_time);
                    if (abs(gnss_stamp.seconds() - imu_time) < TIME_EPSILON_S)
                        gnss_msg = gnss_buffer_.front();
                    gnss_buffer_.pop();
                    
                }
            }

            // If new GNSS, do update

            if (gnss_msg) {

                rclcpp::Time gnss_stamp(gnss_msg->header.stamp);
                ScalarType gnss_time = gnss_stamp.seconds();
                RCLCPP_INFO(this->get_logger(), "GNSS Update: %f", gnss_time);

                PosMeasEcef gnss_meas;
                gnss_meas.time = gnss_time;
                gnss_meas.r_eb_e = nedToEcef(NavSolutionNed{0,gnss_msg->latitude, gnss_msg->longitude, gnss_msg->altitude, Vector3::Zero(), Matrix3::Identity()}).r_eb_e;
                gnss_meas.cov_mat = Matrix3::Identity() * std::pow(2.5,2); // Covariance matrix

                state_est_ecef_ = lcUpdateKFPosEcef(gnss_meas, state_est_ecef_);
                
            }

            // Publish estimated pose
            publishPose();

            // Log nav sol to file
            if(write_log_){
                NavSolutionNed est_nav_ned = ecefToNed(state_est_ecef_.nav_sol);
                log_writer_->writeProfileRow(est_nav_ned);
            }
        }
    }

    void publishPose() {

        // ===== Pose + Cov stamped =====

        auto pose_msg = geometry_msgs::msg::PoseWithCovarianceStamped();
        pose_msg.header.stamp = rclcpp::Time(static_cast<int64_t>(state_est_ecef_.nav_sol.time * 1e9));
        pose_msg.header.frame_id = "ecef";

        // Fill in position
        pose_msg.pose.pose.position.x = state_est_ecef_.nav_sol.r_eb_e[0];
        pose_msg.pose.pose.position.y = state_est_ecef_.nav_sol.r_eb_e[1];
        pose_msg.pose.pose.position.z = state_est_ecef_.nav_sol.r_eb_e[2];

        // Fill in orientation (convert rotation matrix to quaternion)
        Eigen::Quaterniond q(state_est_ecef_.nav_sol.C_b_e);
        pose_msg.pose.pose.orientation.x = q.x();
        pose_msg.pose.pose.orientation.y = q.y();
        pose_msg.pose.pose.orientation.z = q.z();
        pose_msg.pose.pose.orientation.w = q.w();

        Eigen::Matrix<ScalarType, 6, 6> pose_cov;
        pose_cov.block<3,3>(0,0) = state_est_ecef_.P_matrix.block<3,3>(6,6);  // position covariance
        pose_cov.block<3,3>(0,3) = state_est_ecef_.P_matrix.block<3,3>(6,0);  // cross-covariance
        pose_cov.block<3,3>(3,0) = state_est_ecef_.P_matrix.block<3,3>(0,6);  // cross-covariance
        pose_cov.block<3,3>(3,3) = state_est_ecef_.P_matrix.block<3,3>(0,0);  // attitude covariance

        // Copy the 6x6 matrix into the covariance array (which is in row-major order)
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                pose_msg.pose.covariance[i * 6 + j] = pose_cov(i, j);
            }
        }

        pose_pub_->publish(pose_msg);

        // ====== TF =======

        geometry_msgs::msg::TransformStamped tf_msg;
        tf_msg.header.stamp = pose_msg.header.stamp;
        tf_msg.header.frame_id = "ecef";  // Parent frame
        tf_msg.child_frame_id = "base";  // or whatever child frame makes sense

        tf_msg.transform.translation.x = state_est_ecef_.nav_sol.r_eb_e[0];
        tf_msg.transform.translation.y = state_est_ecef_.nav_sol.r_eb_e[1];
        tf_msg.transform.translation.z = state_est_ecef_.nav_sol.r_eb_e[2];
        tf_msg.transform.rotation.x = q.x();
        tf_msg.transform.rotation.y = q.y();
        tf_msg.transform.rotation.z = q.z();
        tf_msg.transform.rotation.w = q.w();

        tf_broadcaster_->sendTransform(tf_msg);

        // ====== Path ======

        // Erase oldest pose if array too big
        if (path_msg_.poses.size() > 1000) {
            path_msg_.poses.erase(path_msg_.poses.begin());
        }

        // Publish only every nth pose
        path_counter_ ++;
        if (path_counter_ >= 100) {
            path_counter_ = 0;
            geometry_msgs::msg::PoseStamped pose_stamped;
            pose_stamped.header = pose_msg.header;
            pose_stamped.pose = pose_msg.pose.pose;
            path_msg_.header = pose_msg.header;
            path_msg_.poses.push_back(pose_stamped);
            path_pub_->publish(path_msg_);
        }
    }
};

int main(int argc, char *argv[]) {
    rclcpp::init(argc, argv);
    rclcpp::spin(std::make_shared<InsGnssNode>());
    rclcpp::shutdown();
    return 0;
}
