#ifndef SIMULATION_H
#define SIMULATION_H

#include <eigen3/Eigen/Dense>

#ifdef BUILD_VISION

#include <osgViewer/Viewer>
#include <osgViewer/config/SingleWindow>
#include <osgEarth/ExampleResources>
#include <osgEarth/MapNode>
#include "egl_graphics_window_embedded/egl_graphics_window_embedded.h"

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>

using namespace osgEarth;
using namespace osgEarth::Util;

#endif

#include "constants_types.h"
#include "helpers.h"

namespace intnavlib {

    /// \defgroup simulation Simulation
    /// IMU, GNSS, and generic position / pose sensor behavioral models and utilities.
    /// @{

    template<typename T>
    struct Simulation {

        using Vector3 = Eigen::Matrix<T,3,1>;
        using Vector2 = Eigen::Matrix<T,2,1>;
        using Matrix3 = Eigen::Matrix<T,3,3>;

        using ImuMeasurements = typename Types<T>::ImuMeasurements;
        using PosMeasEcef = typename Types<T>::PosMeasEcef;
        using PosRotMeasEcef = typename Types<T>::PosRotMeasEcef;
        using GnssPosVelMeasEcef = typename Types<T>::GnssPosVelMeasEcef;
        using ImuErrors = typename Types<T>::ImuErrors;
        using NavSolutionNed = typename Types<T>::NavSolutionNed;
        using NavSolutionEcef = typename Types<T>::NavSolutionEcef;
        using StateEstEcef = typename Types<T>::StateEstEcef;
        using ErrorsNed = typename Types<T>::ErrorsNed;
        using EvalDataEcef = typename Types<T>::EvalDataEcef;
        using KfConfig = typename Types<T>::KfConfig;
        using GnssConfig = typename Types<T>::GnssConfig;
        using SatPosVel = typename Types<T>::SatPosVel;
        using GnssMeasurements = typename Types<T>::GnssMeasurements;
        using GnssLsPosVelClock = typename Types<T>::GnssLsPosVelClock;

        static constexpr auto kMaxGnssSatellites = Constants<T>::kMaxGnssSatellites;
        static constexpr auto kEpsilon = Constants<T>::kEpsilon;
        static constexpr auto kR0 = Constants<T>::kR0;
        static constexpr auto kEccentricity = Constants<T>::kEccentricity;
        static constexpr auto kOmega_ie = Constants<T>::kOmega_ie;
        static constexpr auto kGravConst = Constants<T>::kGravConst;
        static constexpr auto kJ2 = Constants<T>::kJ2;
        static constexpr auto kC = Constants<T>::kC;
        static constexpr auto kDegToRad = Constants<T>::kDegToRad;
        static constexpr auto kRadToDeg = Constants<T>::kRadToDeg;
        static constexpr auto kMuGToMetersPerSecondSquared = Constants<T>::kMuGToMetersPerSecondSquared;

        /// \brief Calculates true IMU measurements (specific force and angular velocity) from kinematic states.
        /// This function determines the IMU outputs that would be observed given the change in navigation states.
        /// \param[in] new_nav The navigation solution at the current time step.
        /// \param[in] old_nav The navigation solution at the previous time step.
        /// \return An ImuMeasurements structure containing the true specific force and angular velocity.
        static ImuMeasurements kinematicsEcef(const NavSolutionEcef & new_nav,
                                        const NavSolutionEcef & old_nav);

        /// \brief Models IMU sensor outputs by adding errors (biases, noise, quantization) to true IMU measurements.
        /// \param[in] true_imu_meas The true IMU measurements (specific force and angular velocity).
        /// \param[in] old_imu_meas The IMU measurements from the previous time step, used for quantization residuals.
        /// \param[in] imu_errors Structure containing IMU error parameters (biases, scale factors, noise PSDs, quantization levels).
        /// \param[in] tor_i The integration time interval (delta_t) in seconds.
        /// \param[in] gen A Mersenne Twister random number generator for noise generation.
        /// \return An ImuMeasurements structure containing the modeled IMU outputs with errors.
        static ImuMeasurements imuModel(const ImuMeasurements & true_imu_meas, 
                                const ImuMeasurements & old_imu_meas,
                                const ImuErrors & imu_errors,
                                const T & tor_i,
                                std::mt19937 & gen);

        /// \brief Models a generic ECEF position sensor by adding Gaussian noise to the true position.
        /// \param[in] true_nav The true navigation solution in ECEF frame.
        /// \param[in] pos_sigma The standard deviation of the position measurement noise (in meters).
        /// \param[in] gen A Mersenne Twister random number generator for noise generation.
        /// \return A PosMeasEcef structure containing the noisy position measurement and its covariance.
        static PosMeasEcef genericPosSensModel(const NavSolutionEcef & true_nav,
                                    const T & pos_sigma,
                                    std::mt19937 & gen);

        /// \brief Models a generic ECEF pose (position and rotation) sensor by adding Gaussian noise to the true position and rotation.
        /// \param[in] true_nav The true navigation solution in ECEF frame.
        /// \param[in] pos_sigma The standard deviation of the position measurement noise (in meters).
        /// \param[in] rot_sigma The standard deviation of the rotation measurement noise (in radians).
        /// \param[in] gen A Mersenne Twister random number generator for noise generation.
        /// \return A PosRotMeasEcef structure containing the noisy position and rotation measurement and its covariance.
        static PosRotMeasEcef genericPosRotSensModel(const NavSolutionEcef & true_nav,
                                            const T & pos_sigma,
                                            const T & rot_sigma,
                                            std::mt19937 & gen);

        /// \brief Generates satellite positions and velocities in the ECEF frame based on a simplified circular orbit model.
        /// \param[in] time Current simulation time in seconds.
        /// \param[in] gnss_config Configuration parameters for the GNSS constellation.
        /// \return A SatPosVel structure containing the ECEF Cartesian positions and velocities of all satellites.
        static SatPosVel satellitePositionsAndVelocities(const T & time, const GnssConfig & gnss_config);

        /// \brief Generates GNSS pseudo-range and pseudo-range rate measurements for satellites above the elevation mask angle.
        /// This function also includes satellite positions and velocities in the output.
        /// \param[in] time Current simulation time in seconds.
        /// \param[in] gnss_pos_vel Satellite positions and velocities.
        /// \param[in] true_nav_ned True navigation solution in NED frame.
        /// \param[in] true_nav_ecef True navigation solution in ECEF frame.
        /// \param[in] gnss_biases Pre-calculated GNSS range biases for each satellite.
        /// \param[in] gnss_config Configuration parameters for the GNSS receiver and constellation.
        /// \param[in] gen A Mersenne Twister random number generator for noise generation.
        /// \return A GnssMeasurements structure containing the generated pseudo-range and pseudo-range rate measurements,
        ///         along with corresponding satellite positions and velocities.
        static GnssMeasurements generateGnssMeasurements(const T & time,
                                            const SatPosVel & gnss_pos_vel,
                                            const NavSolutionNed& true_nav_ned,
                                            const NavSolutionEcef& true_nav_ecef,
                                            const Eigen::Matrix<T, Eigen::Dynamic, 1, 0, kMaxGnssSatellites> & gnss_biases, 
                                            const GnssConfig& gnss_config,
                                            std::mt19937 & gen);

        /// \brief Initializes GNSS range biases due to signal-in-space, ionosphere, and troposphere errors.
        /// These biases are calculated based on satellite elevation angles and configured error standard deviations.
        /// \param[in] true_nav_ecef True navigation solution in ECEF frame.
        /// \param[in] true_nav_ned True navigation solution in NED frame.
        /// \param[in] gnss_pos_vel Satellite positions and velocities.
        /// \param[in] gnss_config Configuration parameters for the GNSS constellation and error models.
        /// \param[in] gen A Mersenne Twister random number generator for noise generation.
        /// \return An Eigen::Matrix containing the initialized GNSS range biases for each satellite.
        Eigen::Matrix<T, Eigen::Dynamic, 1, 0, kMaxGnssSatellites>
        static initializeGnssBiases(const NavSolutionEcef & true_nav_ecef,
                                        const NavSolutionNed & true_nav_ned,
                                        const SatPosVel & gnss_pos_vel,
                                        const GnssConfig& gnss_config,
                                        std::mt19937 & gen);

        
        #ifdef BUILD_VISION

        /// \brief Simple OSG manipulator
        class SimpleManipulator : public osgGA::CameraManipulator {
            public:
                inline void setByMatrix(const osg::Matrixd& matrix) override {
                    matrix_ = matrix;
                }

                inline void setByInverseMatrix(const osg::Matrixd& matrix) override {
                    osg::Matrixd inv;
                    inv.invert(matrix);
                    matrix_ = inv;
                }
            
                inline osg::Matrixd getMatrix() const override {
                    return matrix_;
                }
            
                inline osg::Matrixd getInverseMatrix() const override {
                    osg::Matrixd inv;
                    inv.invert(matrix_);
                    return inv;
                }
            
            private:
                osg::Matrixd matrix_;
        };

        /// \brief An interface to the OSG viewer that allows us to set view matrix, proj matrix, get snapshot. It's threadsafe
        class ViewerInterface {

        public:

            ViewerInterface() = delete;

            ViewerInterface(const std::string& earth_file_path,
                            double fx, 
                            double fy, 
                            double cx, 
                            double cy, 
                            unsigned int width, 
                            unsigned int height, 
                            double near, 
                            double far, 
                            unsigned int delay_frames,
                            unsigned int init_delay_frames);
                            
            ~ViewerInterface() {
                viewer->setDone(true);
                if(renderThread.joinable()) {
                    renderThread.join();
                }
            }
            
            // Get a screenshot from the renderer as a cv::Mat
            cv::Mat getSnapShot();
            
            // Set View matrix (camera pose)
            void setPose(const Eigen::Matrix3d &rotation, const Eigen::Vector3d &translation);
            
            // Set projection matrix (camera intrinsics)
            void setProj(const double & fx, 
                        const double & fy, 
                        const double & cx, 
                        const double & cy, 
                        const double & w, 
                        const double & h, 
                        const double & near, 
                        const double & far);

            // Start render thread
            inline void start() {
                renderThread = std::thread(&ViewerInterface::renderLoop, this);
            }

            // starts render thread, but first sets view matrix and waits for scene to load
            void startWithInit(const Eigen::Matrix3d &rotation, 
                                const Eigen::Vector3d &translation);

        private:

            double fx;
            double fy;
            double cx;
            double cy;
            unsigned int width;
            unsigned int height;
            double near;
            double far;
            unsigned int delay_frames;
            unsigned int init_delay_frames;
            
            std::unique_ptr<osg::ArgumentParser> arguments;
            osg::ref_ptr<osg::GraphicsContext> gc;
            osg::ref_ptr<osg::GraphicsContext::Traits> traits;
            osg::ref_ptr<SimpleManipulator> manip;
            osg::ref_ptr<osgViewer::Viewer> viewer;
            osg::ref_ptr<osgViewer::SingleWindow> window;
            osg::ref_ptr<osg::Image> color_image;
            osg::ref_ptr<osg::Image> depth_image;

            std::thread renderThread;

            // mutex & cv to set up producer-consumer model to get renderer snapshot
            std::mutex viewer_mutex;
            std::condition_variable viewer_cv;
            bool viewer_flag;

            void renderLoop ();
        };

        #endif

    };
    
/// @}

// Convenience wrappers

using Simulationd = Simulation<double>;

using Simulationf = Simulation<float>;

// double

inline Typesd::ImuMeasurements kinematicsEcef(const Typesd::NavSolutionEcef & new_nav,
                                const Typesd::NavSolutionEcef & old_nav) {
    return Simulationd::kinematicsEcef(new_nav, old_nav);
}

inline Typesd::ImuMeasurements imuModel(const Typesd::ImuMeasurements & true_imu_meas, 
                        const Typesd::ImuMeasurements & old_imu_meas,
                        const Typesd::ImuErrors & imu_errors,
                        const double & tor_i,
                        std::mt19937 & gen) {
    return Simulationd::imuModel(true_imu_meas, old_imu_meas, imu_errors, tor_i, gen);
}

inline Typesd::PosMeasEcef genericPosSensModel(const Typesd::NavSolutionEcef & true_nav, 
                            const double & pos_sigma,
                            std::mt19937 & gen) {
    return Simulationd::genericPosSensModel(true_nav, pos_sigma, gen);
}

inline Typesd::PosRotMeasEcef genericPosRotSensModel(const Typesd::NavSolutionEcef & true_nav, 
                                    const double & pos_sigma,
                                    const double & rot_sigma,
                                    std::mt19937 & gen) {
    return Simulationd::genericPosRotSensModel(true_nav, pos_sigma, rot_sigma, gen);
}

inline Typesd::SatPosVel satellitePositionsAndVelocities(const double & time, const Typesd::GnssConfig & gnss_config) {
    return Simulationd::satellitePositionsAndVelocities(time, gnss_config);
}

inline Typesd::GnssMeasurements generateGnssMeasurements(const double & time,
                                    const Typesd::SatPosVel & gnss_pos_vel,
                                    const Typesd::NavSolutionNed& true_nav_ned,
                                    const Typesd::NavSolutionEcef& true_nav_ecef,
                                    const Eigen::Matrix<double, Eigen::Dynamic, 1, 0, Constantsd::kMaxGnssSatellites> & gnss_biases, 
                                    const Typesd::GnssConfig& gnss_config,
                                    std::mt19937 & gen) {
    return Simulationd::generateGnssMeasurements(time, gnss_pos_vel, true_nav_ned, true_nav_ecef, gnss_biases, gnss_config, gen);
}

inline Eigen::Matrix<double, Eigen::Dynamic, 1, 0, Constantsd::kMaxGnssSatellites>
initializeGnssBiases(const Typesd::NavSolutionEcef & true_nav_ecef,
                        const Typesd::NavSolutionNed & true_nav_ned,
                        const Typesd::SatPosVel & gnss_pos_vel,
                        const Typesd::GnssConfig& gnss_config,
                        std::mt19937 & gen) {
    return Simulationd::initializeGnssBiases(true_nav_ecef, true_nav_ned, gnss_pos_vel, gnss_config, gen);
}

// float

inline Typesf::ImuMeasurements kinematicsEcef(const Typesf::NavSolutionEcef & new_nav,
                                const Typesf::NavSolutionEcef & old_nav) {
    return Simulationf::kinematicsEcef(new_nav, old_nav);
}

inline Typesf::ImuMeasurements imuModel(const Typesf::ImuMeasurements & true_imu_meas, 
                        const Typesf::ImuMeasurements & old_imu_meas,
                        const Typesf::ImuErrors & imu_errors,
                        const float & tor_i,
                        std::mt19937 & gen) {
    return Simulationf::imuModel(true_imu_meas, old_imu_meas, imu_errors, tor_i, gen);
}

inline Typesf::PosMeasEcef genericPosSensModel(const Typesf::NavSolutionEcef & true_nav, 
                            const float & pos_sigma,
                            std::mt19937 & gen) {
    return Simulationf::genericPosSensModel(true_nav, pos_sigma, gen);
}

inline Typesf::PosRotMeasEcef genericPosRotSensModel(const Typesf::NavSolutionEcef & true_nav, 
                                    const float & pos_sigma,
                                    const float & rot_sigma,
                                    std::mt19937 & gen) {
    return Simulationf::genericPosRotSensModel(true_nav, pos_sigma, rot_sigma, gen);
}

inline Typesf::SatPosVel satellitePositionsAndVelocities(const float & time, const Typesf::GnssConfig & gnss_config) {
    return Simulationf::satellitePositionsAndVelocities(time, gnss_config);
}

inline Typesf::GnssMeasurements generateGnssMeasurements(const float & time,
                                    const Typesf::SatPosVel & gnss_pos_vel,
                                    const Typesf::NavSolutionNed& true_nav_ned,
                                    const Typesf::NavSolutionEcef& true_nav_ecef,
                                    const Eigen::Matrix<float, Eigen::Dynamic, 1, 0, Constantsf::kMaxGnssSatellites> & gnss_biases, 
                                    const Typesf::GnssConfig& gnss_config,
                                    std::mt19937 & gen) {
    return Simulationf::generateGnssMeasurements(time, gnss_pos_vel, true_nav_ned, true_nav_ecef, gnss_biases, gnss_config, gen);
}

inline Eigen::Matrix<float, Eigen::Dynamic, 1, 0, Constantsf::kMaxGnssSatellites>
initializeGnssBiases(const Typesf::NavSolutionEcef & true_nav_ecef,
                        const Typesf::NavSolutionNed & true_nav_ned,
                        const Typesf::SatPosVel & gnss_pos_vel,
                        const Typesf::GnssConfig& gnss_config,
                        std::mt19937 & gen) {
    return Simulationf::initializeGnssBiases(true_nav_ecef, true_nav_ned, gnss_pos_vel, gnss_config, gen);
}

extern template struct Simulation<double>;
extern template struct Simulation<float>;

};

#endif // SIMULATION_H
