#ifndef NAVIGATION_H
#define NAVIGATION_H

#include <boost/math/distributions/chi_squared.hpp>

#include "constants_types.h"
#include "helpers.h"

#ifdef BUILD_VISION

#include <glog/logging.h>
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/video/tracking.hpp>
#include <opencv2/features2d.hpp>
#include <opencv2/calib3d.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/core/eigen.hpp>

#endif

namespace intnavlib {

/// \defgroup navigation Navigation
/// Navigation utilities.
/// @{

template<typename T>
struct Navigation {

    using Vector3 = Eigen::Matrix<T,3,1>;
    using Vector4 = Eigen::Matrix<T,4,1>;
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

    /// \brief Predicts the state estimation error covariance matrix.
    /// This function propagates the error covariance matrix using the discrete-time Kalman filter prediction equation.
    /// \tparam n_x The dimension of the state vector.
    /// \param[in] Phi_matrix The state transition matrix.
    /// \param[in] Q_matrix The system noise covariance matrix.
    /// \param[in] P_matrix_old The error covariance matrix from the previous time step.
    /// \param[out] P_matrix The predicted error covariance matrix for the current time step.
    template<int n_x>
    static void predictKF(const Eigen::Matrix<T, n_x, n_x> & Phi_matrix, 
                    const Eigen::Matrix<T, n_x, n_x> & Q_matrix,
                    const Eigen::Matrix<T, n_x, n_x> & P_matrix_old, 
                    Eigen::Matrix<T, n_x, n_x> & P_matrix);

    /// \brief Generic Kalman filter update
    /// This function performs a generic Kalman filter update step. It calculates the Kalman gain,
    /// updates the error state estimates, and propagates the error covariance matrix.
    /// It also includes a real-time consistency check using the Chi-squared test.
    /// \tparam n_x The dimension of the state vector.
    /// \tparam n_z The dimension of the measurement vector.
    /// \param[in] delta_z The measurement innovation vector.
    /// \param[in] P_matrix The prior error covariance matrix.
    /// \param[in] H_matrix The measurement matrix.
    /// \param[in] R_matrix The measurement noise covariance matrix.
    /// \param[out] S_matrix The innovation covariance matrix.
    /// \param[out] x_est_new The updated error state estimates.
    /// \param[out] P_matrix_post The updated (posterior) error covariance matrix.
    template<int n_x, int n_z, int max_n_z = 2* kMaxGnssSatellites>
    static bool updateKF(const Eigen::Matrix<T, n_z, 1, 0, max_n_z, 1> & delta_z,
                const Eigen::Matrix<T, n_x, n_x> & P_matrix, 
                const Eigen::Matrix<T, n_z, n_x, 0, max_n_z, n_x> & H_matrix,
                const Eigen::Matrix<T, n_z, n_z> & R_matrix,
                const T & p_value,
                Eigen::Matrix<T, n_z, n_z, 0, max_n_z, max_n_z> & S_matrix,
                Eigen::Matrix<T, n_x, 1> & x_est_new,
                Eigen::Matrix<T, n_x, n_x> & P_matrix_post);

    /// \brief Solves the navigation equations in the ECEF frame.
    /// This function propagates the navigation solution (position, velocity, and attitude)
    /// from a previous time step to the current time step using IMU measurements.
    /// The navigation solution update is performed as follows:
    /*!

    \f[
    \mathbf{C}_b^e(+)=\left(\begin{array}{ccc}
    \cos \omega_{i e} \tau_i & \sin \omega_{i e} \tau_i & 0 \\
    -\sin \omega_{i e} \tau_i & \cos \omega_{i e} \tau_i & 0 \\
    0 & 0 & 1
    \end{array}\right) \mathbf{C}_b^e(-) \mathbf{C}_{b+}^{b-} ,
    \f]

    \f[
    \mathbf{f}_{i b}^e=\overline{\mathbf{C}}_b^e \mathbf{f}_{i b}^b, \quad \overline{\mathbf{C}}_b^e=\mathbf{C}_b^e(-) \mathbf{C}_b^{b-}-\frac{1}{2} \boldsymbol{\Omega}_e^e \mathbf{C}_b^e(-) \tau_i ,
    \f]

    \f[
    \mathbf{v}_{e b}^e(+) \approx \mathbf{v}_{e b}^e(-)+\left(\mathbf{f}_{i b}^e+\mathbf{g}_b^e\left(\mathbf{r}_{e b}^e(-)\right)-2 \boldsymbol{\Omega}_{i e}^e \mathbf{v}_{e b}^e(-)\right) \tau_i ,
    \f]

    \f[
    \mathbf{r}_{e b}^e(+)=\mathbf{r}_{e b}^e(-)+\left(\mathbf{v}_{e b}^e(-)+\mathbf{v}_{e b}^e(+)\right) \frac{\tau_i}{2}.
    \f]

    */
    /// \param[in] old_nav The navigation solution at the previous time step.
    /// \param[in] imu_meas IMU measurements (specific force and angular velocity) for the current interval.
    /// \param[in] tor_i The integration time interval (delta_t) in seconds.
    /// \return The updated navigation solution in the ECEF frame.
    static NavSolutionEcef navEquationsEcef(const NavSolutionEcef & old_nav, 
                                const ImuMeasurements & imu_meas,
                                const T & tor_i);

    /// \brief Utility to perform both state prediction and uncertainty propagation (Loose)
    static StateEstEcef lcPredictKF(const StateEstEcef & state_est_old, 
                                const ImuMeasurements & imu_meas,
                                const KfConfig & kf_config,
                                const T & tor_i);

    /// \brief Utility to perform both state prediction and uncertainty propagation (Tight)
    static StateEstEcef tcPredictKF(const StateEstEcef & state_est_old, 
                                const ImuMeasurements & imu_meas,
                                const KfConfig & kf_config,
                                const T & tor_i);

    /// \brief Propagates the error-state uncertainty for a Loosely Coupled (LC) Kalman Filter.
    /// This function updates the error covariance matrix based on the system dynamics
    /// and noise characteristics over a given integration time.
    /*!

    Loosely-coupled ECEF Error state:

    \f[
    \mathbf{x}_{I N S}^e=\left(\begin{array}{c}
    \delta \boldsymbol{\Psi}_{e b}^e \\
    \delta \mathbf{v}_{e b}^e \\
    \delta \mathbf{r}_{e b}^e \\
    \mathbf{b}_a \\
    \mathbf{b}_g
    \end{array}\right)
    \f]

    Loosely-coupled ECEF Error state first-order linear system model:

    \f[
    \boldsymbol{\Phi}_{\text {INS }}^e \approx\left[\begin{array}{ccccc}
    \mathbf{I}_3-\boldsymbol{\Omega}_e^e \tau_s & 0_3 & 0_3 & 0_3 & \hat{\mathbf{C}}_b^e \tau_s \\
    \mathbf{F}_{12}^e \tau_s & \mathbf{I}_3-2 \boldsymbol{\Omega}_e^e \tau_s & \mathbf{F}_{23}^e \tau_s & \hat{\mathbf{C}}_b^e \tau_s & 0_3 \\
    0_3 & \mathbf{I}_3 \tau_s & \mathbf{I}_3 & 0_3 & 0_3 \\
    0_3 & 0_3 & 0_3 & \mathbf{I}_3 & 0_3 \\
    0_3 & 0_3 & 0_3 & 0_3 & \mathbf{I}_3
    \end{array}\right]
    \f]

    \f[
    \mathbf{Q}_{I N S}^e \approx  \mathbf{Q'}_{I N S}^e =\left(\begin{array}{ccccc}
    S_{r g} \mathbf{I}_3 & 0_3 & 0_3 & 0_3 & 0_3 \\
    0_3 & S_{r a} \mathbf{I}_3 & 0_3 & 0_3 & 0_3 \\
    0_3 & 0_3 & 0_3 & 0_3 & 0_3 \\
    0_3 & 0_3 & 0_3 & S_{b a d} \mathbf{I}_3 & 0_3 \\
    0_3 & 0_3 & 0_3 & 0_3 & S_{b g d} \mathbf{I}_3
    \end{array}\right) \tau_s
    \f]

    Loosely-coupled ECEF Error state covariance propagation:

    \f[
    \mathbf{P}_k^{-} \approx \boldsymbol{\Phi}_{k-1}\left(\mathbf{P}_{k-1}^{+}+\frac{1}{2} \mathbf{Q}_{k-1}^{\prime}\right) \boldsymbol{\Phi}_{k-1}^{\mathrm{T}}+\frac{1}{2} \mathbf{Q}_{k-1}^{\prime} .
    \f]

    */
    /// \param[in] P_matrix_old The error covariance matrix from the previous time step.
    /// \param[in] old_nav_est_ecef The estimated navigation solution in ECEF frame at the previous time step.
    /// \param[in] old_nav_est_ned The estimated navigation solution in NED frame at the previous time step.
    /// \param[in] imu_meas IMU measurements used for propagation.
    /// \param[in] lc_kf_config Configuration parameters for the LC Kalman Filter,
    /// including noise PSDs for gyro, accelerometer, and biases.
    /// \param[in] tor_s Integration time (delta_t) in seconds.
    /// \return An Eigen::Matrix<T,15,15> representing the propagated error covariance matrix.
    static Eigen::Matrix<T,15,15> lcPropUnc(const Eigen::Matrix<T,15,15> & P_matrix_old, 
                                            const NavSolutionEcef & old_nav_est_ecef,
                                            const NavSolutionNed & old_nav_est_ned,
                                            const ImuMeasurements & imu_meas,
                                            const KfConfig & lc_kf_config,
                                            const T & tor_s);

    /// \brief Propagates the error-state uncertainty for a Tightly Coupled (TC) Kalman Filter.
    /// This function updates the error covariance matrix based on the system dynamics
    /// and noise characteristics over a given integration time, including clock states.
    /*!

    Tightly-coupled ECEF Error state:

    \f[
    \mathbf{x}_{I N S}^e=\left(\begin{array}{c}
    \delta \boldsymbol{\Psi}_{e b}^e \\
    \delta \mathbf{v}_{e b}^e \\
    \delta \mathbf{r}_{e b}^e \\
    \mathbf{b}_a \\
    \mathbf{b}_g \\
    \delta \rho^a_c \\
    \delta \dot{\rho}^a_c
    \end{array}\right)
    \f]

    Tightly-coupled ECEF Error state first-order linear system model:

    \f[
    \boldsymbol{\Phi}_{\text{INS}}^e \approx \left[\begin{array}{ccccccc}
    \mathbf{I}_3 - \boldsymbol{\Omega}_e^e \tau_s & 0_3 & 0_3 & 0_3 & \hat{\mathbf{C}}_b^e \tau_s & 0_{3 \times 1} & 0_{3 \times 1} \\
    \mathbf{F}_{12}^e \tau_s & \mathbf{I}_3 - 2\boldsymbol{\Omega}_e^e \tau_s & \mathbf{F}_{23}^e \tau_s & \hat{\mathbf{C}}_b^e \tau_s & 0_3 & 0_{3 \times 1} & 0_{3 \times 1} \\
    0_3 & \mathbf{I}_3 \tau_s & \mathbf{I}_3 & 0_3 & 0_3 & 0_{3 \times 1} & 0_{3 \times 1} \\
    0_3 & 0_3 & 0_3 & \mathbf{I}_3 & 0_3 & 0_{3 \times 1} & 0_{3 \times 1} \\
    0_3 & 0_3 & 0_3 & 0_3 & \mathbf{I}_3 & 0_{3 \times 1} & 0_{3 \times 1} \\
    0_{1 \times 3} & 0_{1 \times 3} & 0_{1 \times 3} & 0_{1 \times 3} & 0_{1 \times 3} & 1 & \tau_s \\
    0_{1 \times 3} & 0_{1 \times 3} & 0_{1 \times 3} & 0_{1 \times 3} & 0_{1 \times 3} & 0 & 1 \\
    \end{array}
    \right]
    \f]

    \f[
    \mathbf{Q}_{INS}^e \approx \mathbf{Q'}_{I N S}^e =
    \left(\begin{array}{ccccccc}
    S_{rg} \mathbf{I}_3 & 0_3 & 0_3 & 0_3 & 0_3 & 0_{3 \times 1} & 0_{3 \times 1} \\
    0_3 & S_{ra} \mathbf{I}_3 & 0_3 & 0_3 & 0_3 & 0_{3 \times 1} & 0_{3 \times 1} \\
    0_3 & 0_3 & 0_3 & 0_3 & 0_3 & 0_{3 \times 1} & 0_{3 \times 1} \\
    0_3 & 0_3 & 0_3 & S_{bad} \mathbf{I}_3 & 0_3 & 0_{3 \times 1} & 0_{3 \times 1} \\
    0_3 & 0_3 & 0_3 & 0_3 & S_{bgd} \mathbf{I}_3 & 0_{3 \times 1} & 0_{3 \times 1} \\
    0_{1 \times 3} & 0_{1 \times 3} & 0_{1 \times 3} & 0_{1 \times 3} & 0_{1 \times 3} & S_{cb} & 0 \\
    0_{1 \times 3} & 0_{1 \times 3} & 0_{1 \times 3} & 0_{1 \times 3} & 0_{1 \times 3} & 0 & S_{cd} \\
    \end{array}\right) \tau_s
    \f]

    Tightly-coupled ECEF Error state covariance propagation:

    \f[
    \mathbf{P}_k^{-} \approx \boldsymbol{\Phi}_{k-1}\left(\mathbf{P}_{k-1}^{+}+\frac{1}{2} \mathbf{Q}_{k-1}^{\prime}\right) \boldsymbol{\Phi}_{k-1}^{\mathrm{T}}+\frac{1}{2} \mathbf{Q}_{k-1}^{\prime} .
    \f]
    */
    /// \param[in] P_matrix_old The error covariance matrix from the previous time step.
    /// \param[in] old_nav_est_ecef The estimated navigation solution in ECEF frame at the previous time step.
    /// \param[in] old_nav_est_ned The estimated navigation solution in NED frame at the previous time step.
    /// \param[in] imu_meas IMU measurements used for propagation.
    /// \param[in] tc_kf_config Configuration parameters for the TC Kalman Filter,
    /// including noise PSDs for gyro, accelerometer, biases, and clock states.
    /// \param[in] tor_s Integration time (delta_t) in seconds.
    /// \return An Eigen::Matrix<T,17,17> representing the propagated error covariance matrix.
    static Eigen::Matrix<T,17,17> tcPropUnc(const Eigen::Matrix<T,17,17> & P_matrix_old, 
                                            const NavSolutionEcef & old_nav_est_ecef,
                                            const NavSolutionNed & old_nav_est_ned,
                                            const ImuMeasurements & imu_meas,
                                            const KfConfig & tc_kf_config,
                                            const T & tor_s);

    /// \brief Performs a Loosely Coupled (LC) Kalman Filter update using an ECEF position measurement.
    /// This is a closed-loop error-state KF update that corrects the navigation solution
    /// and IMU bias estimates by integrating a position measurement.
    /// \param[in] pos_meas The position measurement in ECEF frame.
    /// \param[in] state_est_prior The prior estimated state (navigation solution and biases)
    /// and its covariance matrix.
    /// \param[in] p_value P-value threshold for Chi squared consistency check: accept update only if residual in the confidence interval.
    /// \return A StateEstEcefLc structure containing the updated state estimation.
    static StateEstEcef lcUpdateKFPosEcef (const PosMeasEcef & pos_meas, 
                                        const StateEstEcef & state_est_old,
                                        const T & p_value = 0.99);

    /// \brief Performs a Loosely Coupled (LC) Kalman Filter update using a GNSS position and velocity measurement.
    /// This is a closed-loop error-state KF update that corrects the navigation solution
    /// and IMU bias estimates by integrating a GNSS position and velocity measurement.
    /*!

    Error state:

    \f[
    \mathbf{x}_{I N S}^e=\left(\begin{array}{c}
    \delta \boldsymbol{\Psi}_{e b}^e \\
    \delta \mathbf{v}_{e b}^e \\
    \delta \mathbf{r}_{e b}^e \\
    \mathbf{b}_a \\
    \mathbf{b}_g
    \end{array}\right)
    \f]

    Observation matrix:

    \f[
    \mathbf{H}_{G, k}^e=\left(\begin{array}{ccccc}
    \mathbf{H}_{r 1}^e & 0_3 & -\mathbf{I}_3 & 0_3 & 0_3 \\
    \mathbf{H}_{v 1}^e & -\mathbf{I}_3 & 0_3 & 0_3 & \mathbf{H}_{v 5}^e
    \end{array}\right)_k
    \f]

    Innovation:

    \f[
    \delta \mathbf{z}_{G, k}^{e-} =
    \begin{bmatrix}
    \hat{\mathbf{r}}_{e a G}^e - \hat{\mathbf{r}}_{e b}^e - \hat{\mathbf{C}}_b^e \mathbf{l}_{b a}^b \\
    \hat{\mathbf{v}}_{e a G}^e - \hat{\mathbf{v}}_{e b}^e - \hat{\mathbf{C}}_b^e \left( \hat{\boldsymbol{\omega}}_{i b}^b \wedge \mathbf{l}_{b a}^b \right) + \boldsymbol{\Omega}_e^e \hat{\mathbf{C}}_b^e \mathbf{l}_{b a}^b
    \end{bmatrix}_k^{-}
    \f]

    Error state update:

    \f[
    \mathbf{K}_k=\mathbf{P}_k^{-} \mathbf{H}_k^{\mathrm{T}}\left(\mathbf{H}_k \mathbf{P}_k^{-} \mathbf{H}_k^{\mathrm{T}}+\mathbf{R}_k\right)^{-1}
    \f]

    \f[
    \begin{aligned}
    \hat{\mathbf{x}}_k^{+} & =\hat{\mathbf{x}}_k^{-}+\mathbf{K}_k\left(\mathbf{z}_k-\mathbf{H}_k \hat{\mathbf{x}}_k^{-}\right) \\
    & =\hat{\mathbf{x}}_k^{-}+\mathbf{K}_k \delta \mathbf{z}_k^{-}
    \end{aligned}
    \f]

    \f[
    \mathbf{P}_k^{+}=\left(\mathbf{I}-\mathbf{K}_k \mathbf{H}_k\right) \mathbf{P}_k^{-}
    \f]

    Navigation solution update

    \f[
    \begin{gathered}
    \hat{\mathbf{C}}_b^\gamma(+)=\delta \hat{\mathbf{C}}_b^{\gamma^{\mathrm{T}}} \hat{\mathbf{C}}_b^\gamma(-) \approx\left(\mathbf{I}_3-\left[\delta \hat{\boldsymbol{\psi}}_{\gamma b}^\gamma \wedge\right]\right) \hat{\mathbf{C}}_b^\gamma(-) \\
    \hat{\mathbf{v}}_{\beta b}^\gamma(+)=\hat{\mathbf{v}}_{\beta b}^\gamma(-)-\delta \hat{\mathbf{v}}_{\beta b}^\gamma \\
    \hat{\mathbf{r}}_{\beta b}^\gamma(+)=\hat{\mathbf{r}}_{\beta b}^\gamma(-)-\delta \hat{\mathbf{r}}_{\beta b}^\gamma
    \end{gathered}
    \f]

    */
    /// \param[in] gnss_meas The GNSS measurements, including pseudo-ranges, pseudo-range rates,
    /// satellite positions, and satellite velocities.
    /// \param[in] state_est_prior The prior estimated state (navigation solution and biases) and its covariance matrix.
    /// \param[in] gnss_config GNSS configuration.
    /// \param[in] p_value P-value threshold for Chi squared consistency check: accept update only if residual in the confidence interval.
    /// \return A StateEstEcefLc structure containing the updated state estimation.
    static StateEstEcef lcUpdateKFGnssEcef (const GnssMeasurements & gnss_meas, 
                                        const StateEstEcef & state_est_old,
                                        const GnssConfig & gnss_config,
                                        const T & p_value = 0.99);

    /// \brief Performs a Tightly Coupled (TC) Kalman Filter update using GNSS pseudo-range and pseudo-range rate measurements.
    /// This is a closed-loop error-state KF update that corrects the navigation solution,
    /// IMU bias estimates, and receiver clock states by integrating raw GNSS measurements.
    /*!

    Tightly-coupled ECEF Error state:

    \f[
    \mathbf{x}_{I N S}^e=\left(\begin{array}{c}
    \delta \boldsymbol{\Psi}_{e b}^e \\
    \delta \mathbf{v}_{e b}^e \\
    \delta \mathbf{r}_{e b}^e \\
    \mathbf{b}_a \\
    \mathbf{b}_g \\
    \delta \rho^a_c \\
    \delta \dot{\rho}^a_c
    \end{array}\right)
    \f]

    Observation matrix:

    \f[
    \mathbf{H}_{G, k}^\gamma=\left(\begin{array}{ccccccc}
    \frac{\partial \mathbf{z}_\rho}{\partial \delta \boldsymbol{\psi}_{\gamma b}^\gamma} & 0_{m, 3} & \frac{\partial \mathbf{z}_\rho}{\partial \delta \mathbf{r}_{\gamma b}^\gamma} & 0_{m, 3} & 0_{m, 3} & \frac{\partial \mathbf{z}_\rho}{\partial \delta \rho_c^a} & 0_{m, 1} \\
    \frac{\partial \mathbf{z}_r}{\partial \delta \boldsymbol{\psi}_{\gamma b}^\gamma} & \frac{\partial \mathbf{z}_r}{\partial \delta \mathbf{v}_{\gamma b}^\gamma} & \frac{\partial \mathbf{z}_r}{\partial \delta \mathbf{r}_{\gamma b}^\gamma} & 0_{m, 3} & \frac{\partial \mathbf{z}_\rho}{\partial \mathbf{b}_g} & 0_{m, 1} & \frac{\partial \mathbf{z}_r}{\partial \delta \dot{\rho}_c^a}
    \end{array}\right)_{\mathbf{x}=\hat{\mathbf{x}}_{\bar{k}}}
    \f]

    Innovation:

    \f[
    \delta \mathbf{z}_{\bar{G}, k}=\binom{\delta \mathbf{z}_{\rho, k}^{-}}{\delta \mathbf{z}_{r, k}^{-}}, where \quad \delta \mathbf{z}_{\rho, k}^{-}=\left(\tilde{\rho}_{a, C}^1-\hat{\rho}_{a, C}^{1-}, \tilde{\rho}_{a, C}^2-\hat{\rho}_{a, C}^{2-}, \cdots \tilde{\rho}_{a, C}^m-\hat{\rho}_{a, C}^{m-}\right)_k , \quad \delta \mathbf{z}_{r, k}^{-}=\left(\tilde{\dot{\rho}}_{a, C}^1-\hat{\hat{\rho}}_{a, C}^{1-}, \tilde{\hat{\rho}}_{a, C}^2-\hat{\hat{\rho}}_{a, C}^{2-}, \cdots \tilde{\hat{\rho}}_{a, C}^m-\hat{\hat{\rho}}_{a, C}^{m-}\right)_k .
    \f]

    Error state update:

    \f[
    \mathbf{K}_k=\mathbf{P}_k^{-} \mathbf{H}_k^{\mathrm{T}}\left(\mathbf{H}_k \mathbf{P}_k^{-} \mathbf{H}_k^{\mathrm{T}}+\mathbf{R}_k\right)^{-1}
    \f]

    \f[
    \begin{aligned}
    \hat{\mathbf{x}}_k^{+} & =\hat{\mathbf{x}}_k^{-}+\mathbf{K}_k\left(\mathbf{z}_k-\mathbf{H}_k \hat{\mathbf{x}}_k^{-}\right) \\
    & =\hat{\mathbf{x}}_k^{-}+\mathbf{K}_k \delta \mathbf{z}_k^{-}
    \end{aligned}
    \f]

    \f[
    \mathbf{P}_k^{+}=\left(\mathbf{I}-\mathbf{K}_k \mathbf{H}_k\right) \mathbf{P}_k^{-}
    \f]

    Navigation solution update

    \f[
    \begin{gathered}
    \hat{\mathbf{C}}_b^\gamma(+)=\delta \hat{\mathbf{C}}_b^{\gamma^{\mathrm{T}}} \hat{\mathbf{C}}_b^\gamma(-) \approx\left(\mathbf{I}_3-\left[\delta \hat{\boldsymbol{\psi}}_{\gamma b}^\gamma \wedge\right]\right) \hat{\mathbf{C}}_b^\gamma(-) \\
    \hat{\mathbf{v}}_{\beta b}^\gamma(+)=\hat{\mathbf{v}}_{\beta b}^\gamma(-)-\delta \hat{\mathbf{v}}_{\beta b}^\gamma \\
    \hat{\mathbf{r}}_{\beta b}^\gamma(+)=\hat{\mathbf{r}}_{\beta b}^\gamma(-)-\delta \hat{\mathbf{r}}_{\beta b}^\gamma
    \end{gathered}
    \f]

    */
    /// \param[in] gnss_meas The GNSS measurements, including pseudo-ranges, pseudo-range rates,
    /// satellite positions, and satellite velocities.
    /// \param[in] state_est_prior The prior estimated state (navigation solution, biases, and clock states)
    /// and its covariance matrix.
    /// \param[in] tor_s Elapsed time since last update in seconds.
    /// \param[in] p_value P-value threshold for Chi squared consistency check: accept update only if residual in the confidence interval.
    /// \return A StateEstEcefTc structure containing the updated state estimation.
    static StateEstEcef tcUpdateKFGnssEcef(const GnssMeasurements & gnss_meas, 
                                        const StateEstEcef & state_est_prior,
                                        const T & tor_s,
                                        const T & p_value = 0.99);

    /// \brief Performs a Loosely Coupled (LC) Kalman Filter update using an ECEF position and rotation measurement.
    /// This is a closed-loop error-state KF update that corrects the navigation solution
    /// and IMU bias estimates by integrating a position and rotation measurement.
    /// \param[in] pos_rot_meas The position and rotation measurement in ECEF frame.
    /// \param[in] state_est_prior The prior estimated state (navigation solution and biases)
    /// and its covariance matrix.
    /// \param[in] p_value P-value threshold for Chi squared consistency check: accept update only if residual in the confidence interval.
    /// \return A StateEstEcefLc structure containing the updated state estimation.
    static StateEstEcef lcUpdateKFPosRotEcef(const PosRotMeasEcef & pos_rot_meas, 
                                        const StateEstEcef & state_est_old,
                                        const T & p_value = 0.99);

    /// \brief Performs a Loosely Coupled (LC) Kalman filter optical flow update.
    /// Inspired by https://hilandtom.com/tombotterill/Hide-Botterill-Andreotti-ION10.pdf
    /// \param[in] img_0 The first image (previous frame).
    /// \param[in] img_1 The second image (current frame).
    /// \param[in] K The camera intrinsic matrix.
    /// \param[in] C_c_b The rotation matrix from body frame to camera frame.
    /// \param[in] state_est_prior The prior estimated state (navigation solution and biases)
    /// and its covariance matrix.
    /// \param[in] p_value P-value threshold for Chi squared consistency check: accept update only if residual in the confidence interval.
    /// \return A StateEstEcefLc structure containing the updated state estimation.
    static StateEstEcef lcUpdateKFOptFlow(const cv::Mat & img_0,
                                    const cv::Mat & img_1,
                                    const Matrix3 & K,
                                    const Matrix3 & C_c_b,
                                    const StateEstEcef & state_est_prior,
                                    const T & p_value = 0.99);

    /// \brief Calculates position, velocity, clock offset, and clock drift rate using unweighted iterated least squares.
    /// \param[in] gnss_meas Structure containing GNSS pseudo-range and pseudo-range rate measurements,
    ///                              along with satellite positions and velocities.
    /// \param[in] prior_r_ea_e Prior estimated antenna position in ECEF frame.
    /// \param[in] prior_v_ea_e Prior estimated antenna velocity in ECEF frame.
    /// \return Structure containing the estimated antenna position, velocity, clock offset and drift rate.
    static GnssLsPosVelClock gnssLsPositionVelocityClock(const GnssMeasurements & gnss_meas,
                                                const Vector3 & prior_r_ea_e,
                                                const Vector3 & prior_v_ea_e);

    /// \brief Calculates position, velocity using unweighted iterated least squares.
    /// \param[in] gnss_meas Structure containing GNSS pseudo-range and pseudo-range rate measurements,
    ///                              along with satellite positions and velocities.
    /// \param[in] prior_r_ea_e Prior estimated antenna position in ECEF frame.
    /// \param[in] prior_v_ea_e Prior estimated antenna velocity in ECEF frame.
    /// \param[in] gnss_config GNSS config parameters
    /// \return Structure containing the estimated antenna position, velocity.
    static GnssPosVelMeasEcef gnssLsPositionVelocity(const GnssMeasurements & gnss_meas,
                                    const Vector3 & prior_r_ea_e,
                                    const Vector3 & prior_v_ea_e,
                                    const GnssConfig & gnss_config);

    /// \brief Initializes the error covariance matrix for a navigation Kalman Filter.
    /// \param[in] tc_kf_config Configuration parameters for the Kalman Filter,
    /// including initial uncertainties for attitude, velocity, position, biases,
    /// and clock offset/drift.
    /// \return An Eigen::Matrix<T,17,17> representing the initialized error covariance matrix.
    static Eigen::Matrix<T,17,17> initializePMmatrix(const KfConfig & tc_kf_config);

    /// \brief Init navigation filter state from ground truth.
    /// \param[in] true_nav_ecef Ground truth navigation solution.
    /// \param[in] kf_config Navigation filter configuration.
    /// \param[in] gnss_meas GNSS measurements.
    /// \param[in] gen Random source.
    /// \return Navigation filter state.
    static StateEstEcef initStateFromGroundTruth(const NavSolutionEcef & true_nav_ecef, const KfConfig & kf_config, const GnssMeasurements & gnss_meas, std::mt19937 & gen);

    /// \brief Nav Kalman filter class.
    class NavKF {
        private:
            /// \brief State estimate
            StateEstEcef state_est_;
            /// \brief KF config.
            KfConfig kf_config_;
        public:
            /// \brief Deleted default constructor.
            NavKF() = delete;
            NavKF(const NavKF & ) = default;
            NavKF& operator=(const NavKF &) = default;
            /// \brief Construct from prior state estimate and KF config.
            NavKF(const StateEstEcef & state_est, const KfConfig & kf_config) { state_est_ = state_est; kf_config_ = kf_config;}
            /// \brief Get state estimate.
            inline StateEstEcef getStateEst() const { return state_est_;}
            /// \brief Get time.
            inline T getTime() const  { return state_est_.nav_sol.time;}
            /// \brief Predict (loose).
            inline void lcPredict(const ImuMeasurements & imu_meas, const T & tor_i) { state_est_ = lcPredictKF(state_est_, imu_meas, kf_config_, tor_i);}
            /// \brief Predict (tight).
            inline void tcPredict(const ImuMeasurements & imu_meas, const T & tor_i) { state_est_ = tcPredictKF(state_est_, imu_meas, kf_config_, tor_i);}
            /// \brief Update with position.
            inline void lcUpdatePosEcef(const PosMeasEcef & pos_meas) { state_est_ = lcUpdateKFPosEcef(pos_meas, state_est_);}
            /// \brief Update with position and rotation.
            inline void lcUpdatePosRotEcef(const PosRotMeasEcef & pos_rot_meas) { state_est_ = lcUpdateKFPosRotEcef(pos_rot_meas, state_est_);}
            /// \brief Update with GNSS (loose).
            inline void lcUpdateGnssEcef(const GnssMeasurements & gnss_meas, const GnssConfig & gnss_config) { state_est_ = lcUpdateKFGnssEcef(gnss_meas, state_est_, gnss_config);}
            /// \brief Update with GNSS (tight).
            inline void tcUpdateGnssEcef(const GnssMeasurements & gnss_meas, const T & tor_s) { state_est_ = tcUpdateKFGnssEcef(gnss_meas, state_est_, tor_s);}
            
            #ifdef BUILD_VISION
            /// \brief Update with optical flow
            inline void lcUpdateOptFlow(const cv::Mat & img_0,
                                    const cv::Mat & img_1,
                                    const Matrix3 & K,
                                    const Matrix3 & C_c_b) { 
                                        
                                        state_est_ = lcUpdateKFOptFlow(img_0,
                                                                    img_1,
                                                                    K,
                                                                    C_c_b,
                                                                    state_est_);}
            #endif
    };
};

// Convenience wrappers

using Navigationd = Navigation<double>;
using Navigationf = Navigation<float>;

// double

template<int n_x>
inline void predictKF(const Eigen::Matrix<double, n_x, n_x> & Phi_matrix,
                const Eigen::Matrix<double, n_x, n_x> & Q_matrix,
                const Eigen::Matrix<double, n_x, n_x> & P_matrix_old,
                Eigen::Matrix<double, n_x, n_x> & P_matrix) {
    Navigationd::predictKF<n_x>(Phi_matrix, Q_matrix, P_matrix_old, P_matrix);
}

template<int n_x, int n_z, int max_n_z = 2* Constantsd::kMaxGnssSatellites>
inline bool updateKF(const Eigen::Matrix<double, n_z, 1, 0, max_n_z, 1> & delta_z,
            const Eigen::Matrix<double, n_x, n_x> & P_matrix,
            const Eigen::Matrix<double, n_z, n_x, 0, max_n_z, n_x> & H_matrix,
            const Eigen::Matrix<double, n_z, n_z> & R_matrix,
            const double & p_value,
            Eigen::Matrix<double, n_z, n_z, 0, max_n_z, max_n_z> & S_matrix,
            Eigen::Matrix<double, n_x, 1> & x_est_new,
            Eigen::Matrix<double, n_x, n_x> & P_matrix_post) {
    return Navigationd::updateKF<n_x, n_z, max_n_z>(delta_z, P_matrix, H_matrix, R_matrix, p_value, S_matrix, x_est_new, P_matrix_post);
}

inline Typesd::NavSolutionEcef navEquationsEcef(const Typesd::NavSolutionEcef & old_nav,
                            const Typesd::ImuMeasurements & imu_meas,
                            const double & tor_i) {
    return Navigationd::navEquationsEcef(old_nav, imu_meas, tor_i);
}

inline Typesd::StateEstEcef lcPredictKF(const Typesd::StateEstEcef & state_est_old,
                            const Typesd::ImuMeasurements & imu_meas,
                            const Typesd::KfConfig & kf_config,
                            const double & tor_i) {
    return Navigationd::lcPredictKF(state_est_old, imu_meas, kf_config, tor_i);
}

inline Typesd::StateEstEcef tcPredictKF(const Typesd::StateEstEcef & state_est_old,
                            const Typesd::ImuMeasurements & imu_meas,
                            const Typesd::KfConfig & kf_config,
                            const double & tor_i) {
    return Navigationd::tcPredictKF(state_est_old, imu_meas, kf_config, tor_i);
}

inline Eigen::Matrix<double,15,15> lcPropUnc(const Eigen::Matrix<double,15,15> & P_matrix_old,
                                        const Typesd::NavSolutionEcef & old_nav_est_ecef,
                                        const Typesd::NavSolutionNed & old_nav_est_ned,
                                        const Typesd::ImuMeasurements & imu_meas,
                                        const Typesd::KfConfig & lc_kf_config,
                                        const double & tor_s) {
    return Navigationd::lcPropUnc(P_matrix_old, old_nav_est_ecef, old_nav_est_ned, imu_meas, lc_kf_config, tor_s);
}

inline Eigen::Matrix<double,17,17> tcPropUnc(const Eigen::Matrix<double,17,17> & P_matrix_old,
                                        const Typesd::NavSolutionEcef & old_nav_est_ecef,
                                        const Typesd::NavSolutionNed & old_nav_est_ned,
                                        const Typesd::ImuMeasurements & imu_meas,
                                        const Typesd::KfConfig & tc_kf_config,
                                        const double & tor_s) {
    return Navigationd::tcPropUnc(P_matrix_old, old_nav_est_ecef, old_nav_est_ned, imu_meas, tc_kf_config, tor_s);
}

inline Typesd::StateEstEcef lcUpdateKFPosEcef (const Typesd::PosMeasEcef & pos_meas,
                                    const Typesd::StateEstEcef & state_est_old,
                                    const double & p_value = 0.99) {
    return Navigationd::lcUpdateKFPosEcef(pos_meas, state_est_old, p_value);
}

inline Typesd::StateEstEcef lcUpdateKFGnssEcef (const Typesd::GnssMeasurements & gnss_meas,
                                    const Typesd::StateEstEcef & state_est_old,
                                    const Typesd::GnssConfig & gnss_config,
                                    const double & p_value = 0.99) {
    return Navigationd::lcUpdateKFGnssEcef(gnss_meas, state_est_old, gnss_config, p_value);
}

inline Typesd::StateEstEcef tcUpdateKFGnssEcef (const Typesd::GnssMeasurements & gnss_meas,
                                    const Typesd::StateEstEcef & state_est_prior,
                                    const double & tor_s,
                                    const double & p_value = 0.99) {
    return Navigationd::tcUpdateKFGnssEcef(gnss_meas, state_est_prior, tor_s, p_value);
}

inline Typesd::StateEstEcef lcUpdateKFPosRotEcef (const Typesd::PosRotMeasEcef & pos_rot_meas,
                                    const Typesd::StateEstEcef & state_est_old,
                                    const double & p_value = 0.99) {
    return Navigationd::lcUpdateKFPosRotEcef(pos_rot_meas, state_est_old, p_value);
}

inline Typesd::GnssLsPosVelClock gnssLsPositionVelocityClock(const Typesd::GnssMeasurements & gnss_meas,
                                            const Eigen::Matrix<double,3,1> & prior_r_ea_e,
                                            const Eigen::Matrix<double,3,1> & prior_v_ea_e) {
    return Navigationd::gnssLsPositionVelocityClock(gnss_meas, prior_r_ea_e, prior_v_ea_e);
}

inline Typesd::GnssPosVelMeasEcef gnssLsPositionVelocity(const Typesd::GnssMeasurements & gnss_meas,
                                const Eigen::Matrix<double,3,1> & prior_r_ea_e,
                                const Eigen::Matrix<double,3,1> & prior_v_ea_e,
                                const Typesd::GnssConfig & gnss_config) {
    return Navigationd::gnssLsPositionVelocity(gnss_meas, prior_r_ea_e, prior_v_ea_e, gnss_config);
}

inline Eigen::Matrix<double,17,17> initializePMmatrix(const Typesd::KfConfig & tc_kf_config) {
    return Navigationd::initializePMmatrix(tc_kf_config);
}

inline Typesd::StateEstEcef initStateFromGroundTruth(const Typesd::NavSolutionEcef & true_nav_ecef, 
                                            const Typesd::KfConfig & kf_config, 
                                            const Typesd::GnssMeasurements & gnss_meas, 
                                            std::mt19937 & gen) {
    return Navigationd::initStateFromGroundTruth(true_nav_ecef, kf_config, gnss_meas, gen);
}

// float 

template<int n_x>
inline void predictKF(const Eigen::Matrix<float, n_x, n_x> & Phi_matrix,
                const Eigen::Matrix<float, n_x, n_x> & Q_matrix,
                const Eigen::Matrix<float, n_x, n_x> & P_matrix_old,
                Eigen::Matrix<float, n_x, n_x> & P_matrix) {
    Navigationf::predictKF<n_x>(Phi_matrix, Q_matrix, P_matrix_old, P_matrix);
}

template<int n_x, int n_z, int max_n_z = 2* Constantsf::kMaxGnssSatellites>
inline bool updateKF(const Eigen::Matrix<float, n_z, 1, 0, max_n_z, 1> & delta_z,
            const Eigen::Matrix<float, n_x, n_x> & P_matrix,
            const Eigen::Matrix<float, n_z, n_x, 0, max_n_z, n_x> & H_matrix,
            const Eigen::Matrix<float, n_z, n_z> & R_matrix,
            const float & p_value,
            Eigen::Matrix<float, n_z, n_z, 0, max_n_z, max_n_z> & S_matrix,
            Eigen::Matrix<float, n_x, 1> & x_est_new,
            Eigen::Matrix<float, n_x, n_x> & P_matrix_post) {
    return Navigationf::updateKF<n_x, n_z, max_n_z>(delta_z, P_matrix, H_matrix, R_matrix, p_value, S_matrix, x_est_new, P_matrix_post);
}

inline Typesf::NavSolutionEcef navEquationsEcef(const Typesf::NavSolutionEcef & old_nav,
                            const Typesf::ImuMeasurements & imu_meas,
                            const float & tor_i) {
    return Navigationf::navEquationsEcef(old_nav, imu_meas, tor_i);
}

inline Typesf::StateEstEcef lcPredictKF(const Typesf::StateEstEcef & state_est_old,
                            const Typesf::ImuMeasurements & imu_meas,
                            const Typesf::KfConfig & kf_config,
                            const float & tor_i) {
    return Navigationf::lcPredictKF(state_est_old, imu_meas, kf_config, tor_i);
}

inline Typesf::StateEstEcef tcPredictKF(const Typesf::StateEstEcef & state_est_old,
                            const Typesf::ImuMeasurements & imu_meas,
                            const Typesf::KfConfig & kf_config,
                            const float & tor_i) {
    return Navigationf::tcPredictKF(state_est_old, imu_meas, kf_config, tor_i);
}

inline Eigen::Matrix<float,15,15> lcPropUnc(const Eigen::Matrix<float,15,15> & P_matrix_old,
                                        const Typesf::NavSolutionEcef & old_nav_est_ecef,
                                        const Typesf::NavSolutionNed & old_nav_est_ned,
                                        const Typesf::ImuMeasurements & imu_meas,
                                        const Typesf::KfConfig & lc_kf_config,
                                        const float & tor_s) {
    return Navigationf::lcPropUnc(P_matrix_old, old_nav_est_ecef, old_nav_est_ned, imu_meas, lc_kf_config, tor_s);
}

inline Eigen::Matrix<float,17,17> tcPropUnc(const Eigen::Matrix<float,17,17> & P_matrix_old,
                                        const Typesf::NavSolutionEcef & old_nav_est_ecef,
                                        const Typesf::NavSolutionNed & old_nav_est_ned,
                                        const Typesf::ImuMeasurements & imu_meas,
                                        const Typesf::KfConfig & tc_kf_config,
                                        const float & tor_s) {
    return Navigationf::tcPropUnc(P_matrix_old, old_nav_est_ecef, old_nav_est_ned, imu_meas, tc_kf_config, tor_s);
}

inline Typesf::StateEstEcef lcUpdateKFPosEcef (const Typesf::PosMeasEcef & pos_meas,
                                    const Typesf::StateEstEcef & state_est_old,
                                    const float & p_value = 0.99) {
    return Navigationf::lcUpdateKFPosEcef(pos_meas, state_est_old, p_value);
}

inline Typesf::StateEstEcef lcUpdateKFGnssEcef (const Typesf::GnssMeasurements & gnss_meas,
                                    const Typesf::StateEstEcef & state_est_old,
                                    const Typesf::GnssConfig & gnss_config,
                                    const float & p_value = 0.99) {
    return Navigationf::lcUpdateKFGnssEcef(gnss_meas, state_est_old, gnss_config, p_value);
}

inline Typesf::StateEstEcef tcUpdateKFGnssEcef (const Typesf::GnssMeasurements & gnss_meas,
                                    const Typesf::StateEstEcef & state_est_prior,
                                    const float & tor_s,
                                    const float & p_value = 0.99) {
    return Navigationf::tcUpdateKFGnssEcef(gnss_meas, state_est_prior, tor_s, p_value);
}

inline Typesf::StateEstEcef lcUpdateKFPosRotEcef (const Typesf::PosRotMeasEcef & pos_rot_meas,
                                    const Typesf::StateEstEcef & state_est_old,
                                    const float & p_value = 0.99) {
    return Navigationf::lcUpdateKFPosRotEcef(pos_rot_meas, state_est_old, p_value);
}

inline Typesf::GnssLsPosVelClock gnssLsPositionVelocityClock(const Typesf::GnssMeasurements & gnss_meas,
                                            const Eigen::Matrix<float,3,1> & prior_r_ea_e,
                                            const Eigen::Matrix<float,3,1> & prior_v_ea_e) {
    return Navigationf::gnssLsPositionVelocityClock(gnss_meas, prior_r_ea_e, prior_v_ea_e);
}

inline Typesf::GnssPosVelMeasEcef gnssLsPositionVelocity(const Typesf::GnssMeasurements & gnss_meas,
                                const Eigen::Matrix<float,3,1> & prior_r_ea_e,
                                const Eigen::Matrix<float,3,1> & prior_v_ea_e,
                                const Typesf::GnssConfig & gnss_config) {
    return Navigationf::gnssLsPositionVelocity(gnss_meas, prior_r_ea_e, prior_v_ea_e, gnss_config);
}

inline Eigen::Matrix<float,17,17> initializePMmatrix(const Typesf::KfConfig & tc_kf_config) {
    return Navigationf::initializePMmatrix(tc_kf_config);
}

inline Typesf::StateEstEcef initStateFromGroundTruth(const Typesf::NavSolutionEcef & true_nav_ecef, 
                                            const Typesf::KfConfig & kf_config, 
                                            const Typesf::GnssMeasurements & gnss_meas, 
                                            std::mt19937 & gen) {
    return Navigationf::initStateFromGroundTruth(true_nav_ecef, kf_config, gnss_meas, gen);
}

extern template struct Navigation<double>;
extern template struct Navigation<float>;

};

/// @}

#endif // NAVIGATION_H
