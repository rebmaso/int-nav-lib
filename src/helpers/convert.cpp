#include "helpers.h"

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

template struct Helpers<float>;
template struct Helpers<double>;

};