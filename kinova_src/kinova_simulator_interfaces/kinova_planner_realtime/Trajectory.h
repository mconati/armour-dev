#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include "PZsparse.h"

// These values are specifically corresponded with a Bezier curve parameterization
#define QDD_DES_K_DEP_MAXIMA (0.5 - sqrt(3) / 6)
#define QDD_DES_K_DEP_MINIMA (0.5 + sqrt(3) / 6)

// 5th order Bezier curve
// The initial position/velocity/acceleration is equal to q0/qd0/qdd0
// The end position is equal to q0 + k
// The end velocity/acceleration is equal to 0
//
// NOTE:
// This is just a simplified implementation!!!
// t is automatically set to range from 0 to 1, so you don't have to scale.
// Everything has to be changed if the range of t is not [0,1]
//
// -1 <= k <= 1
// k_actual = k * k_range
//
// q_des   = t^3*(6*t^2 - 15*t + 10) * k_actual + 
//           q0 + qd0*t - 6*qd0*t^3 + 8*qd0*t^4 - 3*qd0*t^5 + (qdd0*t^2)/2 - (3*qdd0*t^3)/2 + (3*qdd0*t^4)/2 - (qdd0*t^5)/2
//
// qd_des  = 30*t^2*(t - 1)^2 * k_actual + 
//           ((t - 1)^2*(2*qd0 + 4*qd0*t + 2*qdd0*t - 30*qd0*t^2 - 5*qdd0*t^2))/2
//
// qdd_des = 60*t*(2*t^2 - 3*t + 1) * k_actual + 
//           -(t - 1)*(qdd0 - 36*qd0*t - 8*qdd0*t + 60*qd0*t^2 + 10*qdd0*t^2)
//
class BezierCurve{
public:
    TYPE* q0 = nullptr;
    TYPE* qd0 = nullptr;
    TYPE* qdd0 = nullptr;

    TYPE q_des_k_indep_extrema_1[NUM_FACTORS] = {0.0};
    TYPE q_des_k_indep_extrema_2[NUM_FACTORS] = {0.0};
    TYPE q_des_k_indep_extremum_1[NUM_FACTORS] = {0.0};
    TYPE q_des_k_indep_extremum_2[NUM_FACTORS] = {0.0};

    TYPE qd_des_k_indep_extrema_1[NUM_FACTORS] = {0.0};
    TYPE qd_des_k_indep_extrema_2[NUM_FACTORS] = {0.0};
    TYPE qd_des_k_indep_extremum_1[NUM_FACTORS] = {0.0};
    TYPE qd_des_k_indep_extremum_2[NUM_FACTORS] = {0.0};

    TYPE qdd_des_k_indep_extrema_1[NUM_FACTORS] = {0.0};
    TYPE qdd_des_k_indep_extrema_2[NUM_FACTORS] = {0.0};
    TYPE qdd_des_k_indep_extremum_1[NUM_FACTORS] = {0.0};
    TYPE qdd_des_k_indep_extremum_2[NUM_FACTORS] = {0.0};

    TYPE dt;

    PZsparse cos_q_des[NUM_TIME_STEPS * NUM_FACTORS];
    PZsparse sin_q_des[NUM_TIME_STEPS * NUM_FACTORS];
    PZsparse qd_des[NUM_TIME_STEPS * NUM_FACTORS];
    PZsparse qda_des[NUM_TIME_STEPS * NUM_FACTORS];
    PZsparse qdda_des[NUM_TIME_STEPS * NUM_FACTORS];

    BezierCurve() {};

    BezierCurve(TYPE* q0_inp, TYPE* qd0_inp, TYPE* qdd0_inp);

    ~BezierCurve() {};

    // convert to polynomial zonotope using 1st/2nd order Taylor expansion
    void makePolyZono(int t_ind);

    // return the min and max of the joint position throughout the whole desired trajectory
    void returnJointPositionExtremum(TYPE* extremum, const TYPE* k);

    // return the gradient of min and max of the joint position throughout the whole desired trajectory
    void returnJointPositionExtremumGradient(TYPE* extremumGradient, const TYPE* k);

    // return the min and max of the joint velocity throughout the whole desired trajectory
    void returnJointVelocityExtremum(TYPE* extremum, const TYPE* k);

    // return the gradient of min and max of the joint velocity throughout the whole desired trajectory
    void returnJointVelocityExtremumGradient(TYPE* extremumGradient, const TYPE* k);
};

// helper functions
// q0, qd0, qdd0, k here are scalars since all joints are using the same Bezier curve representation
TYPE q_des_func(TYPE q0, TYPE qd0, TYPE qdd0, TYPE k, TYPE t);

TYPE qd_des_func(TYPE q0, TYPE qd0, TYPE qdd0, TYPE k, TYPE t);

TYPE qdd_des_func(TYPE q0, TYPE qd0, TYPE qdd0, TYPE k, TYPE t);

// derivative of the second extrema of q_des (when qd_des = 0) w.r.t k 
TYPE q_des_extrema2_k_derivative(TYPE q0, TYPE qd0, TYPE qdd0, TYPE k);

// derivative of the third extrema of q_des (when qd_des = 0) w.r.t k 
TYPE q_des_extrema3_k_derivative(TYPE q0, TYPE qd0, TYPE qdd0, TYPE k);

// derivative of the second extrema of qd_des (when qd_des = 0) w.r.t k 
TYPE qd_des_extrema2_k_derivative(TYPE q0, TYPE qd0, TYPE qdd0, TYPE k);

// derivative of the third extrema of qd_des (when qd_des = 0) w.r.t k 
TYPE qd_des_extrema3_k_derivative(TYPE q0, TYPE qd0, TYPE qdd0, TYPE k);

// k-independent part of q_des
TYPE q_des_k_indep(TYPE q0, TYPE qd0, TYPE qdd0, TYPE t);

// k-independent part of qd_des
TYPE qd_des_k_indep(TYPE q0, TYPE qd0, TYPE qdd0, TYPE t);

// k-independent part of qdd_des
TYPE qdd_des_k_indep(TYPE q0, TYPE qd0, TYPE qdd0, TYPE t);

#endif