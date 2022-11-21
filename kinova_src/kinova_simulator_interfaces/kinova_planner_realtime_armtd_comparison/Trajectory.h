#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include "PZsparse.h"

// A simple quadratic trajectory whose acceleration is parameterized
//
// -1 <= k <= 1
// k_actual = k * k_range
//
// q_des   = q0 + qd0 * t + 0.5 * k_actual * t^2
//
// qd_des  = qd0 + k_actual * t
//
// qdd_des = k_actual
//
class ConstantAccelerationCurve{
public:
    TYPE* q0 = nullptr;
    TYPE* qd0 = nullptr;

    TYPE dt;

    TYPE* c_cos_q_des = nullptr; // center of zonotope
    TYPE* g_cos_q_des = nullptr; // k-dependent generator of zonotope
    TYPE* r_cos_q_des = nullptr; // radius of zonotope

    TYPE* c_sin_q_des = nullptr; // center of zonotope
    TYPE* g_sin_q_des = nullptr; // k-dependent generator of zonotope
    TYPE* r_sin_q_des = nullptr; // radius of zonotope

    PZsparse cos_q_des[NUM_TIME_STEPS * NUM_FACTORS];
    PZsparse sin_q_des[NUM_TIME_STEPS * NUM_FACTORS];

    TYPE* k_range = nullptr;

    ConstantAccelerationCurve() {};

    ConstantAccelerationCurve(TYPE* q0_inp, TYPE* qd0_inp, 
                              TYPE* c_cos_q_des_inp, TYPE* g_cos_q_des_inp, TYPE* r_cos_q_des_inp,
                              TYPE* c_sin_q_des_inp, TYPE* g_sin_q_des_inp, TYPE* r_sin_q_des_inp,
                              TYPE* k_range_inp);

    ~ConstantAccelerationCurve() {};

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

#endif