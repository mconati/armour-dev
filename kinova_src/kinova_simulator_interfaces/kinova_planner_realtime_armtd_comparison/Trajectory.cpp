#ifndef TRAJECTORY_CPP
#define TRAJECTORY_CPP

#include "Trajectory.h"

ConstantAccelerationCurve::ConstantAccelerationCurve(TYPE* q0_inp, TYPE* qd0_inp, 
                                                     TYPE* c_cos_q_des_inp, TYPE* g_cos_q_des_inp, TYPE* r_cos_q_des_inp,
                                                     TYPE* c_sin_q_des_inp, TYPE* g_sin_q_des_inp, TYPE* r_sin_q_des_inp,
                                                     TYPE* k_range_inp) {
    q0 = q0_inp;
    qd0 = qd0_inp; 

    c_cos_q_des = c_cos_q_des_inp;
    g_cos_q_des = g_cos_q_des_inp;
    r_cos_q_des = r_cos_q_des_inp;

    c_sin_q_des = c_sin_q_des_inp;
    g_sin_q_des = g_sin_q_des_inp;
    r_sin_q_des = r_sin_q_des_inp;

    k_range = k_range_inp;

    dt = 1.0 / NUM_TIME_STEPS;
}

void ConstantAccelerationCurve::makePolyZono(int t_ind) {
    assert(t_ind < NUM_TIME_STEPS);

    for (int i = 0; i < NUM_FACTORS; i++) {
        const TYPE k_range_elt = k_range[i];
        const TYPE cos_q0 = cos(q0[i]);
        const TYPE sin_q0 = sin(q0[i]);

        // cos(q_des)
        TYPE cos_q_des_center = cos_q0 * c_cos_q_des[i * NUM_TIME_STEPS + t_ind] - sin_q0 * c_sin_q_des[i * NUM_TIME_STEPS + t_ind];
        TYPE cos_q_des_coeff[2];
        cos_q_des_coeff[0] = cos_q0 * g_cos_q_des[i * NUM_TIME_STEPS + t_ind] - sin_q0 * g_sin_q_des[i * NUM_TIME_STEPS + t_ind];
        cos_q_des_coeff[1] = fabs(cos_q0) * r_cos_q_des[i * NUM_TIME_STEPS + t_ind] + fabs(sin_q0) * r_sin_q_des[i * NUM_TIME_STEPS + t_ind];

        cos_q_des_coeff[1] *= 5.0;

        unsigned int cos_q_des_degree[2][NUM_FACTORS * 3] = {0};
        cos_q_des_degree[0][i] = 1; // k
        cos_q_des_degree[1][i + NUM_FACTORS * 1] = 1; // cosqe
        cos_q_des[t_ind * NUM_FACTORS + i] = PZsparse(cos_q_des_center, cos_q_des_coeff, cos_q_des_degree, 2);

        // sin(q_des)
        TYPE sin_q_des_center = cos_q0 * c_sin_q_des[i * NUM_TIME_STEPS + t_ind] + sin_q0 * c_cos_q_des[i * NUM_TIME_STEPS + t_ind];
        TYPE sin_q_des_coeff[2];
        sin_q_des_coeff[0] = cos_q0 * g_sin_q_des[i * NUM_TIME_STEPS + t_ind] + sin_q0 * g_cos_q_des[i * NUM_TIME_STEPS + t_ind];
        sin_q_des_coeff[1] = fabs(cos_q0) * r_sin_q_des[i * NUM_TIME_STEPS + t_ind] + fabs(sin_q0) * r_cos_q_des[i * NUM_TIME_STEPS + t_ind];

        sin_q_des_coeff[1] *= 5.0;

        unsigned int sin_q_des_degree[2][NUM_FACTORS * 3] = {0};
        sin_q_des_degree[0][i] = 1; // k
        sin_q_des_degree[1][i + NUM_FACTORS * 2] = 1; // sinqe
        sin_q_des[t_ind * NUM_FACTORS + i] = PZsparse(sin_q_des_center, sin_q_des_coeff, sin_q_des_degree, 2);
    }
}

void ConstantAccelerationCurve::returnJointPositionExtremum(TYPE* extremum, const TYPE* k) {
    for (int i = 0; i < NUM_FACTORS; i++) {
        // k[i] range is [-1,1] since it is defined for PZ, the following is the actual k
        TYPE k_actual = k_range[i] * k[i];

        extremum[i] = 0;
        extremum[i + NUM_FACTORS] = 0;
    }
}

void ConstantAccelerationCurve::returnJointPositionExtremumGradient(TYPE* extremumGradient, const TYPE* k) {
    for (int i = 0; i < NUM_FACTORS; i++) {
        // k[i] range is [-1,1] since it is defined for PZ, the following is the actual k
        TYPE k_actual = k_range[i] * k[i];

        for (int j = 0; j < NUM_FACTORS; j++) {
            extremumGradient[(i              ) * NUM_FACTORS + j] = 0.0;
            extremumGradient[(i + NUM_FACTORS) * NUM_FACTORS + j] = 0.0;
        }
    }
}

void ConstantAccelerationCurve::returnJointVelocityExtremum(TYPE* extremum, const TYPE* k) {
    for (int i = 0; i < NUM_FACTORS; i++) {
        // k[i] range is [-1,1] since it is defined for PZ, the following is the actual k
        TYPE k_actual = k_range[i] * k[i];

        extremum[i] = 0;
        extremum[i + NUM_FACTORS] = 0;
    }
}

void ConstantAccelerationCurve::returnJointVelocityExtremumGradient(TYPE* extremumGradient, const TYPE* k) {
    for (int i = 0; i < NUM_FACTORS; i++) {
        // k[i] range is [-1,1] since it is defined for PZ, the following is the actual k
        TYPE k_actual = k_range[i] * k[i];

        for (int j = 0; j < NUM_FACTORS; j++) {
            extremumGradient[(i              ) * NUM_FACTORS + j] = 0.0;
            extremumGradient[(i + NUM_FACTORS) * NUM_FACTORS + j] = 0.0;
        }
    }
}

#endif
