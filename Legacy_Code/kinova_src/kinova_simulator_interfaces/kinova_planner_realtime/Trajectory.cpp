#ifndef TRAJECTORY_CPP
#define TRAJECTORY_CPP

#include "Trajectory.h"

BezierCurve::BezierCurve(TYPE* q0_inp, TYPE* qd0_inp, TYPE* qdd0_inp) {
    q0 = q0_inp;
    qd0 = qd0_inp;
    qdd0 = qdd0_inp;    

    // initialize the extrema of the k independent part of q_des
    for (int i = 0; i < NUM_FACTORS; i++) {
        q_des_k_indep_extrema_1[i] = (2*qd0[i] + qdd0[i] + sqrt(64*pow(qd0[i],2) + 14*qd0[i]*qdd0[i] + pow(qdd0[i],2)))/(5*(6*qd0[i] + qdd0[i]));
        q_des_k_indep_extrema_2[i] = (2*qd0[i] + qdd0[i] - sqrt(64*pow(qd0[i],2) + 14*qd0[i]*qdd0[i] + pow(qdd0[i],2)))/(5*(6*qd0[i] + qdd0[i]));
        q_des_k_indep_extremum_1[i] = q_des_k_indep(q0[i], qd0[i], qdd0[i], q_des_k_indep_extrema_1[i]);
        q_des_k_indep_extremum_2[i] = q_des_k_indep(q0[i], qd0[i], qdd0[i], q_des_k_indep_extrema_2[i]);
    }

    // initialize the extrema of the k independent part of qd_des
    for (int i = 0; i < NUM_FACTORS; i++) {
        qd_des_k_indep_extrema_1[i] = (18*qd0[i] + 4*qdd0[i] + sqrt(6*(54*pow(qd0[i],2) + 14*qd0[i]*qdd0[i] + pow(qdd0[i],2))))/(10*(6*qd0[i] + qdd0[i]));
        qd_des_k_indep_extrema_2[i] = (18*qd0[i] + 4*qdd0[i] - sqrt(6*(54*pow(qd0[i],2) + 14*qd0[i]*qdd0[i] + pow(qdd0[i],2))))/(10*(6*qd0[i] + qdd0[i]));
        qd_des_k_indep_extremum_1[i] = qd_des_k_indep(q0[i], qd0[i], qdd0[i], qd_des_k_indep_extrema_1[i]);
        qd_des_k_indep_extremum_2[i] = qd_des_k_indep(q0[i], qd0[i], qdd0[i], qd_des_k_indep_extrema_2[i]);
    }

    // initialize the extrema of the k independent part of qdd_des
    for (int i = 0; i < NUM_FACTORS; i++) {
        qdd_des_k_indep_extrema_1[i] = (32*qd0[i] + 6*qdd0[i] + sqrt(2*(152*pow(qd0[i],2) + 42*qd0[i]*qdd0[i] + 3*pow(qdd0[i],2))))/(10*(6*qd0[i] + qdd0[i]));
        qdd_des_k_indep_extrema_2[i] = (32*qd0[i] + 6*qdd0[i] - sqrt(2*(152*pow(qd0[i],2) + 42*qd0[i]*qdd0[i] + 3*pow(qdd0[i],2))))/(10*(6*qd0[i] + qdd0[i]));
        qdd_des_k_indep_extremum_1[i] = qdd_des_k_indep(q0[i], qd0[i], qdd0[i], qdd_des_k_indep_extrema_1[i]);
        qdd_des_k_indep_extremum_2[i] = qdd_des_k_indep(q0[i], qd0[i], qdd0[i], qdd_des_k_indep_extrema_2[i]);
    }

    dt = 1.0 / NUM_TIME_STEPS;
}

void BezierCurve::makePolyZono(int t_ind) {
    assert(t_ind < NUM_TIME_STEPS);

    const TYPE t_lb = t_ind * dt;
    const TYPE t_ub = (t_ind + 1) * dt;

    const Interval t_int(t_lb, t_ub);

    for (int i = 0; i < NUM_FACTORS; i++) {
        const TYPE k_range_elt = k_range[i];

        // Part 1: q_des
        TYPE k_dep_coeff_lb = pow(t_lb,3) * (6 * pow(t_lb,2) - 15 * t_lb + 10);
        TYPE k_dep_coeff_ub = pow(t_ub,3) * (6 * pow(t_ub,2) - 15 * t_ub + 10);
        TYPE k_dep_coeff_center = (k_dep_coeff_ub + k_dep_coeff_lb) * 0.5;
        TYPE k_dep_coeff_radius = (k_dep_coeff_ub - k_dep_coeff_lb) * 0.5 * k_range_elt;
        
        TYPE k_indep_lb = q_des_k_indep(q0[i], qd0[i], qdd0[i], t_lb);
        TYPE k_indep_ub = q_des_k_indep(q0[i], qd0[i], qdd0[i], t_ub);
        if (k_indep_lb > k_indep_ub) {
            swap(k_indep_lb, k_indep_ub);
        }
        if (t_lb < q_des_k_indep_extrema_1[i] && q_des_k_indep_extrema_1[i] < t_ub) {
            k_indep_lb = min(k_indep_lb, q_des_k_indep_extremum_1[i]);
            k_indep_ub = max(k_indep_ub, q_des_k_indep_extremum_1[i]);
        }
        if (t_lb < q_des_k_indep_extrema_2[i] && q_des_k_indep_extrema_2[i] < t_ub) {
            k_indep_lb = min(k_indep_lb, q_des_k_indep_extremum_2[i]);
            k_indep_ub = max(k_indep_ub, q_des_k_indep_extremum_2[i]);
        }
        TYPE k_indep_radius = (k_indep_ub - k_indep_lb) * 0.5;
        TYPE q_des_center = (k_indep_lb + k_indep_ub) * 0.5;
        
        // q_des_k_dep = k_dep_coeff_center * k;
        Interval q_des_radius_int(-k_dep_coeff_radius - k_indep_radius - qe, k_dep_coeff_radius + k_indep_radius + qe);
        
        // q_des_int = q_des_center + q_des_k_dep + q_des_radius_int;

        // first order Taylor expansion
        // Part 1.a: cos(q_des) 
        TYPE cos_q_des_center = cos(q_des_center);
        Interval cos_q_des_radius_int = - q_des_radius_int * sin(q_des_center) 
                                        - 0.5 * cos(q_des_center + k_dep_coeff_center * Interval(-k_range_elt, k_range_elt) + q_des_radius_int) 
                                            * pow(q_des_radius_int + k_dep_coeff_center * Interval(-k_range_elt, k_range_elt), 2);

        cos_q_des_center += getCenter(cos_q_des_radius_int);
        cos_q_des_radius_int = cos_q_des_radius_int - getCenter(cos_q_des_radius_int);
        TYPE cos_q_des_coeff[] = {-k_dep_coeff_center * k_range_elt * sin(q_des_center), getRadius(cos_q_des_radius_int)}; 

        // cos_q_des_int = cos_q_des_center + cos_q_des_coeff[0] * k + cos_q_des_coeff[1] * cosqe;
        uint64_t cos_q_des_degree[2][NUM_FACTORS * 6] = {0};
        cos_q_des_degree[0][i] = 1; // k
        cos_q_des_degree[1][i + NUM_FACTORS * 1] = 1; // cosqe

        cos_q_des[t_ind * NUM_FACTORS + i] = PZsparse(cos_q_des_center, cos_q_des_coeff, cos_q_des_degree, 2);

        // Part 1.b: sin(q_des) 
        TYPE sin_q_des_center = sin(q_des_center);
        Interval sin_q_des_radius_int = q_des_radius_int * cos(q_des_center) 
                                        - 0.5 * sin(q_des_center + k_dep_coeff_center * Interval(-k_range_elt, k_range_elt) + q_des_radius_int) 
                                            * pow(q_des_radius_int + k_dep_coeff_center * Interval(-k_range_elt, k_range_elt), 2);

        sin_q_des_center += getCenter(sin_q_des_radius_int);
        sin_q_des_radius_int = sin_q_des_radius_int - getCenter(sin_q_des_radius_int);
        TYPE sin_q_des_coeff[] = {k_dep_coeff_center * k_range_elt * cos(q_des_center), getRadius(sin_q_des_radius_int)};

        uint64_t sin_q_des_degree[2][NUM_FACTORS * 6] = {0};
        sin_q_des_degree[0][i] = 1; // k
        sin_q_des_degree[1][i + NUM_FACTORS * 2] = 1; // sinqe

        // sin_q_des_int = sin_q_des_center + sin_q_des_coeff[0] * k + sin_q_des_coeff[1] * sinqe;
        sin_q_des[t_ind * NUM_FACTORS + i] = PZsparse(sin_q_des_center, sin_q_des_coeff, sin_q_des_degree, 2);
        
        // Part 2: qd_des
        // NOTE:
        // This is just a simplified implementation!!!
        // 30*t^2*(t - 1)^2 in qd_des is just a function with one maxima at t = 0.5
        // So as long as NUM_TIME_STEPS is even number, the following bounding trick holds!
        k_dep_coeff_lb = 30 * pow(t_lb,2) * pow(t_lb - 1,2);
        k_dep_coeff_ub = 30 * pow(t_ub,2) * pow(t_ub - 1,2);
        if (k_dep_coeff_ub < k_dep_coeff_lb) { // we are at t >= 0.5, which is a monotonically decreasing region 
            swap(k_dep_coeff_lb, k_dep_coeff_ub);
        }

        k_dep_coeff_center = (k_dep_coeff_ub + k_dep_coeff_lb) * 0.5 * k_range_elt; // Have to scale to [-1,1] in order to fit in PZsparse
        k_dep_coeff_radius = (k_dep_coeff_ub - k_dep_coeff_lb) * 0.5 * k_range_elt;

        k_indep_lb = qd_des_k_indep(q0[i], qd0[i], qdd0[i], t_lb);
        k_indep_ub = qd_des_k_indep(q0[i], qd0[i], qdd0[i], t_ub);
        if (k_indep_lb > k_indep_ub) {
            swap(k_indep_lb, k_indep_ub);
        }
        if (t_lb < qd_des_k_indep_extrema_1[i] && qd_des_k_indep_extrema_1[i] < t_ub) {
            k_indep_lb = min(k_indep_lb, qd_des_k_indep_extremum_1[i]);
            k_indep_ub = max(k_indep_ub, qd_des_k_indep_extremum_1[i]);
        }
        if (t_lb < qd_des_k_indep_extrema_2[i] && qd_des_k_indep_extrema_2[i] < t_ub) {
            k_indep_lb = min(k_indep_lb, qd_des_k_indep_extremum_2[i]);
            k_indep_ub = max(k_indep_ub, qd_des_k_indep_extremum_2[i]);
        }
        k_indep_radius = (k_indep_ub - k_indep_lb) * 0.5;
        TYPE qd_des_center = (k_indep_lb + k_indep_ub) * 0.5;

        TYPE qd_des_coeff[] = {k_dep_coeff_center, k_dep_coeff_radius + k_indep_radius + qde};

        uint64_t qd_des_degree[2][NUM_FACTORS * 6] = {0};
        qd_des_degree[0][i] = 1; // k
        qd_des_degree[1][i + NUM_FACTORS * 3] = 1; // qde

        // qd_des_int = qd_des_center + qd_des_coeff[0] * k + qd_des_coeff[1] * qde;
        qd_des[t_ind * NUM_FACTORS + i] = PZsparse(qd_des_center, qd_des_coeff, qd_des_degree, 2);

        TYPE qda_des_coeff[] = {k_dep_coeff_center, k_dep_coeff_radius + k_indep_radius + qdae};

        uint64_t qda_des_degree[2][NUM_FACTORS * 6] = {0};
        qda_des_degree[0][i] = 1; // k
        qda_des_degree[1][i + NUM_FACTORS * 4] = 1; // qdae

        // qda_des_int = qd_des_center + qda_des_coeff[0] * k + qda_des_coeff[1] * qdae;
        qda_des[t_ind * NUM_FACTORS + i] = PZsparse(qd_des_center, qda_des_coeff, qda_des_degree, 2);

        // Part 3: qdd_des
        TYPE temp_lb = 60 * t_lb * (2 * pow(t_lb,2) - 3 * t_lb + 1);
        TYPE temp_ub = 60 * t_ub * (2 * pow(t_ub,2) - 3 * t_ub + 1);
        if (t_ub <= QDD_DES_K_DEP_MAXIMA) { // monotonically increasing region
            k_dep_coeff_lb = temp_lb;
            k_dep_coeff_ub = temp_ub;
        }
        else if (t_lb <= QDD_DES_K_DEP_MAXIMA) { // maxima lives inside
            k_dep_coeff_lb = min(temp_lb, temp_ub);
            k_dep_coeff_ub = 60 * QDD_DES_K_DEP_MAXIMA * (2 * pow(QDD_DES_K_DEP_MAXIMA,2) - 3 * QDD_DES_K_DEP_MAXIMA + 1);
        }
        else if (t_ub <= QDD_DES_K_DEP_MINIMA) { // monotonically decreasing region
            k_dep_coeff_lb = temp_ub;
            k_dep_coeff_ub = temp_lb;
        }
        else if (t_lb <= QDD_DES_K_DEP_MINIMA) { // minima lives inside
            k_dep_coeff_lb = 60 * QDD_DES_K_DEP_MINIMA * (2 * pow(QDD_DES_K_DEP_MINIMA,2) - 3 * QDD_DES_K_DEP_MINIMA + 1);
            k_dep_coeff_ub = max(temp_lb, temp_ub);
        }
        else { // monotonically increasing region
            k_dep_coeff_lb = temp_lb;
            k_dep_coeff_ub = temp_ub;
        }
        
        k_dep_coeff_center = (k_dep_coeff_ub + k_dep_coeff_lb) * 0.5 * k_range_elt; // Have to scale to [-1,1] in order to fit in PZsparse
        k_dep_coeff_radius = (k_dep_coeff_ub - k_dep_coeff_lb) * 0.5 * k_range_elt;

        k_indep_lb = qdd_des_k_indep(q0[i], qd0[i], qdd0[i], t_lb);
        k_indep_ub = qdd_des_k_indep(q0[i], qd0[i], qdd0[i], t_ub);
        if (k_indep_lb > k_indep_ub) {
            swap(k_indep_lb, k_indep_ub);
        }
        if (t_lb < qdd_des_k_indep_extrema_1[i] && qdd_des_k_indep_extrema_1[i] < t_ub) {
            k_indep_lb = min(k_indep_lb, qdd_des_k_indep_extremum_1[i]);
            k_indep_ub = max(k_indep_ub, qdd_des_k_indep_extremum_1[i]);
        }
        if (t_lb < qdd_des_k_indep_extrema_2[i] && qdd_des_k_indep_extrema_2[i] < t_ub) {
            k_indep_lb = min(k_indep_lb, qdd_des_k_indep_extremum_2[i]);
            k_indep_ub = max(k_indep_ub, qdd_des_k_indep_extremum_2[i]);
        }
        k_indep_radius = (k_indep_ub - k_indep_lb) * 0.5;
        TYPE qdd_des_center = (k_indep_lb + k_indep_ub) * 0.5;

        TYPE qdd_des_coeff[] = {k_dep_coeff_center, k_dep_coeff_radius + k_indep_radius + qddae};

        uint64_t qdd_des_degree[2][NUM_FACTORS * 6] = {0};
        qdd_des_degree[0][i] = 1; // k
        qdd_des_degree[1][i + NUM_FACTORS * 5] = 1; // qadae

        // qdd_des_int = qdd_des_center + qdd_des_coeff[0] * k + qdd_des_coeff[1] * qdde;
        qdda_des[t_ind * NUM_FACTORS + i] = PZsparse(qdd_des_center, qdd_des_coeff, qdd_des_degree, 2);
    }
}

void BezierCurve::returnJointPositionExtremum(TYPE* extremum, const TYPE* k) {
    for (int i = 0; i < NUM_FACTORS; i++) {
        // k[i] range is [-1,1] since it is defined for PZ, the following is the actual k
        TYPE k_actual = k_range[i] * k[i];

        // list all possible extremas
        TYPE extrema1 = 0;
        TYPE extrema2 = (2*qd0[i] + qdd0[i] + sqrt(64*pow(qd0[i],2) + 14*qd0[i]*qdd0[i] - 120*k_actual*qd0[i] + pow(qdd0[i],2))) / (5*(6*qd0[i] - 12*k_actual + qdd0[i]));
        TYPE extrema3 = (2*qd0[i] + qdd0[i] - sqrt(64*pow(qd0[i],2) + 14*qd0[i]*qdd0[i] - 120*k_actual*qd0[i] + pow(qdd0[i],2))) / (5*(6*qd0[i] - 12*k_actual + qdd0[i]));
        TYPE extrema4 = 1;

        // get extremums of all extremas
        TYPE extremum1 = q_des_func(q0[i], qd0[i], qdd0[i], k_actual, extrema1);
        TYPE extremum2 = q_des_func(q0[i], qd0[i], qdd0[i], k_actual, extrema2);
        TYPE extremum3 = q_des_func(q0[i], qd0[i], qdd0[i], k_actual, extrema3);
        TYPE extremum4 = q_des_func(q0[i], qd0[i], qdd0[i], k_actual, extrema4);

        // find the min and max values
        TYPE minPosition = min(extremum1, extremum4);
        TYPE maxPosition = max(extremum1, extremum4);
        if (0 <= extrema2 && extrema2 <= 1) { // check if this extrema is inside the time range [0,1]
            minPosition = min(minPosition, extremum2);
            maxPosition = max(maxPosition, extremum2);
        }
        if (0 <= extrema3 && extrema3 <= 1) { // check if this extrema is inside the time range [0,1]
            minPosition = min(minPosition, extremum3);
            maxPosition = max(maxPosition, extremum3);
        }

        extremum[i              ] = minPosition;
        extremum[i + NUM_FACTORS] = maxPosition;
    }
}

void BezierCurve::returnJointPositionExtremumGradient(TYPE* extremumGradient, const TYPE* k) {
    for (int i = 0; i < NUM_FACTORS; i++) {
        // k[i] range is [-1,1] since it is defined for PZ, the following is the actual k
        TYPE k_actual = k_range[i] * k[i];

        // list all possible extremas
        TYPE extrema1 = 0;
        TYPE extrema2 = (2*qd0[i] + qdd0[i] + sqrt(64*pow(qd0[i],2) + 14*qd0[i]*qdd0[i] - 120*k_actual*qd0[i] + pow(qdd0[i],2))) / (5*(6*qd0[i] - 12*k_actual + qdd0[i]));
        TYPE extrema3 = (2*qd0[i] + qdd0[i] - sqrt(64*pow(qd0[i],2) + 14*qd0[i]*qdd0[i] - 120*k_actual*qd0[i] + pow(qdd0[i],2))) / (5*(6*qd0[i] - 12*k_actual + qdd0[i]));
        TYPE extrema4 = 1;

        // get extremums of all extremas
        TYPE extremum1 = q_des_func(q0[i], qd0[i], qdd0[i], k_actual, extrema1);
        TYPE extremum2 = q_des_func(q0[i], qd0[i], qdd0[i], k_actual, extrema2);
        TYPE extremum3 = q_des_func(q0[i], qd0[i], qdd0[i], k_actual, extrema3);
        TYPE extremum4 = q_des_func(q0[i], qd0[i], qdd0[i], k_actual, extrema4);

        // find the min and max values
        TYPE minPosition;
        int minId;
        TYPE maxPosition;
        int maxId;

        if (extremum1 < extremum4) {
            minPosition = extremum1;
            minId = 1;

            maxPosition = extremum4;
            maxId = 4;
        }
        else {
            minPosition = extremum4;
            minId = 4;

            maxPosition = extremum1;
            maxId = 1;
        }

        if (0 <= extrema2 && extrema2 <= 1) { // check if this extrema is inside the time range [0,1]
            if (extremum2 < minPosition) {
                minPosition = extremum2;
                minId = 2;
            }
            if (maxPosition < extremum2) {
                maxPosition = extremum2;
                maxId = 2;
            }
        }
        if (0 <= extrema3 && extrema3 <= 1) { // check if this extrema is inside the time range [0,1]
            if (extremum3 < minPosition) {
                minPosition = extremum3;
                minId = 3;
            }
            if (maxPosition < extremum3) {
                maxPosition = extremum3;
                maxId = 3;
            }
        }

        TYPE minPositionGradient;
        TYPE maxPositionGradient;

        switch (minId) {
            case 1: // t = 0
                minPositionGradient = 0.0;
                break;
            case 2: // t = extrema2
                minPositionGradient = q_des_extrema2_k_derivative(q0[i], qd0[i], qdd0[i], k_actual);
                break;
            case 3: // t = extrema3
                minPositionGradient = q_des_extrema3_k_derivative(q0[i], qd0[i], qdd0[i], k_actual);
                break;
            case 4: // t = 1
                minPositionGradient = 1.0;
                break;
            default:
                break;
        }

        switch (maxId) {
            case 1: // t = 0
                maxPositionGradient = 0.0;
                break;
            case 2: // t = extrema2
                maxPositionGradient = q_des_extrema2_k_derivative(q0[i], qd0[i], qdd0[i], k_actual);
                break;
            case 3: // t = extrema3
                maxPositionGradient = q_des_extrema3_k_derivative(q0[i], qd0[i], qdd0[i], k_actual);
                break;
            case 4: // t = 1
                maxPositionGradient = 1.0;
                break;
            default:
                break;
        }

        for (int j = 0; j < NUM_FACTORS; j++) {
            if (i == j) {
                extremumGradient[(i              ) * NUM_FACTORS + j] = minPositionGradient * k_range[i];
                extremumGradient[(i + NUM_FACTORS) * NUM_FACTORS + j] = maxPositionGradient * k_range[i];
            }
            else {
                extremumGradient[(i              ) * NUM_FACTORS + j] = 0.0;
                extremumGradient[(i + NUM_FACTORS) * NUM_FACTORS + j] = 0.0;
            }
        }
    }
}

void BezierCurve::returnJointVelocityExtremum(TYPE* extremum, const TYPE* k) {
    for (int i = 0; i < NUM_FACTORS; i++) {
        // k[i] range is [-1,1] since it is defined for PZ, the following is the actual k
        TYPE k_actual = k_range[i] * k[i];

        // list all possible extremas
        TYPE extrema1 = 0;
        TYPE extrema2 = (18*qd0[i] - 30*k_actual + 4*qdd0[i] + sqrt(6*(150*pow(k_actual,2) - 180*k_actual*qd0[i] - 20*k_actual*qdd0[i] + 54*pow(qd0[i],2) + 14*qd0[i]*qdd0[i] + pow(qdd0[i],2))))/(10*(6*qd0[i] - 12*k_actual + qdd0[i]));
        TYPE extrema3 = (18*qd0[i] - 30*k_actual + 4*qdd0[i] - sqrt(6*(150*pow(k_actual,2) - 180*k_actual*qd0[i] - 20*k_actual*qdd0[i] + 54*pow(qd0[i],2) + 14*qd0[i]*qdd0[i] + pow(qdd0[i],2))))/(10*(6*qd0[i] - 12*k_actual + qdd0[i]));
        TYPE extrema4 = 1;

        // get extremums of all extremas
        TYPE extremum1 = qd_des_func(q0[i], qd0[i], qdd0[i], k_actual, extrema1);
        TYPE extremum2 = qd_des_func(q0[i], qd0[i], qdd0[i], k_actual, extrema2);
        TYPE extremum3 = qd_des_func(q0[i], qd0[i], qdd0[i], k_actual, extrema3);
        TYPE extremum4 = qd_des_func(q0[i], qd0[i], qdd0[i], k_actual, extrema4);

        // find the min and max values
        TYPE minVelocity = min(extremum1, extremum4);
        TYPE maxVelocity = max(extremum1, extremum4);
        if (0 <= extrema2 && extrema2 <= 1) { // check if this extrema is inside the time range [0,1]
            minVelocity = min(minVelocity, extremum2); 
            maxVelocity = max(maxVelocity, extremum2);
        }
        if (0 <= extrema3 && extrema3 <= 1) { // check if this extrema is inside the time range [0,1]
            minVelocity = min(minVelocity, extremum3);
            maxVelocity = max(maxVelocity, extremum3);
        }

        extremum[i              ] = minVelocity;
        extremum[i + NUM_FACTORS] = maxVelocity;
    }
}

void BezierCurve::returnJointVelocityExtremumGradient(TYPE* extremumGradient, const TYPE* k) {
    for (int i = 0; i < NUM_FACTORS; i++) {
        // k[i] range is [-1,1] since it is defined for PZ, the following is the actual k
        TYPE k_actual = k_range[i] * k[i];

        // list all possible extremas
        TYPE extrema1 = 0;
        TYPE extrema2 = (18*qd0[i] - 30*k_actual + 4*qdd0[i] + sqrt(6*(150*pow(k_actual,2) - 180*k_actual*qd0[i] - 20*k_actual*qdd0[i] + 54*pow(qd0[i],2) + 14*qd0[i]*qdd0[i] + pow(qdd0[i],2))))/(10*(6*qd0[i] - 12*k_actual + qdd0[i]));
        TYPE extrema3 = (18*qd0[i] - 30*k_actual + 4*qdd0[i] - sqrt(6*(150*pow(k_actual,2) - 180*k_actual*qd0[i] - 20*k_actual*qdd0[i] + 54*pow(qd0[i],2) + 14*qd0[i]*qdd0[i] + pow(qdd0[i],2))))/(10*(6*qd0[i] - 12*k_actual + qdd0[i]));
        TYPE extrema4 = 1;

        // get extremums of all extremas
        TYPE extremum1 = qd_des_func(q0[i], qd0[i], qdd0[i], k_actual, extrema1);
        TYPE extremum2 = qd_des_func(q0[i], qd0[i], qdd0[i], k_actual, extrema2);
        TYPE extremum3 = qd_des_func(q0[i], qd0[i], qdd0[i], k_actual, extrema3);
        TYPE extremum4 = qd_des_func(q0[i], qd0[i], qdd0[i], k_actual, extrema4);

        // find the min and max values
        TYPE minVelocity;
        int minId;
        TYPE maxVelocity;
        int maxId;

        if (extremum1 < extremum4) {
            minVelocity = extremum1;
            minId = 1;

            maxVelocity = extremum4;
            maxId = 4;
        }
        else {
            minVelocity = extremum4;
            minId = 4;

            maxVelocity = extremum1;
            maxId = 1;
        }

        if (0 <= extrema2 && extrema2 <= 1) { // check if this extrema is inside the time range [0,1]
            if (extremum2 < minVelocity) {
                minVelocity = extremum2;
                minId = 2;
            }
            if (maxVelocity < extremum2) {
                maxVelocity = extremum2;
                maxId = 2;
            }
        }
        if (0 <= extrema3 && extrema3 <= 1) { // check if this extrema is inside the time range [0,1]
            if (extremum3 < minVelocity) {
                minVelocity = extremum3;
                minId = 3;
            }
            if (maxVelocity < extremum3) {
                maxVelocity = extremum3;
                maxId = 3;
            }
        }

        TYPE minVelocityGradient;
        TYPE maxVelocityGradient;

        switch (minId) {
            case 1: // t = 0
                minVelocityGradient = 0.0;
                break;
            case 2: // t = extrema2
                minVelocityGradient = qd_des_extrema2_k_derivative(q0[i], qd0[i], qdd0[i], k_actual);
                break;
            case 3: // t = extrema3
                minVelocityGradient = qd_des_extrema3_k_derivative(q0[i], qd0[i], qdd0[i], k_actual);
                break;
            case 4: // t = 1
                minVelocityGradient = 1.0;
                break;
            default:
                break;
        }

        switch (maxId) {
            case 1: // t = 0
                maxVelocityGradient = 0.0;
                break;
            case 2: // t = extrema2
                maxVelocityGradient = qd_des_extrema2_k_derivative(q0[i], qd0[i], qdd0[i], k_actual);
                break;
            case 3: // t = extrema3
                maxVelocityGradient = qd_des_extrema3_k_derivative(q0[i], qd0[i], qdd0[i], k_actual);
                break;
            case 4: // t = 1
                maxVelocityGradient = 1.0;
                break;
            default:
                break;
        }

        for (int j = 0; j < NUM_FACTORS; j++) {
            if (i == j) {
                extremumGradient[(i              ) * NUM_FACTORS + j] = minVelocityGradient * k_range[i];
                extremumGradient[(i + NUM_FACTORS) * NUM_FACTORS + j] = maxVelocityGradient * k_range[i];
            }
            else {
                extremumGradient[(i              ) * NUM_FACTORS + j] = 0.0;
                extremumGradient[(i + NUM_FACTORS) * NUM_FACTORS + j] = 0.0;
            }
        }
    }
}

TYPE q_des_func(TYPE q0, TYPE qd0, TYPE qdd0, TYPE k, TYPE t) {
    TYPE B0 = -pow(t - 1,5);
    TYPE B1 = 5*t*pow(t - 1,4);
    TYPE B2 = -10*pow(t,2)*pow(t - 1,3);
    TYPE B3 = 10*pow(t,3)*pow(t - 1,2);
    TYPE B4 = -5*pow(t,4)*(t - 1);
    TYPE B5 = pow(t,5);
    TYPE beta0 = q0;
    TYPE beta1 = q0 + qd0/5;
    TYPE beta2 = q0 + (2*qd0)/5 + qdd0/20;
    TYPE beta3 = q0 + k;
    TYPE beta4 = q0 + k;
    TYPE beta5 = q0 + k;
    return B0 * beta0 + B1 * beta1 + B2 * beta2 + B3 * beta3 + B4 * beta4 + B5 * beta5;
}

TYPE qd_des_func(TYPE q0, TYPE qd0, TYPE qdd0, TYPE k, TYPE t) {
    TYPE dB0 = pow(t-1.0,4.0)*-5.0;
    TYPE dB1 = t*pow(t-1.0,3.0)*2.0E+1+pow(t-1.0,4.0)*5.0;
    TYPE dB2 = t*pow(t-1.0,3.0)*-2.0E+1-(t*t)*pow(t-1.0,2.0)*3.0E+1;
    TYPE dB3 = pow(t,3.0)*(t*2.0-2.0)*1.0E+1+(t*t)*pow(t-1.0,2.0)*3.0E+1;
    TYPE dB4 = pow(t,3.0)*(t-1.0)*-2.0E+1-pow(t,4.0)*5.0;
    TYPE dB5 = pow(t,4.0)*5.0;
    TYPE beta0 = q0;
    TYPE beta1 = q0 + qd0/5;
    TYPE beta2 = q0 + (2*qd0)/5 + qdd0/20;
    TYPE beta3 = q0 + k;
    TYPE beta4 = q0 + k;
    TYPE beta5 = q0 + k;
    return dB0 * beta0 + dB1 * beta1 + dB2 * beta2 + dB3 * beta3 + dB4 * beta4 + dB5 * beta5;
}

TYPE qdd_des_func(TYPE q0, TYPE qd0, TYPE qdd0, TYPE k, TYPE t) {
    TYPE t2 = t*2.0;
    TYPE t3 = t*t;
    TYPE t4 = t*t*t;
    TYPE t5 = t-1.0;
    TYPE t6 = t2-2.0;
    TYPE t7 = t4*2.0E+1;
    TYPE t8 = t5*t5;
    TYPE t9 = t5*t5*t5;
    TYPE t10 = t9*2.0E+1;
    TYPE t11 = t*t8*6.0E+1;
    TYPE t12 = -t10;
    TYPE ddB0 = t12;
    TYPE ddB1 = t9*4.0E+1+t11;
    TYPE ddB2 = t12-t*t8*1.2E+2-t3*t6*3.0E+1;
    TYPE ddB3 = t7+t11+t3*t6*6.0E+1;
    TYPE ddB4 = t4*-4.0E+1-t3*t5*6.0E+1;
    TYPE ddB5 = t7;
    TYPE beta0 = q0;
    TYPE beta1 = q0 + qd0/5;
    TYPE beta2 = q0 + (2*qd0)/5 + qdd0/20;
    TYPE beta3 = q0 + k;
    TYPE beta4 = q0 + k;
    TYPE beta5 = q0 + k;
    return ddB0 * beta0 + ddB1 * beta1 + ddB2 * beta2 + ddB3 * beta3 + ddB4 * beta4 + ddB5 * beta5;
}

TYPE q_des_extrema2_k_derivative(TYPE q0, TYPE qd0, TYPE qdd0, TYPE k) {
    TYPE t2 = k+q0;
    TYPE t3 = qd0*2.0;
    TYPE t4 = qd0*6.0;
    TYPE t5 = qd0*qd0;
    TYPE t6 = qdd0*qdd0;
    TYPE t7 = k*1.2E+1;
    TYPE t8 = qd0*qdd0*1.4E+1;
    TYPE t10 = qd0/5.0;
    TYPE t11 = qd0*(2.0/5.0);
    TYPE t12 = k*qd0*1.2E+2;
    TYPE t13 = qdd0/2.0E+1;
    TYPE t9 = -t7;
    TYPE t14 = t5*6.4E+1;
    TYPE t15 = -t12;
    TYPE t16 = q0+t10;
    TYPE t18 = q0+t11+t13;
    TYPE t17 = qdd0+t4+t9;
    TYPE t24 = t6+t8+t14+t15;
    TYPE t19 = 1.0/t17;
    TYPE t25 = sqrt(t24);
    TYPE t20 = t19*t19;
    TYPE t21 = t19*t19*t19;
    TYPE t23 = t19*t19*t19*t19*t19;
    TYPE t26 = 1.0/t25;
    TYPE t27 = qdd0+t3+t25;
    TYPE t22 = t20*t20;
    TYPE t28 = t27*t27;
    TYPE t29 = t27*t27*t27;
    TYPE t31 = t27*t27*t27*t27*t27;
    TYPE t32 = qd0*t19*t26*1.2E+1;
    TYPE t34 = (t19*t27)/5.0;
    TYPE t35 = t20*t27*(1.2E+1/5.0);
    TYPE t30 = t28*t28;
    TYPE t33 = -t32;
    TYPE t36 = t34-1.0;
    TYPE t37 = t36*t36;
    TYPE t38 = t36*t36*t36;
    TYPE t40 = t33+t35;
    TYPE t39 = t37*t37;
    return (t23*t31)/3.125E+3+t2*(t20*t20*t20)*t31*(1.2E+1/6.25E+2)+t21*t29*t37*(2.0/2.5E+1)-(t22*t30*t36)/1.25E+2+q0*t39*(t32-t35)*5.0+t2*t22*t29*t37*(7.2E+1/2.5E+1)-t2*t23*t30*t36*(4.8E+1/1.25E+2)+t16*t20*t27*t39*1.2E+1-t18*t21*t28*t38*(4.8E+1/5.0)+(t2*t22*t30*(t32-t35))/1.25E+2-qd0*t2*t23*t26*t30*(1.2E+1/1.25E+2)-qd0*t16*t19*t26*t39*6.0E+1-t2*t21*t29*t36*(t32-t35)*(4.0/2.5E+1)-t16*t19*t27*t38*(t32-t35)*4.0+t18*t20*t28*t37*(t32-t35)*(6.0/5.0)-qd0*t2*t21*t26*t28*t37*(7.2E+1/5.0)+qd0*t2*t22*t26*t29*t36*(4.8E+1/2.5E+1)+qd0*t18*t20*t26*t27*t38*4.8E+1;
}

TYPE q_des_extrema3_k_derivative(TYPE q0, TYPE qd0, TYPE qdd0, TYPE k) {
    TYPE t2 = k+q0;
    TYPE t3 = qd0*2.0;
    TYPE t4 = qd0*6.0;
    TYPE t5 = qd0*qd0;
    TYPE t6 = qdd0*qdd0;
    TYPE t7 = k*1.2E+1;
    TYPE t8 = qd0*qdd0*1.4E+1;
    TYPE t10 = qd0/5.0;
    TYPE t11 = qd0*(2.0/5.0);
    TYPE t12 = k*qd0*1.2E+2;
    TYPE t13 = qdd0/2.0E+1;
    TYPE t9 = -t7;
    TYPE t14 = t5*6.4E+1;
    TYPE t15 = -t12;
    TYPE t16 = q0+t10;
    TYPE t18 = q0+t11+t13;
    TYPE t17 = qdd0+t4+t9;
    TYPE t24 = t6+t8+t14+t15;
    TYPE t19 = 1.0/t17;
    TYPE t25 = sqrt(t24);
    TYPE t20 = t19*t19;
    TYPE t21 = t19*t19*t19;
    TYPE t23 = t19*t19*t19*t19*t19;
    TYPE t26 = 1.0/t25;
    TYPE t27 = -t25;
    TYPE t22 = t20*t20;
    TYPE t28 = qdd0+t3+t27;
    TYPE t33 = qd0*t19*t26*1.2E+1;
    TYPE t29 = t28*t28;
    TYPE t30 = t28*t28*t28;
    TYPE t32 = t28*t28*t28*t28*t28;
    TYPE t34 = (t19*t28)/5.0;
    TYPE t35 = t20*t28*(1.2E+1/5.0);
    TYPE t31 = t29*t29;
    TYPE t36 = t34-1.0;
    TYPE t40 = t33+t35;
    TYPE t37 = t36*t36;
    TYPE t38 = t36*t36*t36;
    TYPE t39 = t37*t37;
    return (t23*t32)/3.125E+3+t2*(t20*t20*t20)*t32*(1.2E+1/6.25E+2)-q0*t39*t40*5.0+t21*t30*t37*(2.0/2.5E+1)-(t22*t31*t36)/1.25E+2+t2*t22*t30*t37*(7.2E+1/2.5E+1)-t2*t23*t31*t36*(4.8E+1/1.25E+2)-(t2*t22*t31*t40)/1.25E+2+t16*t20*t28*t39*1.2E+1-t18*t21*t29*t38*(4.8E+1/5.0)+qd0*t2*t23*t26*t31*(1.2E+1/1.25E+2)+qd0*t16*t19*t26*t39*6.0E+1+t2*t21*t30*t36*t40*(4.0/2.5E+1)+t16*t19*t28*t38*t40*4.0-t18*t20*t29*t37*t40*(6.0/5.0)+qd0*t2*t21*t26*t29*t37*(7.2E+1/5.0)-qd0*t2*t22*t26*t30*t36*(4.8E+1/2.5E+1)-qd0*t18*t20*t26*t28*t38*4.8E+1;
}

TYPE qd_des_extrema2_k_derivative(TYPE q0, TYPE qd0, TYPE qdd0, TYPE k) {
    TYPE t2 = k+q0;
    TYPE t3 = k*k;
    TYPE t4 = qd0*6.0;
    TYPE t5 = qdd0*4.0;
    TYPE t6 = qd0*qd0;
    TYPE t7 = qdd0*qdd0;
    TYPE t8 = k*1.2E+1;
    TYPE t9 = k*3.0E+1;
    TYPE t10 = qd0*1.8E+1;
    TYPE t11 = qdd0*2.0E+1;
    TYPE t13 = qd0*qdd0*1.4E+1;
    TYPE t14 = sqrt(6.0);
    TYPE t17 = k*3.0E+2;
    TYPE t18 = qd0*1.8E+2;
    TYPE t21 = k*qdd0*-2.0E+1;
    TYPE t24 = k*qd0*-1.8E+2;
    TYPE t12 = k*t11;
    TYPE t15 = -t8;
    TYPE t16 = -t9;
    TYPE t19 = t6*5.4E+1;
    TYPE t20 = k*t18;
    TYPE t22 = -t17;
    TYPE t23 = t3*1.5E+2;
    TYPE t25 = qdd0+t4+t15;
    TYPE t31 = t11+t18+t22;
    TYPE t32 = t7+t13+t19+t21+t23+t24;
    TYPE t26 = 1.0/t25;
    TYPE t33 = sqrt(t32);
    TYPE t27 = t26*t26;
    TYPE t28 = t26*t26*t26;
    TYPE t30 = t26*t26*t26*t26*t26;
    TYPE t34 = 1.0/t33;
    TYPE t35 = t14*t33;
    TYPE t29 = t27*t27;
    TYPE t36 = t5+t10+t16+t35;
    TYPE t40 = (t14*t31*t34)/2.0;
    TYPE t37 = t36*t36;
    TYPE t38 = t36*t36*t36;
    TYPE t41 = t40+3.0E+1;
    TYPE t42 = (t26*t36)/5.0;
    TYPE t43 = t27*t36*(6.0/5.0);
    TYPE t44 = (t26*t36)/1.0E+1;
    TYPE t39 = t37*t37;
    TYPE t45 = t42-2.0;
    TYPE t46 = t44-1.0;
    TYPE t49 = (t26*t41)/1.0E+1;
    TYPE t47 = t46*t46;
    TYPE t48 = t46*t46*t46;
    TYPE t50 = -t49;
    TYPE t51 = t27*t36*t48*2.4E+1;
    TYPE t52 = t28*t37*t47*(3.6E+1/5.0);
    TYPE t53 = t43+t50;
    TYPE t54 = t26*t41*t48*2.0;
    TYPE t56 = t27*t36*t41*t47*(3.0/5.0);
    TYPE t55 = -t54;
    TYPE t57 = -t56;
    TYPE t58 = t26*t36*t47*t53*6.0;
    TYPE t59 = t27*t37*t46*t53*(3.0/5.0);
    return (q0+qd0/5.0)*(t51+t55+t58+t48*t53*2.0E+1)+t2*(t52+t57+t59+(t28*t38*(t27*t36*(1.2E+1/5.0)-(t26*t41)/5.0))/1.0E+2+t29*t38*t45*(9.0/2.5E+1)-t28*t37*t41*t45*(3.0/1.0E+2))-t2*(t30*t39*(3.0/1.25E+2)-(t29*t38*t41)/5.0E+2+t29*t38*t46*(1.8E+1/2.5E+1)+(t28*t38*t53)/5.0E+1-t28*t37*t41*t46*(3.0/5.0E+1))-(q0+qd0*(2.0/5.0)+qdd0/2.0E+1)*(t51+t52+t55+t57+t58+t59)-q0*t48*t53*2.0E+1+t2*t30*t39*(3.0/1.25E+2)+t27*t37*t47*(3.0/1.0E+1)+(t28*t38*t45)/1.0E+2-(t28*t38*t46)/5.0E+1-(t2*t29*t38*t41)/5.0E+2;
}

TYPE qd_des_extrema3_k_derivative(TYPE q0, TYPE qd0, TYPE qdd0, TYPE k) {
    TYPE t2 = k+q0;
    TYPE t3 = k*k;
    TYPE t4 = qd0*6.0;
    TYPE t5 = qdd0*4.0;
    TYPE t6 = qd0*qd0;
    TYPE t7 = qdd0*qdd0;
    TYPE t8 = k*1.2E+1;
    TYPE t9 = k*3.0E+1;
    TYPE t10 = qd0*1.8E+1;
    TYPE t12 = qdd0*2.0E+1;
    TYPE t14 = qd0*qdd0*1.4E+1;
    TYPE t15 = sqrt(6.0);
    TYPE t17 = k*3.0E+2;
    TYPE t19 = qd0*1.8E+2;
    TYPE t22 = k*qdd0*-2.0E+1;
    TYPE t25 = k*qd0*-1.8E+2;
    TYPE t11 = -t5;
    TYPE t13 = k*t12;
    TYPE t16 = -t8;
    TYPE t18 = -t10;
    TYPE t20 = t6*5.4E+1;
    TYPE t21 = k*t19;
    TYPE t23 = -t17;
    TYPE t24 = t3*1.5E+2;
    TYPE t26 = qdd0+t4+t16;
    TYPE t32 = t12+t19+t23;
    TYPE t33 = t7+t14+t20+t22+t24+t25;
    TYPE t27 = 1.0/t26;
    TYPE t34 = sqrt(t33);
    TYPE t28 = t27*t27;
    TYPE t29 = t27*t27*t27;
    TYPE t31 = t27*t27*t27*t27*t27;
    TYPE t35 = 1.0/t34;
    TYPE t36 = t15*t34;
    TYPE t30 = t28*t28;
    TYPE t37 = t9+t11+t18+t36;
    TYPE t38 = pow(t5-t9+t10-t36,2.0);
    TYPE t39 = -pow(t5-t9+t10-t36,3.0);
    TYPE t41 = (t15*t32*t35)/2.0;
    TYPE t43 = t27*(t5-t9+t10-t36)*(-1.0/5.0);
    TYPE t44 = t28*(t5-t9+t10-t36)*(-6.0/5.0);
    TYPE t45 = t27*(t5-t9+t10-t36)*(-1.0/1.0E+1);
    TYPE t48 = pow((t27*(t5-t9+t10-t36))/1.0E+1-1.0,2.0);
    TYPE t49 = -pow((t27*(t5-t9+t10-t36))/1.0E+1-1.0,3.0);
    TYPE t52 = t28*pow((t27*(t5-t9+t10-t36))/1.0E+1-1.0,3.0)*(t5-t9+t10-t36)*2.4E+1;
    TYPE t40 = t38*t38;
    TYPE t42 = t41-3.0E+1;
    TYPE t46 = t43+2.0;
    TYPE t47 = t45+1.0;
    TYPE t53 = t29*t38*t48*(3.6E+1/5.0);
    TYPE t50 = (t27*t42)/1.0E+1;
    TYPE t55 = t27*t42*pow((t27*(t5-t9+t10-t36))/1.0E+1-1.0,3.0)*-2.0;
    TYPE t56 = t27*t42*pow((t27*(t5-t9+t10-t36))/1.0E+1-1.0,3.0)*2.0;
    TYPE t57 = t28*t42*t48*(t5-t9+t10-t36)*(-3.0/5.0);
    TYPE t58 = t28*t42*t48*(t5-t9+t10-t36)*(3.0/5.0);
    TYPE t51 = -t50;
    TYPE t59 = t27*t48*(t50+t28*(t5-t9+t10-t36)*(6.0/5.0))*(t5-t9+t10-t36)*6.0;
    TYPE t60 = t28*t38*((t27*(t5-t9+t10-t36))/1.0E+1-1.0)*(t50+t28*(t5-t9+t10-t36)*(6.0/5.0))*(3.0/5.0);
    TYPE t54 = t44+t51;
    return (q0+qd0/5.0)*(t52+t56+t59+pow((t27*(t5-t9+t10-t36))/1.0E+1-1.0,3.0)*(t50+t28*(t5-t9+t10-t36)*(6.0/5.0))*2.0E+1)+t2*(t53+t58+t60+t30*((t27*(t5-t9+t10-t36))/5.0-2.0)*pow(t5-t9+t10-t36,3.0)*(9.0/2.5E+1)+(t29*((t27*t42)/5.0+t28*(t5-t9+t10-t36)*(1.2E+1/5.0))*pow(t5-t9+t10-t36,3.0))/1.0E+2+t29*t38*t42*((t27*(t5-t9+t10-t36))/5.0-2.0)*(3.0/1.0E+2))-t2*(t31*t40*(3.0/1.25E+2)+t30*((t27*(t5-t9+t10-t36))/1.0E+1-1.0)*pow(t5-t9+t10-t36,3.0)*(1.8E+1/2.5E+1)+(t30*t42*pow(t5-t9+t10-t36,3.0))/5.0E+2+(t29*(t50+t28*(t5-t9+t10-t36)*(6.0/5.0))*pow(t5-t9+t10-t36,3.0))/5.0E+1+t29*t38*t42*((t27*(t5-t9+t10-t36))/1.0E+1-1.0)*(3.0/5.0E+1))-(q0+qd0*(2.0/5.0)+qdd0/2.0E+1)*(t52+t53+t56+t58+t59+t60)+(t29*((t27*(t5-t9+t10-t36))/5.0-2.0)*pow(t5-t9+t10-t36,3.0))/1.0E+2-(t29*((t27*(t5-t9+t10-t36))/1.0E+1-1.0)*pow(t5-t9+t10-t36,3.0))/5.0E+1+t2*t31*t40*(3.0/1.25E+2)+t28*t38*t48*(3.0/1.0E+1)-q0*pow((t27*(t5-t9+t10-t36))/1.0E+1-1.0,3.0)*(t50+t28*(t5-t9+t10-t36)*(6.0/5.0))*2.0E+1+(t2*t30*t42*pow(t5-t9+t10-t36,3.0))/5.0E+2;
}

TYPE q_des_k_indep(TYPE q0, TYPE qd0, TYPE qdd0, TYPE t) {
    return q0 + qd0*t - 6*qd0*pow(t,3) + 8*qd0*pow(t,4) - 3*qd0*pow(t,5) + (qdd0*pow(t,2))*0.5 - (3*qdd0*pow(t,3))*0.5 + (3*qdd0*pow(t,4))*0.5 - (qdd0*pow(t,5))*0.5;
}

TYPE qd_des_k_indep(TYPE q0, TYPE qd0, TYPE qdd0, TYPE t) {
    return (pow(t - 1,2)*(2*qd0 + 4*qd0*t + 2*qdd0*t - 30*qd0*pow(t,2) - 5*qdd0*pow(t,2)))*0.5;
}

TYPE qdd_des_k_indep(TYPE q0, TYPE qd0, TYPE qdd0, TYPE t) {
    return -(t - 1.0)*(qdd0 - (36*qd0 + 8*qdd0)*t + (60*qd0 + 10*qdd0)*pow(t, 2));
}

#endif