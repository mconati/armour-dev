#ifndef DYNAMICS_CPP
#define DYNAMICS_CPP

#include "Dynamics.h"

KinematicsDynamics::KinematicsDynamics(BezierCurve* traj_input) {
    traj = traj_input;

    // pre-allocate memory
    links = PZsparseArray(NUM_FACTORS * 3, NUM_TIME_STEPS);
    mass_nominal_arr = PZsparseArray(NUM_JOINTS, 1);
    mass_uncertain_arr = PZsparseArray(NUM_JOINTS, 1);
    I_nominal_arr = PZsparseArray(NUM_JOINTS, 1);
    I_uncertain_arr = PZsparseArray(NUM_JOINTS, 1);
    u_nom = PZsparseArray(NUM_FACTORS, NUM_TIME_STEPS);
    u_nom_int = PZsparseArray(NUM_FACTORS, NUM_TIME_STEPS);

    // initialize robot properties
    for (int i = 0; i < NUM_JOINTS; i++) {
        trans_matrix(i, 0) = Eigen::MatrixXd::Zero(3, 1);
        trans_matrix(i, 0)(0) = trans[3 * i];
        trans_matrix(i, 0)(1) = trans[3 * i + 1];
        trans_matrix(i, 0)(2) = trans[3 * i + 2];

        // com_matrix(i, 0) = Eigen::MatrixXd::Zero(3, 1);
        // com_matrix(i, 0)(0) = com[3 * i];
        // com_matrix(i, 0)(1) = com[3 * i + 1];
        // com_matrix(i, 0)(2) = com[3 * i + 2];

        Eigen::MatrixXd mass_matrix(1, 1);
        mass_matrix(0) = mass[i];
        mass_nominal_arr(i) = PZsparse(mass_matrix);
        mass_uncertain_arr(i) = PZsparse(mass_matrix, mass_uncertainty);

        Eigen::Matrix3d inertia_matrix;
        for (int j = 0; j < 9; j++) {
            inertia_matrix(j) = inertia[i * 9 + j]; // This may not be right...
        }
        I_nominal_arr(i) = PZsparse(inertia_matrix);
        I_uncertain_arr(i) = PZsparse(inertia_matrix, inertia_uncertainty);
    }

    trans_matrix(NUM_JOINTS, 0) = Eigen::MatrixXd::Zero(3, 1);
    trans_matrix(NUM_JOINTS, 0)(0) = trans[3 * NUM_JOINTS];
    trans_matrix(NUM_JOINTS, 0)(1) = trans[3 * NUM_JOINTS + 1];
    trans_matrix(NUM_JOINTS, 0)(2) = trans[3 * NUM_JOINTS + 2];

    // define original link PZs
    links = PZsparseArray(NUM_JOINTS, NUM_TIME_STEPS);

    for (int i = 0; i < NUM_JOINTS; i++) {
        PZsparseArray link(3, 1);

        for (int j = 0; j < 3; j++) {
            uint64_t degree[1][NUM_FACTORS * 6] = {0};
            degree[0][NUM_FACTORS * (j + 1)] = 1; // use qde, qdae, qdde for x, y, z generator
            double temp = link_zonotope_generators[i][j];
            link(j, 0) = PZsparse(link_zonotope_center[i][j], &temp, degree, 1);
        }

        links(i, 0) = stack(link);

        for (int j = 1; j < NUM_TIME_STEPS; j++) {
            links(i, j) = links(i, 0);
        }
    }
}

void KinematicsDynamics::fk(uint s_ind) {
    PZsparse FK_R = PZsparse(0, 0, 0); // identity matrix
    PZsparse FK_T(3, 1);
    
    for (int i = 0; i < NUM_JOINTS; i++) {
        PZsparse P(trans_matrix(i, 0));
        
        FK_T = FK_T + FK_R * P;
        FK_R = FK_R * traj->R(i, s_ind);
        
        links(i, s_ind) = FK_R * links(i, s_ind) + FK_T;
    }
}

void KinematicsDynamics::rnea(uint s_ind,
                              PZsparseArray& mass_arr,
                              PZsparseArray& I_arr,
                              PZsparseArray& u,
                              bool setGravity) {
    PZsparse& cq1 = traj->cos_q_des(0, s_ind);
    PZsparse& cq2 = traj->cos_q_des(1, s_ind);
    PZsparse& cq3 = traj->cos_q_des(2, s_ind);
    PZsparse& cq4 = traj->cos_q_des(3, s_ind);
    PZsparse& cq5 = traj->cos_q_des(4, s_ind);
    PZsparse& cq6 = traj->cos_q_des(5, s_ind);
    PZsparse& cq7 = traj->cos_q_des(6, s_ind);

    PZsparse& sq1 = traj->sin_q_des(0, s_ind);
    PZsparse& sq2 = traj->sin_q_des(1, s_ind);
    PZsparse& sq3 = traj->sin_q_des(2, s_ind);
    PZsparse& sq4 = traj->sin_q_des(3, s_ind);
    PZsparse& sq5 = traj->sin_q_des(4, s_ind);
    PZsparse& sq6 = traj->sin_q_des(5, s_ind);
    PZsparse& sq7 = traj->sin_q_des(6, s_ind);

    PZsparse& qd1 = traj->qd_des(0, s_ind);
    PZsparse& qd2 = traj->qd_des(1, s_ind);
    PZsparse& qd3 = traj->qd_des(2, s_ind);
    PZsparse& qd4 = traj->qd_des(3, s_ind);
    PZsparse& qd5 = traj->qd_des(4, s_ind);
    PZsparse& qd6 = traj->qd_des(5, s_ind);
    PZsparse& qd7 = traj->qd_des(6, s_ind);

    PZsparse& qda1 = traj->qda_des(0, s_ind);
    PZsparse& qda2 = traj->qda_des(1, s_ind);
    PZsparse& qda3 = traj->qda_des(2, s_ind);
    PZsparse& qda4 = traj->qda_des(3, s_ind);
    PZsparse& qda5 = traj->qda_des(4, s_ind);
    PZsparse& qda6 = traj->qda_des(5, s_ind);
    PZsparse& qda7 = traj->qda_des(6, s_ind);

    PZsparse& qdd1 = traj->qdda_des(0, s_ind);
    PZsparse& qdd2 = traj->qdda_des(1, s_ind);
    PZsparse& qdd3 = traj->qdda_des(2, s_ind);
    PZsparse& qdd4 = traj->qdda_des(3, s_ind);
    PZsparse& qdd5 = traj->qdda_des(4, s_ind);
    PZsparse& qdd6 = traj->qdda_des(5, s_ind);
    PZsparse& qdd7 = traj->qdda_des(6, s_ind);

    PZsparse w1;
    PZsparse w2;
    PZsparse w3;
    PZsparse w_aux1;
    PZsparse w_aux2;
    PZsparse w_aux3;
    PZsparse wdot1;
    PZsparse wdot2;
    PZsparse wdot3;
    PZsparse linear_acc1;
    PZsparse linear_acc2;
    PZsparse linear_acc3;

    PZsparse w_new1;
    PZsparse w_new2;
    PZsparse w_new3;
    PZsparse w_aux_new1;
    PZsparse w_aux_new2;
    PZsparse w_aux_new3;
    PZsparse wdot_new1;
    PZsparse wdot_new2;
    PZsparse wdot_new3;
    PZsparse linear_acc_new1;
    PZsparse linear_acc_new2;
    PZsparse linear_acc_new3;

    PZsparse t1;
    PZsparse t2;
    PZsparse t3;
    PZsparse t4;
    PZsparse t5;
    PZsparse t6;
    PZsparse t7;
    PZsparse t8;
    PZsparse t9;
    PZsparse t10;
    PZsparse t11;
    PZsparse t12;
    PZsparse t13;
    PZsparse t14;
    PZsparse t15;
    PZsparse t16;
    PZsparse t17;

    // joint 1
    w_new3 = qd1;

    w_aux_new3 = qda1;

    wdot_new3 = qdd1;

    if (setGravity) {
        linear_acc_new3 = PZsparse(-gravity);
    }

    w3 = w_new3;
    w_aux3 = w_aux_new3;
    wdot3 = wdot_new3;
    linear_acc3 = linear_acc_new3;

    t2 = com[0][0]*w_aux2;
    t3 = com[0][0]*w_aux3;
    t4 = com[0][1]*w_aux1;
    t5 = com[0][1]*w_aux3;
    t6 = com[0][2]*w_aux1;
    t7 = com[0][2]*w_aux2;
    t8 = -t4;
    t9 = -t6;
    t10 = -t7;
    t11 = t2+t8;
    t12 = t3+t9;
    t13 = t5+t10;
    PZsparse F1_1 = -mass_arr(0, 0)*(-linear_acc1+com[0][1]*wdot3-com[0][2]*wdot2+t11*w2+t12*w3);
    PZsparse F1_2 = mass_arr(0, 0)*(linear_acc2+com[0][0]*wdot3-com[0][2]*wdot1+t11*w1-t13*w3);
    PZsparse F1_3 = mass_arr(0, 0)*(linear_acc3-com[0][0]*wdot2+com[0][1]*wdot1+t12*w1+t13*w2);

    t2 = I_arr(0,0)*w1;
    t3 = I_arr(0,1)*w2;
    t4 = I_arr(0,2)*w3;
    t5 = I_arr(0,3)*w1;
    t6 = I_arr(0,4)*w2;
    t7 = I_arr(0,5)*w3;
    t8 = I_arr(0,6)*w1;
    t9 = I_arr(0,7)*w2;
    t10 = I_arr(0,8)*w3;
    t11 = t2+t3+t4;
    t12 = t5+t6+t7;
    t13 = t8+t9+t10;
    PZsparse N1_1 = I_arr(0,0)*wdot1+I_arr(0,1)*wdot2+I_arr(0,2)*wdot3-t12*w_aux3+t13*w_aux2;
    PZsparse N1_2 = I_arr(0,3)*wdot1+I_arr(0,4)*wdot2+I_arr(0,5)*wdot3+t11*w_aux3-t13*w_aux1;
    PZsparse N1_3 = I_arr(0,6)*wdot1+I_arr(0,7)*wdot2+I_arr(0,8)*wdot3-t11*w_aux2+t12*w_aux1;

    // joint 2
    w_new1 = cq2*w1+sq2*w3;
    w_new2 = cq2*w3-sq2*w1;
    w_new3 = qd2-w2;

    w_aux_new1 = cq2*w_aux1+sq2*w_aux3;
    w_aux_new2 = cq2*w_aux3-sq2*w_aux1;
    w_aux_new3 = qda2-w_aux2;

    wdot_new1 = cq2*wdot1+sq2*wdot3+qd2*(cq2*w_aux3-sq2*w_aux1);
    wdot_new2 = cq2*wdot3-sq2*wdot1-qd2*(cq2*w_aux1+sq2*w_aux3);
    wdot_new3 = qdd2-wdot2;

    t2 = -linear_acc1;
    t3 = w_aux3*5.375E-3;
    t4 = wdot1*5.375E-3;
    t5 = wdot3*5.375E-3;
    t6 = w_aux1*w2*5.375E-3;
    t8 = w_aux2*1.2838E-1;
    t9 = wdot2*1.2838E-1;
    t10 = w_aux1*w1*1.2838E-1;
    t11 = w_aux1*w3*1.2838E-1;
    t7 = -t6;
    t12 = t3+t8;
    t13 = t12*w2;
    t15 = t2+t5+t7+t9+t11;
    t14 = linear_acc3+t4+t10+t13;
    linear_acc_new1 = -cq2*t15+sq2*t14;
    linear_acc_new2 = cq2*t14+sq2*t15;
    linear_acc_new3 = -linear_acc2-wdot1*1.2838E-1+t12*w3+w_aux1*w1*5.375E-3;

    w1 = w_new1;
    w2 = w_new2;
    w3 = w_new3;
    w_aux1 = w_aux_new1;
    w_aux2 = w_aux_new2;
    w_aux3 = w_aux_new3;
    wdot1 = wdot_new1;
    wdot2 = wdot_new2;
    wdot3 = wdot_new3;
    linear_acc1 = linear_acc_new1;
    linear_acc2 = linear_acc_new2;
    linear_acc3 = linear_acc_new3;

    t2 = com[1][0]*w_aux2;
    t3 = com[1][0]*w_aux3;
    t4 = com[1][1]*w_aux1;
    t5 = com[1][1]*w_aux3;
    t6 = com[1][2]*w_aux1;
    t7 = com[1][2]*w_aux2;
    t8 = -t4;
    t9 = -t6;
    t10 = -t7;
    t11 = t2+t8;
    t12 = t3+t9;
    t13 = t5+t10;
    PZsparse F2_1 = -mass_arr(1, 0)*(-linear_acc1+com[1][1]*wdot3-com[1][2]*wdot2+t11*w2+t12*w3);
    PZsparse F2_2 = mass_arr(1, 0)*(linear_acc2+com[1][0]*wdot3-com[1][2]*wdot1+t11*w1-t13*w3);
    PZsparse F2_3 = mass_arr(1, 0)*(linear_acc3-com[1][0]*wdot2+com[1][1]*wdot1+t12*w1+t13*w2);

    t2 = I_arr(1,0)*w1;
    t3 = I_arr(1,1)*w2;
    t4 = I_arr(1,2)*w3;
    t5 = I_arr(1,3)*w1;
    t6 = I_arr(1,4)*w2;
    t7 = I_arr(1,5)*w3;
    t8 = I_arr(1,6)*w1;
    t9 = I_arr(1,7)*w2;
    t10 = I_arr(1,8)*w3;
    t11 = t2+t3+t4;
    t12 = t5+t6+t7;
    t13 = t8+t9+t10;
    PZsparse N2_1 = I_arr(1,0)*wdot1+I_arr(1,1)*wdot2+I_arr(1,2)*wdot3-t12*w_aux3+t13*w_aux2;
    PZsparse N2_2 = I_arr(1,3)*wdot1+I_arr(1,4)*wdot2+I_arr(1,5)*wdot3+t11*w_aux3-t13*w_aux1;
    PZsparse N2_3 = I_arr(1,6)*wdot1+I_arr(1,7)*wdot2+I_arr(1,8)*wdot3-t11*w_aux2+t12*w_aux1;

    // joint 3
    w_new1 = cq3*w1-sq3*w3;
    w_new2 = -cq3*w3-sq3*w1;
    w_new3 = qd3+w2;

    w_aux_new1 = cq3*w_aux1-sq3*w_aux3;
    w_aux_new2 = -cq3*w_aux3-sq3*w_aux1;
    w_aux_new3 = qda3+w_aux2;

    wdot_new1 = cq3*wdot1-sq3*wdot3-qd3*(cq3*w_aux3+sq3*w_aux1);
    wdot_new2 = -cq3*wdot3-sq3*wdot1-qd3*(cq3*w_aux1-sq3*w_aux3);
    wdot_new3 = qdd3+wdot2;

    t2 = -linear_acc1;
    t3 = w_aux2*6.375E-3;
    t4 = wdot2*6.375E-3;
    t5 = w_aux1*w1*6.375E-3;
    t6 = w_aux1*w3*6.375E-3;
    t7 = w_aux3*2.1038E-1;
    t8 = wdot1*2.1038E-1;
    t9 = wdot3*2.1038E-1;
    t10 = w_aux1*w2*2.1038E-1;
    t11 = -t7;
    t12 = -t8;
    t13 = -t9;
    t14 = t3+t11;
    t16 = t2+t4+t6+t10+t13;
    t15 = t14*w2;
    t17 = linear_acc3+t5+t12+t15;
    linear_acc_new1 = -cq3*t16-sq3*t17;
    linear_acc_new2 = -cq3*t17+sq3*t16;
    linear_acc_new3 = linear_acc2+wdot1*6.375E-3-t14*w3+w_aux1*w1*2.1038E-1;

    w1 = w_new1;
    w2 = w_new2;
    w3 = w_new3;
    w_aux1 = w_aux_new1;
    w_aux2 = w_aux_new2;
    w_aux3 = w_aux_new3;
    wdot1 = wdot_new1;
    wdot2 = wdot_new2;
    wdot3 = wdot_new3;
    linear_acc1 = linear_acc_new1;
    linear_acc2 = linear_acc_new2;
    linear_acc3 = linear_acc_new3;

    t2 = com[2][0]*w_aux2;
    t3 = com[2][0]*w_aux3;
    t4 = com[2][1]*w_aux1;
    t5 = com[2][1]*w_aux3;
    t6 = com[2][2]*w_aux1;
    t7 = com[2][2]*w_aux2;
    t8 = -t4;
    t9 = -t6;
    t10 = -t7;
    t11 = t2+t8;
    t12 = t3+t9;
    t13 = t5+t10;
    PZsparse F3_1 = -mass_arr(2, 0)*(-linear_acc1+com[2][1]*wdot3-com[2][2]*wdot2+t11*w2+t12*w3);
    PZsparse F3_2 = mass_arr(2, 0)*(linear_acc2+com[2][0]*wdot3-com[2][2]*wdot1+t11*w1-t13*w3);
    PZsparse F3_3 = mass_arr(2, 0)*(linear_acc3-com[2][0]*wdot2+com[2][1]*wdot1+t12*w1+t13*w2);

    t2 = I_arr(2,0)*w1;
    t3 = I_arr(2,1)*w2;
    t4 = I_arr(2,2)*w3;
    t5 = I_arr(2,3)*w1;
    t6 = I_arr(2,4)*w2;
    t7 = I_arr(2,5)*w3;
    t8 = I_arr(2,6)*w1;
    t9 = I_arr(2,7)*w2;
    t10 = I_arr(2,8)*w3;
    t11 = t2+t3+t4;
    t12 = t5+t6+t7;
    t13 = t8+t9+t10;
    PZsparse N3_1 = I_arr(2,0)*wdot1+I_arr(2,1)*wdot2+I_arr(2,2)*wdot3-t12*w_aux3+t13*w_aux2;
    PZsparse N3_2 = I_arr(2,3)*wdot1+I_arr(2,4)*wdot2+I_arr(2,5)*wdot3+t11*w_aux3-t13*w_aux1;
    PZsparse N3_3 = I_arr(2,6)*wdot1+I_arr(2,7)*wdot2+I_arr(2,8)*wdot3-t11*w_aux2+t12*w_aux1;

    // joint 4
    w_new1 = cq4*w1+sq4*w3;
    w_new2 = cq4*w3-sq4*w1;
    w_new3 = qd4-w2;

    w_aux_new1 = cq4*w_aux1+sq4*w_aux3;
    w_aux_new2 = cq4*w_aux3-sq4*w_aux1;
    w_aux_new3 = qda4-w_aux2;

    wdot_new1 = cq4*wdot1+sq4*wdot3+qd4*(cq4*w_aux3-sq4*w_aux1);
    wdot_new2 = cq4*wdot3-sq4*wdot1-qd4*(cq4*w_aux1+sq4*w_aux3);
    wdot_new3 = qdd4-wdot2;

    t2 = -linear_acc1;
    t3 = w_aux3*6.375E-3;
    t4 = wdot1*6.375E-3;
    t5 = wdot3*6.375E-3;
    t6 = w_aux1*w2*6.375E-3;
    t8 = w_aux2*2.1038E-1;
    t9 = wdot2*2.1038E-1;
    t10 = w_aux1*w1*2.1038E-1;
    t11 = w_aux1*w3*2.1038E-1;
    t7 = -t6;
    t12 = t3+t8;
    t13 = t12*w2;
    t15 = t2+t5+t7+t9+t11;
    t14 = linear_acc3+t4+t10+t13;
    linear_acc_new1 = -cq4*t15+sq4*t14;
    linear_acc_new2 = cq4*t14+sq4*t15;
    linear_acc_new3 = -linear_acc2-wdot1*2.1038E-1+t12*w3+w_aux1*w1*6.375E-3;

    w1 = w_new1;
    w2 = w_new2;
    w3 = w_new3;
    w_aux1 = w_aux_new1;
    w_aux2 = w_aux_new2;
    w_aux3 = w_aux_new3;
    wdot1 = wdot_new1;
    wdot2 = wdot_new2;
    wdot3 = wdot_new3;
    linear_acc1 = linear_acc_new1;
    linear_acc2 = linear_acc_new2;
    linear_acc3 = linear_acc_new3;

    t2 = com[3][0]*w_aux2;
    t3 = com[3][0]*w_aux3;
    t4 = com[3][1]*w_aux1;
    t5 = com[3][1]*w_aux3;
    t6 = com[3][2]*w_aux1;
    t7 = com[3][2]*w_aux2;
    t8 = -t4;
    t9 = -t6;
    t10 = -t7;
    t11 = t2+t8;
    t12 = t3+t9;
    t13 = t5+t10;
    PZsparse F4_1 = -mass_arr(3, 0)*(-linear_acc1+com[3][1]*wdot3-com[3][2]*wdot2+t11*w2+t12*w3);
    PZsparse F4_2 = mass_arr(3, 0)*(linear_acc2+com[3][0]*wdot3-com[3][2]*wdot1+t11*w1-t13*w3);
    PZsparse F4_3 = mass_arr(3, 0)*(linear_acc3-com[3][0]*wdot2+com[3][1]*wdot1+t12*w1+t13*w2);

    t2 = I_arr(3,0)*w1;
    t3 = I_arr(3,1)*w2;
    t4 = I_arr(3,2)*w3;
    t5 = I_arr(3,3)*w1;
    t6 = I_arr(3,4)*w2;
    t7 = I_arr(3,5)*w3;
    t8 = I_arr(3,6)*w1;
    t9 = I_arr(3,7)*w2;
    t10 = I_arr(3,8)*w3;
    t11 = t2+t3+t4;
    t12 = t5+t6+t7;
    t13 = t8+t9+t10;
    PZsparse N4_1 = I_arr(3,0)*wdot1+I_arr(3,1)*wdot2+I_arr(3,2)*wdot3-t12*w_aux3+t13*w_aux2;
    PZsparse N4_2 = I_arr(3,3)*wdot1+I_arr(3,4)*wdot2+I_arr(3,5)*wdot3+t11*w_aux3-t13*w_aux1;
    PZsparse N4_3 = I_arr(3,6)*wdot1+I_arr(3,7)*wdot2+I_arr(3,8)*wdot3-t11*w_aux2+t12*w_aux1;

    // joint 5
    w_new1 = cq5*w1-sq5*w3;
    w_new2 = -cq5*w3-sq5*w1;
    w_new3 = qd5+w2;

    w_aux_new1 = cq5*w_aux1-sq5*w_aux3;
    w_aux_new2 = -cq5*w_aux3-sq5*w_aux1;
    w_aux_new3 = qda5+w_aux2;

    wdot_new1 = cq5*wdot1-sq5*wdot3-qd5*(cq5*w_aux3+sq5*w_aux1);
    wdot_new2 = -cq5*wdot3-sq5*wdot1-qd5*(cq5*w_aux1-sq5*w_aux3);
    wdot_new3 = qdd5+wdot2;

    t2 = -linear_acc1;
    t3 = w_aux2*6.375E-3;
    t4 = wdot2*6.375E-3;
    t5 = w_aux1*w1*6.375E-3;
    t6 = w_aux1*w3*6.375E-3;
    t7 = w_aux3*2.0843E-1;
    t8 = wdot1*2.0843E-1;
    t9 = wdot3*2.0843E-1;
    t10 = w_aux1*w2*2.0843E-1;
    t11 = -t7;
    t12 = -t8;
    t13 = -t9;
    t14 = t3+t11;
    t16 = t2+t4+t6+t10+t13;
    t15 = t14*w2;
    t17 = linear_acc3+t5+t12+t15;
    linear_acc_new1 = -cq5*t16-sq5*t17;
    linear_acc_new2 = -cq5*t17+sq5*t16;
    linear_acc_new3 = linear_acc2+wdot1*6.375E-3-t14*w3+w_aux1*w1*2.0843E-1;

    w1 = w_new1;
    w2 = w_new2;
    w3 = w_new3;
    w_aux1 = w_aux_new1;
    w_aux2 = w_aux_new2;
    w_aux3 = w_aux_new3;
    wdot1 = wdot_new1;
    wdot2 = wdot_new2;
    wdot3 = wdot_new3;
    linear_acc1 = linear_acc_new1;
    linear_acc2 = linear_acc_new2;
    linear_acc3 = linear_acc_new3;

    t2 = com[4][0]*w_aux2;
    t3 = com[4][0]*w_aux3;
    t4 = com[4][1]*w_aux1;
    t5 = com[4][1]*w_aux3;
    t6 = com[4][2]*w_aux1;
    t7 = com[4][2]*w_aux2;
    t8 = -t4;
    t9 = -t6;
    t10 = -t7;
    t11 = t2+t8;
    t12 = t3+t9;
    t13 = t5+t10;
    PZsparse F5_1 = -mass_arr(4, 0)*(-linear_acc1+com[4][1]*wdot3-com[4][2]*wdot2+t11*w2+t12*w3);
    PZsparse F5_2 = mass_arr(4, 0)*(linear_acc2+com[4][0]*wdot3-com[4][2]*wdot1+t11*w1-t13*w3);
    PZsparse F5_3 = mass_arr(4, 0)*(linear_acc3-com[4][0]*wdot2+com[4][1]*wdot1+t12*w1+t13*w2);

    t2 = I_arr(4,0)*w1;
    t3 = I_arr(4,1)*w2;
    t4 = I_arr(4,2)*w3;
    t5 = I_arr(4,3)*w1;
    t6 = I_arr(4,4)*w2;
    t7 = I_arr(4,5)*w3;
    t8 = I_arr(4,6)*w1;
    t9 = I_arr(4,7)*w2;
    t10 = I_arr(4,8)*w3;
    t11 = t2+t3+t4;
    t12 = t5+t6+t7;
    t13 = t8+t9+t10;
    PZsparse N5_1 = I_arr(4,0)*wdot1+I_arr(4,1)*wdot2+I_arr(4,2)*wdot3-t12*w_aux3+t13*w_aux2;
    PZsparse N5_2 = I_arr(4,3)*wdot1+I_arr(4,4)*wdot2+I_arr(4,5)*wdot3+t11*w_aux3-t13*w_aux1;
    PZsparse N5_3 = I_arr(4,6)*wdot1+I_arr(4,7)*wdot2+I_arr(4,8)*wdot3-t11*w_aux2+t12*w_aux1;

    // joint 6
    w_new1 = cq6*w1+sq6*w3;
    w_new2 = cq6*w3-sq6*w1;
    w_new3 = qd6-w2;

    w_aux_new1 = cq6*w_aux1+sq6*w_aux3;
    w_aux_new2 = cq6*w_aux3-sq6*w_aux1;
    w_aux_new3 = qda6-w_aux2;

    wdot_new1 = cq6*wdot1+sq6*wdot3+qd6*(cq6*w_aux3-sq6*w_aux1);
    wdot_new2 = cq6*wdot3-sq6*wdot1-qd6*(cq6*w_aux1+sq6*w_aux3);
    wdot_new3 = qdd6-wdot2;

    t2 = -linear_acc1;
    t3 = w_aux2*1.0593E-1;
    t4 = wdot2*1.0593E-1;
    t5 = w_aux1*w1*1.0593E-1;
    t6 = w_aux1*w3*1.0593E-1;
    t7 = w_aux3*1.750499999999995E-4;
    t8 = wdot1*1.750499999999995E-4;
    t9 = wdot3*1.750499999999995E-4;
    t10 = w_aux1*w2*1.750499999999995E-4;
    t11 = -t10;
    t12 = t3+t7;
    t13 = t12*w2;
    t15 = t2+t4+t6+t9+t11;
    t14 = linear_acc3+t5+t8+t13;
    linear_acc_new1 = -cq6*t15+sq6*t14;
    linear_acc_new2 = cq6*t14+sq6*t15;
    linear_acc_new3 = -linear_acc2-wdot1*1.0593E-1+t12*w3+w_aux1*w1*1.750499999999995E-4;

    w1 = w_new1;
    w2 = w_new2;
    w3 = w_new3;
    w_aux1 = w_aux_new1;
    w_aux2 = w_aux_new2;
    w_aux3 = w_aux_new3;
    wdot1 = wdot_new1;
    wdot2 = wdot_new2;
    wdot3 = wdot_new3;
    linear_acc1 = linear_acc_new1;
    linear_acc2 = linear_acc_new2;
    linear_acc3 = linear_acc_new3;

    t2 = com[5][0]*w_aux2;
    t3 = com[5][0]*w_aux3;
    t4 = com[5][1]*w_aux1;
    t5 = com[5][1]*w_aux3;
    t6 = com[5][2]*w_aux1;
    t7 = com[5][2]*w_aux2;
    t8 = -t4;
    t9 = -t6;
    t10 = -t7;
    t11 = t2+t8;
    t12 = t3+t9;
    t13 = t5+t10;
    PZsparse F6_1 = -mass_arr(5, 0)*(-linear_acc1+com[5][1]*wdot3-com[5][2]*wdot2+t11*w2+t12*w3);
    PZsparse F6_2 = mass_arr(5, 0)*(linear_acc2+com[5][0]*wdot3-com[5][2]*wdot1+t11*w1-t13*w3);
    PZsparse F6_3 = mass_arr(5, 0)*(linear_acc3-com[5][0]*wdot2+com[5][1]*wdot1+t12*w1+t13*w2);

    t2 = I_arr(5,0)*w1;
    t3 = I_arr(5,1)*w2;
    t4 = I_arr(5,2)*w3;
    t5 = I_arr(5,3)*w1;
    t6 = I_arr(5,4)*w2;
    t7 = I_arr(5,5)*w3;
    t8 = I_arr(5,6)*w1;
    t9 = I_arr(5,7)*w2;
    t10 = I_arr(5,8)*w3;
    t11 = t2+t3+t4;
    t12 = t5+t6+t7;
    t13 = t8+t9+t10;
    PZsparse N6_1 = I_arr(5,0)*wdot1+I_arr(5,1)*wdot2+I_arr(5,2)*wdot3-t12*w_aux3+t13*w_aux2;
    PZsparse N6_2 = I_arr(5,3)*wdot1+I_arr(5,4)*wdot2+I_arr(5,5)*wdot3+t11*w_aux3-t13*w_aux1;
    PZsparse N6_3 = I_arr(5,6)*wdot1+I_arr(5,7)*wdot2+I_arr(5,8)*wdot3-t11*w_aux2+t12*w_aux1;

    // joint 7
    w_new1 = cq7*w1-sq7*w3;
    w_new2 = -cq7*w3-sq7*w1;
    w_new3 = qd7+w2;

    w_aux_new1 = cq7*w_aux1-sq7*w_aux3;
    w_aux_new2 = -cq7*w_aux3-sq7*w_aux1;
    w_aux_new3 = qda7+w_aux2;

    wdot_new1 = cq7*wdot1-sq7*wdot3-qd7*(cq7*w_aux3+sq7*w_aux1);
    wdot_new2 = -cq7*wdot3-sq7*wdot1-qd7*(cq7*w_aux1-sq7*w_aux3);
    wdot_new3 = qdd7+wdot2;

    t2 = -linear_acc1;
    t3 = w_aux3*1.0593E-1;
    t4 = wdot1*1.0593E-1;
    t5 = wdot3*1.0593E-1;
    t6 = w_aux1*w2*1.0593E-1;
    t10 = w_aux2*1.75050000000003E-4;
    t11 = wdot2*1.75050000000003E-4;
    t12 = w_aux1*w1*1.75050000000003E-4;
    t13 = w_aux1*w3*1.75050000000003E-4;
    t7 = -t3;
    t8 = -t4;
    t9 = -t5;
    t15 = -w2*(t3-t10);
    t14 = t7+t10;
    t16 = linear_acc3+t8+t12+t15;
    t17 = t2+t6+t9+t11+t13;
    linear_acc_new1 = -cq7*t17-sq7*t16;
    linear_acc_new2 = -cq7*t16+sq7*t17;
    linear_acc_new3 = linear_acc2+wdot1*1.75050000000003E-4+w_aux1*w1*1.0593E-1+w3*(t3-t10);

    w1 = w_new1;
    w2 = w_new2;
    w3 = w_new3;
    w_aux1 = w_aux_new1;
    w_aux2 = w_aux_new2;
    w_aux3 = w_aux_new3;
    wdot1 = wdot_new1;
    wdot2 = wdot_new2;
    wdot3 = wdot_new3;
    linear_acc1 = linear_acc_new1;
    linear_acc2 = linear_acc_new2;
    linear_acc3 = linear_acc_new3;

    t2 = com[6][0]*w_aux2;
    t3 = com[6][0]*w_aux3;
    t4 = com[6][1]*w_aux1;
    t5 = com[6][1]*w_aux3;
    t6 = com[6][2]*w_aux1;
    t7 = com[6][2]*w_aux2;
    t8 = -t4;
    t9 = -t6;
    t10 = -t7;
    t11 = t2+t8;
    t12 = t3+t9;
    t13 = t5+t10;
    PZsparse F7_1 = -mass_arr(6, 0)*(-linear_acc1+com[6][1]*wdot3-com[6][2]*wdot2+t11*w2+t12*w3);
    PZsparse F7_2 = mass_arr(6, 0)*(linear_acc2+com[6][0]*wdot3-com[6][2]*wdot1+t11*w1-t13*w3);
    PZsparse F7_3 = mass_arr(6, 0)*(linear_acc3-com[6][0]*wdot2+com[6][1]*wdot1+t12*w1+t13*w2);

    t2 = I_arr(6,0)*w1;
    t3 = I_arr(6,1)*w2;
    t4 = I_arr(6,2)*w3;
    t5 = I_arr(6,3)*w1;
    t6 = I_arr(6,4)*w2;
    t7 = I_arr(6,5)*w3;
    t8 = I_arr(6,6)*w1;
    t9 = I_arr(6,7)*w2;
    t10 = I_arr(6,8)*w3;
    t11 = t2+t3+t4;
    t12 = t5+t6+t7;
    t13 = t8+t9+t10;
    PZsparse N7_1 = I_arr(6,0)*wdot1+I_arr(6,1)*wdot2+I_arr(6,2)*wdot3-t12*w_aux3+t13*w_aux2;
    PZsparse N7_2 = I_arr(6,3)*wdot1+I_arr(6,4)*wdot2+I_arr(6,5)*wdot3+t11*w_aux3-t13*w_aux1;
    PZsparse N7_3 = I_arr(6,6)*wdot1+I_arr(6,7)*wdot2+I_arr(6,8)*wdot3-t11*w_aux2+t12*w_aux1;

    PZsparse f7_1 = F7_1;
    PZsparse f7_2 = F7_2;
    PZsparse f7_3 = F7_3;

    PZsparse n7_1 = N7_1+F7_3*com[6][1]-F7_2*com[6][2];
    PZsparse n7_2 = N7_2-F7_3*com[6][0]+F7_1*com[6][2];
    PZsparse n7_3 = N7_3+F7_2*com[6][0]-F7_1*com[6][1];

    PZsparse f6_1 = F6_1+cq7*f7_1-f7_2*sq7;
    PZsparse f6_2 = F6_2+f7_3;
    PZsparse f6_3 = F6_3-cq7*f7_2-f7_1*sq7;

    PZsparse  n6_1 = N6_1+f7_3*1.75050000000003E-4+F6_3*com[5][1]-F6_2*com[5][2]+cq7*f7_2*1.0593E-1+cq7*n7_1+f7_1*sq7*1.0593E-1-n7_2*sq7;
    PZsparse  n6_2 = N6_2+n7_3-F6_3*com[5][0]+F6_1*com[5][2]-cq7*f7_1*1.75050000000003E-4+f7_2*sq7*1.75050000000003E-4;
    PZsparse  n6_3 = N6_3+F6_2*com[5][0]-F6_1*com[5][1]+cq7*f7_1*1.0593E-1-cq7*n7_2-f7_2*sq7*1.0593E-1-n7_1*sq7;

    PZsparse f5_1 = F5_1+cq6*f6_1-f6_2*sq6;
    PZsparse f5_2 = F5_2-f6_3;
    PZsparse f5_3 = F5_3+cq6*f6_2+f6_1*sq6;

    PZsparse  n5_1 = N5_1-f6_3*1.0593E-1+F5_3*com[4][1]-F5_2*com[4][2]+cq6*f6_2*1.750499999999995E-4+cq6*n6_1+f6_1*sq6*1.750499999999995E-4-n6_2*sq6;
    PZsparse  n5_2 = N5_2-n6_3-F5_3*com[4][0]+F5_1*com[4][2]-cq6*f6_1*1.0593E-1+f6_2*sq6*1.0593E-1;
    PZsparse  n5_3 = N5_3+F5_2*com[4][0]-F5_1*com[4][1]-cq6*f6_1*1.750499999999995E-4+cq6*n6_2+f6_2*sq6*1.750499999999995E-4+n6_1*sq6;

    PZsparse f4_1 = F4_1+cq5*f5_1-f5_2*sq5;
    PZsparse f4_2 = F4_2+f5_3;
    PZsparse f4_3 = F4_3-cq5*f5_2-f5_1*sq5;

    PZsparse  n4_1 = N4_1+f5_3*6.375E-3+F4_3*com[3][1]-F4_2*com[3][2]+cq5*f5_2*2.0843E-1+cq5*n5_1+f5_1*sq5*2.0843E-1-n5_2*sq5;
    PZsparse  n4_2 = N4_2+n5_3-F4_3*com[3][0]+F4_1*com[3][2]-cq5*f5_1*6.375E-3+f5_2*sq5*6.375E-3;
    PZsparse  n4_3 = N4_3+F4_2*com[3][0]-F4_1*com[3][1]+cq5*f5_1*2.0843E-1-cq5*n5_2-f5_2*sq5*2.0843E-1-n5_1*sq5;

    PZsparse f3_1 = F3_1+cq4*f4_1-f4_2*sq4;
    PZsparse f3_2 = F3_2-f4_3;
    PZsparse f3_3 = F3_3+cq4*f4_2+f4_1*sq4;

    PZsparse  n3_1 = N3_1-f4_3*2.1038E-1+F3_3*com[2][1]-F3_2*com[2][2]+cq4*f4_2*6.375E-3+cq4*n4_1+f4_1*sq4*6.375E-3-n4_2*sq4;
    PZsparse  n3_2 = N3_2-n4_3-F3_3*com[2][0]+F3_1*com[2][2]-cq4*f4_1*2.1038E-1+f4_2*sq4*2.1038E-1;
    PZsparse  n3_3 = N3_3+F3_2*com[2][0]-F3_1*com[2][1]-cq4*f4_1*6.375E-3+cq4*n4_2+f4_2*sq4*6.375E-3+n4_1*sq4;

    PZsparse f2_1 = F2_1+cq3*f3_1-f3_2*sq3;
    PZsparse f2_2 = F2_2+f3_3;
    PZsparse f2_3 = F2_3-cq3*f3_2-f3_1*sq3;

    PZsparse  n2_1 = N2_1+f3_3*6.375E-3+F2_3*com[1][1]-F2_2*com[1][2]+cq3*f3_2*2.1038E-1+cq3*n3_1+f3_1*sq3*2.1038E-1-n3_2*sq3;
    PZsparse  n2_2 = N2_2+n3_3-F2_3*com[1][0]+F2_1*com[1][2]-cq3*f3_1*6.375E-3+f3_2*sq3*6.375E-3;
    PZsparse  n2_3 = N2_3+F2_2*com[1][0]-F2_1*com[1][1]+cq3*f3_1*2.1038E-1-cq3*n3_2-f3_2*sq3*2.1038E-1-n3_1*sq3;

    PZsparse f1_1 = F1_1+cq2*f2_1-f2_2*sq2;
    PZsparse f1_2 = F1_2-f2_3;
    PZsparse f1_3 = F1_3+cq2*f2_2+f2_1*sq2;

    PZsparse n1_1 = N1_1-f2_3*1.2838E-1+F1_3*com[0][1]-F1_2*com[0][2]+cq2*f2_2*5.375E-3+cq2*n2_1+f2_1*sq2*5.375E-3-n2_2*sq2;
    PZsparse n1_2 = N1_2-n2_3-F1_3*com[0][0]+F1_1*com[0][2]-cq2*f2_1*1.2838E-1+f2_2*sq2*1.2838E-1;
    PZsparse n1_3 = N1_3+F1_2*com[0][0]-F1_1*com[0][1]-cq2*f2_1*5.375E-3+cq2*n2_2+f2_2*sq2*5.375E-3+n2_1*sq2;

    u(1, 0) = n1_3;
    u(2, 0) = n2_3;
    u(3, 0) = n3_3;
    u(4, 0) = n4_3;
    u(5, 0) = n5_3;
    u(6, 0) = n6_3;
    u(7, 0) = n7_3;
}

// void KinematicsDynamics::rnea(uint s_ind,
//                               PZsparseArray& mass_arr,
//                               PZsparseArray& I_arr,
//                               PZsparseArray& u,
//                               bool setGravity) {
//     PZsparse w(3, 1);
//     PZsparse wdot(3, 1);
//     PZsparse w_aux(3, 1);
//     PZsparse linear_acc(3, 1);

//     PZsparseArray F(NUM_JOINTS, 1);
//     PZsparseArray N(NUM_JOINTS, 1);

//     if (setGravity) { // set gravity
//         // directly modify the center of the PZ instance
//         linear_acc.center(2) = gravity;
//     }

//     // RNEA forward recursion
//     for (int i = 0; i < NUM_JOINTS; i++) {
//         // NOTE:
//         // This is just a simplified implementation!!!
//         // We assume all fixed joints are at the end and the revolute joints are consecutive
//         if (axes[i] != 0) { // revolute joints
//             // line 16
//             linear_acc = traj->R_t(i, s_ind) * (linear_acc 
//                                                  + cross(wdot, trans_matrix(i, 0)) 
//                                                  + cross(w, cross(w_aux, trans_matrix(i, 0))));

//             // line 13
//             w = traj->R_t(i, s_ind) * w;
//             w.addOneDimPZ(traj->qd_des(i, s_ind), abs(axes[i]) - 1, 0);

//             // line 14
//             w_aux = traj->R_t(i, s_ind) * w_aux;

//             // line 15
//             wdot = traj->R_t(i, s_ind) * wdot;

//             PZsparse temp(3, 1); // temp = joint_vel(robot_params.q_index(i))*z(:,i)
//             temp.addOneDimPZ(traj->qd_des(i, s_ind), abs(axes[i]) - 1, 0);

//             wdot = wdot + cross(w_aux, temp);

//             wdot.addOneDimPZ(traj->qdda_des(i, s_ind), abs(axes[i]) - 1, 0);

//             // line 14
//             w_aux.addOneDimPZ(traj->qda_des(i, s_ind), abs(axes[i]) - 1, 0);
//         }
//         else { // fixed joints
//             // line 16
//             linear_acc = traj->R_t(i, s_ind) * (linear_acc 
//                                                  + cross(wdot, trans_matrix(i, 0)) 
//                                                  + cross(w, cross(w_aux, trans_matrix(i, 0))));

//             // line 13
//             w = traj->R_t(i, s_ind) * w;

//             // line 14
//             w_aux = traj->R_t(i, s_ind) * w_aux;

//             // line 15
//             wdot = traj->R_t(i, s_ind) * wdot;
//         }

//         // line 23 & 27
//         F(i, 0) = mass_arr(i, 0) * (linear_acc
//                                      + cross(wdot, com_matrix(i, 0))
//                                      + cross(w, cross(w_aux, com_matrix(i, 0))));

//         // line 29
//         N(i, 0) = I_arr(i, 0) * wdot + cross(w_aux, (I_arr(i, 0) * w));
//     }

//     PZsparse f(3, 1);
//     PZsparse n(3, 1);

//     // RNEA reverse recursion
//     for (int i = NUM_JOINTS - 1; i >= 0; i--) {
//         // line 29
//         n = N(i, 0)
//             + traj->R(i + 1, s_ind) * n
//             + cross(com_matrix(i, 0), F(i, 0))
//             + cross(trans_matrix(i + 1, 0), traj->R(i + 1, s_ind) * f);

//         // line 28
//         f = traj->R(i + 1, s_ind) * f + F(i, 0);

//         if (axes[i] != 0) {
//             u(i, s_ind) = n(abs(axes[i]) - 1, 0);

//             u(i, s_ind) = u(i, s_ind) + armature[i] * traj->qdda_des(i, s_ind);

//             u(i, s_ind) = u(i, s_ind) + damping[i] * traj->qd_des(i, s_ind);

//             // friction is directly cut on the torque limits
//         }
//     }
// }

#endif