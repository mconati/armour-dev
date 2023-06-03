function u = kinova_rnea(q, qd, qda, qdd, mass, com, I)
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

PZsparse& qda1 = traj->qda(0, s_ind);
PZsparse& qda2 = traj->qda(1, s_ind);
PZsparse& qda3 = traj->qda(2, s_ind);
PZsparse& qda4 = traj->qda(3, s_ind);
PZsparse& qda5 = traj->qda(4, s_ind);
PZsparse& qda6 = traj->qda(5, s_ind);
PZsparse& qda7 = traj->qda(6, s_ind);

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
PZsparsewdot_new1;
PZsparsewdot_new2;
PZsparsewdot_new3;
PZsparse linear_acc_new1;
PZsparse linear_acc_new2;
PZsparse linear_acc_new3;

% joint 1
w_new3 = qd1;

w_aux_new3 = qda1;

wdot_new3 = qdd1;

linear_acc_new3 = -9.81E+2/1.0E+2;

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

t2 = com1_1*w_aux2;
t3 = com1_1*w_aux3;
t4 = com2_1*w_aux1;
t5 = com2_1*w_aux3;
t6 = com3_1*w_aux1;
t7 = com3_1*w_aux2;
t8 = -t4;
t9 = -t6;
t10 = -t7;
t11 = t2+t8;
t12 = t3+t9;
t13 = t5+t10;
F1_1 = -mass1*(-linear_acc1+com2_1*wdot3-com3_1*wdot2+t11*w2+t12*w3);
F1_2 = mass1*(linear_acc2+com1_1*wdot3-com3_1*wdot1+t11*w1-t13*w3);
F1_3 = mass1*(linear_acc3-com1_1*wdot2+com2_1*wdot1+t12*w1+t13*w2);

t2 = I1_1_1*w1;
t3 = I1_1_2*w2;
t4 = I1_1_3*w3;
t5 = I1_2_1*w1;
t6 = I1_2_2*w2;
t7 = I1_2_3*w3;
t8 = I1_3_1*w1;
t9 = I1_3_2*w2;
t10 = I1_3_3*w3;
t11 = t2+t3+t4;
t12 = t5+t6+t7;
t13 = t8+t9+t10;
N1_1 = I1_1_1*wdot1+I1_1_2*wdot2+I1_1_3*wdot3-t12*w_aux3+t13*w_aux2;
N1_2 = I1_2_1*wdot1+I1_2_2*wdot2+I1_2_3*wdot3+t11*w_aux3-t13*w_aux1;
N1_3 = I1_3_1*wdot1+I1_3_2*wdot2+I1_3_3*wdot3-t11*w_aux2+t12*w_aux1;

% joint 2
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

t2 = com1_2*w_aux2;
t3 = com1_2*w_aux3;
t4 = com2_2*w_aux1;
t5 = com2_2*w_aux3;
t6 = com3_2*w_aux1;
t7 = com3_2*w_aux2;
t8 = -t4;
t9 = -t6;
t10 = -t7;
t11 = t2+t8;
t12 = t3+t9;
t13 = t5+t10;
F2_1 = -mass2*(-linear_acc1+com2_2*wdot3-com3_2*wdot2+t11*w2+t12*w3);
F2_2 = mass2*(linear_acc2+com1_2*wdot3-com3_2*wdot1+t11*w1-t13*w3);
F2_3 = mass2*(linear_acc3-com1_2*wdot2+com2_2*wdot1+t12*w1+t13*w2);

t2 = I2_1_1*w1;
t3 = I2_1_2*w2;
t4 = I2_1_3*w3;
t5 = I2_2_1*w1;
t6 = I2_2_2*w2;
t7 = I2_2_3*w3;
t8 = I2_3_1*w1;
t9 = I2_3_2*w2;
t10 = I2_3_3*w3;
t11 = t2+t3+t4;
t12 = t5+t6+t7;
t13 = t8+t9+t10;
N2_1 = I2_1_1*wdot1+I2_1_2*wdot2+I2_1_3*wdot3-t12*w_aux3+t13*w_aux2;
N2_2 = I2_2_1*wdot1+I2_2_2*wdot2+I2_2_3*wdot3+t11*w_aux3-t13*w_aux1;
N2_3 = I2_3_1*wdot1+I2_3_2*wdot2+I2_3_3*wdot3-t11*w_aux2+t12*w_aux1;

% joint 3
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

t2 = com1_3*w_aux2;
t3 = com1_3*w_aux3;
t4 = com2_3*w_aux1;
t5 = com2_3*w_aux3;
t6 = com3_3*w_aux1;
t7 = com3_3*w_aux2;
t8 = -t4;
t9 = -t6;
t10 = -t7;
t11 = t2+t8;
t12 = t3+t9;
t13 = t5+t10;
F3_1 = -mass3*(-linear_acc1+com2_3*wdot3-com3_3*wdot2+t11*w2+t12*w3);
F3_2 = mass3*(linear_acc2+com1_3*wdot3-com3_3*wdot1+t11*w1-t13*w3);
F3_3 = mass3*(linear_acc3-com1_3*wdot2+com2_3*wdot1+t12*w1+t13*w2);

t2 = I3_1_1*w1;
t3 = I3_1_2*w2;
t4 = I3_1_3*w3;
t5 = I3_2_1*w1;
t6 = I3_2_2*w2;
t7 = I3_2_3*w3;
t8 = I3_3_1*w1;
t9 = I3_3_2*w2;
t10 = I3_3_3*w3;
t11 = t2+t3+t4;
t12 = t5+t6+t7;
t13 = t8+t9+t10;
N3_1 = I3_1_1*wdot1+I3_1_2*wdot2+I3_1_3*wdot3-t12*w_aux3+t13*w_aux2;
N3_2 = I3_2_1*wdot1+I3_2_2*wdot2+I3_2_3*wdot3+t11*w_aux3-t13*w_aux1;
N3_3 = I3_3_1*wdot1+I3_3_2*wdot2+I3_3_3*wdot3-t11*w_aux2+t12*w_aux1;

% joint 4
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

t2 = com1_4*w_aux2;
t3 = com1_4*w_aux3;
t4 = com2_4*w_aux1;
t5 = com2_4*w_aux3;
t6 = com3_4*w_aux1;
t7 = com3_4*w_aux2;
t8 = -t4;
t9 = -t6;
t10 = -t7;
t11 = t2+t8;
t12 = t3+t9;
t13 = t5+t10;
F4_1 = -mass4*(-linear_acc1+com2_4*wdot3-com3_4*wdot2+t11*w2+t12*w3);
F4_2 = mass4*(linear_acc2+com1_4*wdot3-com3_4*wdot1+t11*w1-t13*w3);
F4_3 = mass4*(linear_acc3-com1_4*wdot2+com2_4*wdot1+t12*w1+t13*w2);

t2 = I4_1_1*w1;
t3 = I4_1_2*w2;
t4 = I4_1_3*w3;
t5 = I4_2_1*w1;
t6 = I4_2_2*w2;
t7 = I4_2_3*w3;
t8 = I4_3_1*w1;
t9 = I4_3_2*w2;
t10 = I4_3_3*w3;
t11 = t2+t3+t4;
t12 = t5+t6+t7;
t13 = t8+t9+t10;
N4_1 = I4_1_1*wdot1+I4_1_2*wdot2+I4_1_3*wdot3-t12*w_aux3+t13*w_aux2;
N4_2 = I4_2_1*wdot1+I4_2_2*wdot2+I4_2_3*wdot3+t11*w_aux3-t13*w_aux1;
N4_3 = I4_3_1*wdot1+I4_3_2*wdot2+I4_3_3*wdot3-t11*w_aux2+t12*w_aux1;

% joint 5
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

t2 = com1_5*w_aux2;
t3 = com1_5*w_aux3;
t4 = com2_5*w_aux1;
t5 = com2_5*w_aux3;
t6 = com3_5*w_aux1;
t7 = com3_5*w_aux2;
t8 = -t4;
t9 = -t6;
t10 = -t7;
t11 = t2+t8;
t12 = t3+t9;
t13 = t5+t10;
F5_1 = -mass5*(-linear_acc1+com2_5*wdot3-com3_5*wdot2+t11*w2+t12*w3);
F5_2 = mass5*(linear_acc2+com1_5*wdot3-com3_5*wdot1+t11*w1-t13*w3);
F5_3 = mass5*(linear_acc3-com1_5*wdot2+com2_5*wdot1+t12*w1+t13*w2);

t2 = I5_1_1*w1;
t3 = I5_1_2*w2;
t4 = I5_1_3*w3;
t5 = I5_2_1*w1;
t6 = I5_2_2*w2;
t7 = I5_2_3*w3;
t8 = I5_3_1*w1;
t9 = I5_3_2*w2;
t10 = I5_3_3*w3;
t11 = t2+t3+t4;
t12 = t5+t6+t7;
t13 = t8+t9+t10;
N5_1 = I5_1_1*wdot1+I5_1_2*wdot2+I5_1_3*wdot3-t12*w_aux3+t13*w_aux2;
N5_2 = I5_2_1*wdot1+I5_2_2*wdot2+I5_2_3*wdot3+t11*w_aux3-t13*w_aux1;
N5_3 = I5_3_1*wdot1+I5_3_2*wdot2+I5_3_3*wdot3-t11*w_aux2+t12*w_aux1;

% joint 6
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

t2 = com1_6*w_aux2;
t3 = com1_6*w_aux3;
t4 = com2_6*w_aux1;
t5 = com2_6*w_aux3;
t6 = com3_6*w_aux1;
t7 = com3_6*w_aux2;
t8 = -t4;
t9 = -t6;
t10 = -t7;
t11 = t2+t8;
t12 = t3+t9;
t13 = t5+t10;
F6_1 = -mass6*(-linear_acc1+com2_6*wdot3-com3_6*wdot2+t11*w2+t12*w3);
F6_2 = mass6*(linear_acc2+com1_6*wdot3-com3_6*wdot1+t11*w1-t13*w3);
F6_3 = mass6*(linear_acc3-com1_6*wdot2+com2_6*wdot1+t12*w1+t13*w2);

t2 = I6_1_1*w1;
t3 = I6_1_2*w2;
t4 = I6_1_3*w3;
t5 = I6_2_1*w1;
t6 = I6_2_2*w2;
t7 = I6_2_3*w3;
t8 = I6_3_1*w1;
t9 = I6_3_2*w2;
t10 = I6_3_3*w3;
t11 = t2+t3+t4;
t12 = t5+t6+t7;
t13 = t8+t9+t10;
N6_1 = I6_1_1*wdot1+I6_1_2*wdot2+I6_1_3*wdot3-t12*w_aux3+t13*w_aux2;
N6_2 = I6_2_1*wdot1+I6_2_2*wdot2+I6_2_3*wdot3+t11*w_aux3-t13*w_aux1;
N6_3 = I6_3_1*wdot1+I6_3_2*wdot2+I6_3_3*wdot3-t11*w_aux2+t12*w_aux1;

% joint 7
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

t2 = com1_7*w_aux2;
t3 = com1_7*w_aux3;
t4 = com2_7*w_aux1;
t5 = com2_7*w_aux3;
t6 = com3_7*w_aux1;
t7 = com3_7*w_aux2;
t8 = -t4;
t9 = -t6;
t10 = -t7;
t11 = t2+t8;
t12 = t3+t9;
t13 = t5+t10;
F7_1 = -mass7*(-linear_acc1+com2_7*wdot3-com3_7*wdot2+t11*w2+t12*w3);
F7_2 = mass7*(linear_acc2+com1_7*wdot3-com3_7*wdot1+t11*w1-t13*w3);
F7_3 = mass7*(linear_acc3-com1_7*wdot2+com2_7*wdot1+t12*w1+t13*w2);

t2 = I7_1_1*w1;
t3 = I7_1_2*w2;
t4 = I7_1_3*w3;
t5 = I7_2_1*w1;
t6 = I7_2_2*w2;
t7 = I7_2_3*w3;
t8 = I7_3_1*w1;
t9 = I7_3_2*w2;
t10 = I7_3_3*w3;
t11 = t2+t3+t4;
t12 = t5+t6+t7;
t13 = t8+t9+t10;
N7_1 = I7_1_1*wdot1+I7_1_2*wdot2+I7_1_3*wdot3-t12*w_aux3+t13*w_aux2;
N7_2 = I7_2_1*wdot1+I7_2_2*wdot2+I7_2_3*wdot3+t11*w_aux3-t13*w_aux1;
N7_3 = I7_3_1*wdot1+I7_3_2*wdot2+I7_3_3*wdot3-t11*w_aux2+t12*w_aux1;

f7_1 = F7_1;
f7_2 = F7_2;
f7_3 = F7_3;

n7_1 = N7_1+F7_3*com2_7-F7_2*com3_7;
n7_2 = N7_2-F7_3*com1_7+F7_1*com3_7;
n7_3 = N7_3+F7_2*com1_7-F7_1*com2_7;

f6_1 = F6_1+cq7*f7_1-f7_2*sq7;
f6_2 = F6_2+f7_3;
f6_3 = F6_3-cq7*f7_2-f7_1*sq7;

n6_1 = N6_1+f7_3*1.75050000000003E-4+F6_3*com2_6-F6_2*com3_6+cq7*f7_2*1.0593E-1+cq7*n7_1+f7_1*sq7*1.0593E-1-n7_2*sq7;
n6_2 = N6_2+n7_3-F6_3*com1_6+F6_1*com3_6-cq7*f7_1*1.75050000000003E-4+f7_2*sq7*1.75050000000003E-4;
n6_3 = N6_3+F6_2*com1_6-F6_1*com2_6+cq7*f7_1*1.0593E-1-cq7*n7_2-f7_2*sq7*1.0593E-1-n7_1*sq7;

f5_1 = F5_1+cq6*f6_1-f6_2*sq6;
f5_2 = F5_2-f6_3;
f5_3 = F5_3+cq6*f6_2+f6_1*sq6;

n5_1 = N5_1-f6_3*1.0593E-1+F5_3*com2_5-F5_2*com3_5+cq6*f6_2*1.750499999999995E-4+cq6*n6_1+f6_1*sq6*1.750499999999995E-4-n6_2*sq6;
n5_2 = N5_2-n6_3-F5_3*com1_5+F5_1*com3_5-cq6*f6_1*1.0593E-1+f6_2*sq6*1.0593E-1;
n5_3 = N5_3+F5_2*com1_5-F5_1*com2_5-cq6*f6_1*1.750499999999995E-4+cq6*n6_2+f6_2*sq6*1.750499999999995E-4+n6_1*sq6;

f4_1 = F4_1+cq5*f5_1-f5_2*sq5;
f4_2 = F4_2+f5_3;
f4_3 = F4_3-cq5*f5_2-f5_1*sq5;

n4_1 = N4_1+f5_3*6.375E-3+F4_3*com2_4-F4_2*com3_4+cq5*f5_2*2.0843E-1+cq5*n5_1+f5_1*sq5*2.0843E-1-n5_2*sq5;
n4_2 = N4_2+n5_3-F4_3*com1_4+F4_1*com3_4-cq5*f5_1*6.375E-3+f5_2*sq5*6.375E-3;
n4_3 = N4_3+F4_2*com1_4-F4_1*com2_4+cq5*f5_1*2.0843E-1-cq5*n5_2-f5_2*sq5*2.0843E-1-n5_1*sq5;

f3_1 = F3_1+cq4*f4_1-f4_2*sq4;
f3_2 = F3_2-f4_3;
f3_3 = F3_3+cq4*f4_2+f4_1*sq4;

n3_1 = N3_1-f4_3*2.1038E-1+F3_3*com2_3-F3_2*com3_3+cq4*f4_2*6.375E-3+cq4*n4_1+f4_1*sq4*6.375E-3-n4_2*sq4;
n3_2 = N3_2-n4_3-F3_3*com1_3+F3_1*com3_3-cq4*f4_1*2.1038E-1+f4_2*sq4*2.1038E-1;
n3_3 = N3_3+F3_2*com1_3-F3_1*com2_3-cq4*f4_1*6.375E-3+cq4*n4_2+f4_2*sq4*6.375E-3+n4_1*sq4;

f2_1 = F2_1+cq3*f3_1-f3_2*sq3;
f2_2 = F2_2+f3_3;
f2_3 = F2_3-cq3*f3_2-f3_1*sq3;

n2_1 = N2_1+f3_3*6.375E-3+F2_3*com2_2-F2_2*com3_2+cq3*f3_2*2.1038E-1+cq3*n3_1+f3_1*sq3*2.1038E-1-n3_2*sq3;
n2_2 = N2_2+n3_3-F2_3*com1_2+F2_1*com3_2-cq3*f3_1*6.375E-3+f3_2*sq3*6.375E-3;
n2_3 = N2_3+F2_2*com1_2-F2_1*com2_2+cq3*f3_1*2.1038E-1-cq3*n3_2-f3_2*sq3*2.1038E-1-n3_1*sq3;

f1_1 = F1_1+cq2*f2_1-f2_2*sq2;
f1_2 = F1_2-f2_3;
f1_3 = F1_3+cq2*f2_2+f2_1*sq2;

n1_1 = N1_1-f2_3*1.2838E-1+F1_3*com2_1-F1_2*com3_1+cq2*f2_2*5.375E-3+cq2*n2_1+f2_1*sq2*5.375E-3-n2_2*sq2;
n1_2 = N1_2-n2_3-F1_3*com1_1+F1_1*com3_1-cq2*f2_1*1.2838E-1+f2_2*sq2*1.2838E-1;
n1_3 = N1_3+F1_2*com1_1-F1_1*com2_1-cq2*f2_1*5.375E-3+cq2*n2_2+f2_2*sq2*5.375E-3+n2_1*sq2;

u = zeros(7,1);
u(1, 0) = n0_3;
u(2, 0) = n1_3;
u(3, 0) = n2_3;
u(4, 0) = n3_3;
u(5, 0) = n4_3;
u(6, 0) = n5_3;
u(7, 0) = n6_3;
