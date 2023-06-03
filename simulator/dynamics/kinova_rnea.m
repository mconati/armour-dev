function u = kinova_rnea(q, qd, qda, qdd)
cq1 = cos(q(1));
cq2 = cos(q(2));
cq3 = cos(q(3));
cq4 = cos(q(4));
cq5 = cos(q(5));
cq6 = cos(q(6));
cq7 = cos(q(7));

sq1 = sin(q(1));
sq2 = sin(q(2));
sq3 = sin(q(3));
sq4 = sin(q(4));
sq5 = sin(q(5));
sq6 = sin(q(6));
sq7 = sin(q(7));

qd1 = qd(1);
qd2 = qd(2);
qd3 = qd(3);
qd4 = qd(4);
qd5 = qd(5);
qd6 = qd(6);
qd7 = qd(7);

qda1 = qda(1);
qda2 = qda(2);
qda3 = qda(3);
qda4 = qda(4);
qda5 = qda(5);
qda6 = qda(6);
qda7 = qda(7);

qdd1 = qdd(1);
qdd2 = qdd(2);
qdd3 = qdd(3);
qdd4 = qdd(4);
qdd5 = qdd(5);
qdd6 = qdd(6);
qdd7 = qdd(7);

w1 = 0;
w2 = 0;
w3 = 0;
w_aux1 = 0;
w_aux2 = 0;
w_aux3 = 0;
wdot1 = 0;
wdot2 = 0;
wdot3 = 0;
linear_acc1 = 0;
linear_acc2 = 0;
linear_acc3 = 0;

w_new1 = 0;
w_new2 = 0;
w_new3 = 0;
w_aux_new1 = 0;
w_aux_new2 = 0;
w_aux_new3 = 0;
wdot_new1 = 0;
wdot_new2 = 0;
wdot_new3 = 0;
linear_acc_new1 = 0;
linear_acc_new2 = 0;
linear_acc_new3 = 0;

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

t2 = w_aux1*7.336E-2;
t3 = w_aux2*7.336E-2;
t4 = w_aux1*1.0364E-2;
t5 = w_aux3*1.0364E-2;
t7 = w_aux2*2.3E-5;
t8 = w_aux3*2.3E-5;
t6 = -t5;
t9 = -t7;
t10 = -t8;
t11 = t3+t6;
t12 = t2+t10;
t13 = t4+t9;
F1_1 = linear_acc1*1.3773-wdot2*1.01038728E-1+wdot3*1.42743372E-2-t12*w3*1.3773-t13*w2*1.3773;
F1_2 = linear_acc2*1.3773+wdot1*1.01038728E-1-wdot3*3.16779E-5-t11*w3*1.3773+t13*w1*1.3773;
F1_3 = linear_acc3*1.3773-wdot1*1.42743372E-2+wdot2*3.16779E-5+t11*w2*1.3773+t12*w1*1.3773;

t2 = w2*4.48E-4;
t3 = w3*4.48E-4;
t4 = w1/1.0E+6;
t5 = w2/1.0E+6;
t6 = w1*4.57E-3;
t7 = w2*4.831000000000001E-3;
t8 = w3*1.409E-3;
t9 = w1*2.0E-6;
t10 = w3*2.0E-6;
t11 = t3+t4+t7;
t12 = t5+t6+t10;
t13 = t2+t8+t9;
N1_1 = wdot1*4.57E-3+wdot2/1.0E+6+wdot3*2.0E-6-t11*w_aux3+t13*w_aux2;
N1_2 = wdot1/1.0E+6+wdot2*4.831000000000001E-3+wdot3*4.48E-4-t13*w_aux1+t12*w_aux3;
N1_3 = wdot1*2.0E-6+wdot2*4.48E-4+wdot3*1.409E-3+t11*w_aux1-t12*w_aux2;

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

t2 = w_aux1*9.958E-2;
t3 = w_aux3*9.958E-2;
t5 = w_aux1*1.3278E-2;
t6 = w_aux2*1.3278E-2;
t7 = w_aux2*4.4E-5;
t8 = w_aux3*4.4E-5;
t4 = -t3;
t9 = -t7;
t10 = -t8;
t11 = t4+t6;
t12 = t2+t9;
t13 = t5+t10;
F2_1 = linear_acc1*1.1636-wdot2*1.54502808E-2+wdot3*1.15871288E-1-t12*w2*1.1636-t13*w3*1.1636;
F2_2 = linear_acc2*1.1636+wdot1*1.54502808E-2-wdot3*5.11984E-5+t12*w1*1.1636+w3*(t3-t6)*1.1636;
F2_3 = linear_acc3*1.1636-wdot1*1.15871288E-1+wdot2*5.11984E-5+t13*w1*1.1636-w2*(t3-t6)*1.1636;

t2 = w2*1.072E-3;
t3 = w1*1.1088E-2;
t4 = w3*1.1255E-2;
t5 = w2*6.910000000000002E-4;
t6 = w3*6.910000000000002E-4;
t7 = w1*5.0E-6;
t8 = w2*5.0E-6;
t9 = -t4;
t10 = -t6;
t11 = t3+t8;
t12 = t5+t9;
t13 = t2+t7+t10;
N2_1 = wdot1*1.1088E-2+wdot2*5.0E-6-t13*w_aux3+w_aux2*(t4-t5);
N2_2 = wdot1*5.0E-6+wdot2*1.072E-3-wdot3*6.910000000000002E-4+t11*w_aux3-w_aux1*(t4-t5);
N2_3 = wdot2*(-6.910000000000002E-4)+wdot3*1.1255E-2-t11*w_aux2+t13*w_aux1;

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

t2 = w_aux1*1.17892E-1;
t3 = w_aux2*1.17892E-1;
t4 = w_aux1*6.641E-3;
t5 = w_aux3*6.641E-3;
t7 = w_aux2*4.4E-5;
t8 = w_aux3*4.4E-5;
t6 = -t5;
t9 = -t7;
t10 = -t8;
t11 = t3+t6;
t12 = t2+t10;
t13 = t4+t9;
F3_1 = linear_acc1*1.1636-wdot2*1.371791312E-1+wdot3*7.7274676E-3-t12*w3*1.1636-t13*w2*1.1636;
F3_2 = linear_acc2*1.1636+wdot1*1.371791312E-1-wdot3*5.11984E-5-t11*w3*1.1636+t13*w1*1.1636;
F3_3 = linear_acc3*1.1636-wdot1*7.7274676E-3+wdot2*5.11984E-5+t11*w2*1.1636+t12*w1*1.1636;

t2 = w2*1.1127E-2;
t3 = w1*1.0932E-2;
t4 = w3*1.043E-3;
t5 = w2*6.059999999999999E-4;
t6 = w3*6.059999999999999E-4;
t7 = w1*7.0E-6;
t8 = w3*7.0E-6;
t9 = -t7;
t10 = -t8;
t11 = t2+t6;
t12 = t3+t10;
t13 = t4+t5+t9;
N3_1 = wdot1*1.0932E-2-wdot3*7.0E-6-t11*w_aux3+t13*w_aux2;
N3_2 = wdot2*1.1127E-2+wdot3*6.059999999999999E-4-t13*w_aux1+t12*w_aux3;
N3_3 = wdot1*(-7.0E-6)+wdot2*6.059999999999999E-4+wdot3*1.043E-3+t11*w_aux1-t12*w_aux2;

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

t2 = w_aux1*7.5478E-2;
t3 = w_aux3*7.5478E-2;
t4 = w_aux1*1.5006E-2;
t5 = w_aux2*1.5006E-2;
t7 = w_aux2*1.8E-5;
t8 = w_aux3*1.8E-5;
t6 = -t3;
t9 = -t7;
t10 = -t8;
t11 = t5+t6;
t12 = t2+t9;
t13 = t4+t10;
F4_1 = linear_acc1*9.302E-1-wdot2*1.39585812E-2+wdot3*7.02096356E-2-t12*w2*9.302E-1-t13*w3*9.302E-1;
F4_2 = linear_acc2*9.302E-1+wdot1*1.39585812E-2-wdot3*1.67436E-5+t12*w1*9.302E-1+w3*(t3-t5)*9.302E-1;
F4_3 = linear_acc3*9.302E-1-wdot1*7.02096356E-2+wdot2*1.67436E-5+t13*w1*9.302E-1-w2*(t3-t5)*9.302E-1;

t2 = w2/2.0E+3;
t3 = w3/2.0E+3;
t4 = w1/1.0E+6;
t5 = w2/1.0E+6;
t7 = w3*8.316E-3;
t8 = w1*8.147000000000001E-3;
t10 = w2*6.310000000000001E-4;
t6 = -t5;
t9 = -t7;
t11 = -t10;
t12 = t2+t9;
t13 = t6+t8;
t14 = t3+t4+t11;
N4_1 = wdot1*8.147000000000001E-3-wdot2/1.0E+6-t12*w_aux2+t14*w_aux3;
N4_2 = wdot1*(-1.0E-6)+wdot2*6.310000000000001E-4-wdot3/2.0E+3+t12*w_aux1-w_aux3*(t5-t8);
N4_3 = wdot2*(-5.0E-4)+wdot3*8.316E-3-t14*w_aux1+w_aux2*(t5-t8);

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

t2 = w_aux2/1.0E+6;
t3 = w_aux3/1.0E+6;
t4 = w_aux1*6.3883E-2;
t5 = w_aux2*6.3883E-2;
t6 = w_aux1*9.431999999999999E-3;
t7 = w_aux3*9.431999999999999E-3;
t8 = -t7;
t9 = t3+t4;
t10 = t2+t6;
t11 = t5+t8;
F5_1 = linear_acc1*6.781E-1-wdot2*4.33190623E-2+wdot3*6.3958392E-3-t9*w3*6.781E-1-t10*w2*6.781E-1;
F5_2 = linear_acc2*6.781E-1+wdot1*4.33190623E-2+wdot3*6.781E-7+t10*w1*6.781E-1-t11*w3*6.781E-1;
F5_3 = linear_acc3*6.781E-1-wdot1*6.3958392E-3-wdot2*6.781E-7+t9*w1*6.781E-1+t11*w2*6.781E-1;

t2 = w2*2.56E-4;
t3 = w3*2.56E-4;
t4 = w2*1.607E-3;
t5 = w3*3.99E-4;
t6 = t3+t4;
t7 = t2+t5;
N5_1 = wdot1*1.596E-3-t6*w_aux3+t7*w_aux2;
N5_2 = wdot2*1.607E-3+wdot3*2.56E-4-t7*w_aux1+w_aux3*w1*1.596E-3;
N5_3 = wdot2*2.56E-4+wdot3*3.99E-4+t6*w_aux1-w_aux2*w1*1.596E-3;

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

t2 = w_aux1*9.65E-3;
t3 = w_aux2*9.65E-3;
t4 = w_aux2/1.0E+6;
t5 = w_aux3/1.0E+6;
t7 = w_aux1*4.5483E-2;
t8 = w_aux3*4.5483E-2;
t6 = t2+t5;
t9 = -t8;
t10 = t4+t7;
t11 = t3+t9;
F6_1 = linear_acc1*6.781E-1-wdot2*6.543665E-3+wdot3*3.08420223E-2-t6*w3*6.781E-1-t10*w2*6.781E-1;
F6_2 = linear_acc2*6.781E-1+wdot1*6.543665E-3+wdot3*6.781E-7+t10*w1*6.781E-1-t11*w3*6.781E-1;
F6_3 = linear_acc3*6.781E-1-wdot1*3.08420223E-2-wdot2*6.781E-7+t6*w1*6.781E-1+t11*w2*6.781E-1;

t2 = w2*4.1E-4;
t3 = w3*1.641E-3;
t5 = w2*2.78E-4;
t6 = w3*2.78E-4;
t4 = -t3;
t7 = -t6;
t8 = t2+t7;
t9 = t4+t5;
N6_1 = wdot1*1.641E-3-t8*w_aux3+w_aux2*(t3-t5);
N6_2 = wdot2*4.1E-4-wdot3*2.78E-4+w_aux3*w1*1.641E-3-w_aux1*(t3-t5);
N6_3 = wdot2*(-2.78E-4)+wdot3*1.641E-3+t8*w_aux1-w_aux2*w1*1.641E-3;

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

t2 = w_aux2*2.81E-4;
t3 = w_aux3*2.81E-4;
t5 = w_aux1*2.9798E-2;
t6 = w_aux2*2.9798E-2;
t7 = w_aux1*1.1402E-2;
t8 = w_aux3*1.1402E-2;
t4 = -t2;
t9 = t3+t5;
t11 = t6+t8;
t10 = t4+t7;
F7_1 = linear_acc1/2.0-wdot2*1.4899E-2-wdot3*5.701E-3-(t9*w3)/2.0-(w2*(t2-t7))/2.0;
F7_2 = linear_acc2/2.0+wdot1*1.4899E-2+wdot3*1.405E-4-(t11*w3)/2.0+(w1*(t2-t7))/2.0;
F7_3 = linear_acc3/2.0+wdot1*5.701E-3-wdot2*1.405E-4+(t9*w1)/2.0+(t11*w2)/2.0;

t2 = w3*6.089999999999999E-4;
t3 = w1*5.87E-4;
t4 = w2*3.69E-4;
t5 = w2*1.18E-4;
t6 = w3*1.18E-4;
t9 = w1*3.0E-6;
t10 = w2*3.0E-6;
t11 = w3*3.0E-6;
t12 = wdot1*3.0E-6;
t7 = -t5;
t8 = -t6;
t15 = t3+t10+t11;
t13 = t2+t7+t9;
t14 = t4+t8+t9;
N7_1 = wdot1*5.87E-4+wdot2*3.0E-6+wdot3*3.0E-6+t13*w_aux2-t14*w_aux3;
N7_2 = t12+wdot2*3.69E-4-wdot3*1.18E-4-t13*w_aux1+t15*w_aux3;
N7_3 = t12-wdot2*1.18E-4+wdot3*6.089999999999999E-4+t14*w_aux1-t15*w_aux2;

f7_1 = F7_1;
f7_2 = F7_2;
f7_3 = F7_3;

n7_1 = F7_2*2.9798E-2+F7_3*1.1402E-2+N7_1;
n7_2 = F7_1*(-2.9798E-2)-F7_3*2.81E-4+N7_2;
n7_3 = F7_1*(-1.1402E-2)+F7_2*2.81E-4+N7_3;

f6_1 = F6_1+cq7*f7_1-f7_2*sq7;
f6_2 = F6_2+f7_3;
f6_3 = F6_3-cq7*f7_2-f7_1*sq7;

n6_1 = F6_2*9.65E-3-F6_3*4.5483E-2+N6_1+f7_3*1.75050000000003E-4+cq7*f7_2*1.0593E-1+cq7*n7_1+f7_1*sq7*1.0593E-1-n7_2*sq7;
n6_2 = F6_1*(-9.65E-3)-F6_3/1.0E+6+N6_2+n7_3-cq7*f7_1*1.75050000000003E-4+f7_2*sq7*1.75050000000003E-4;
n6_3 = F6_1*4.5483E-2+F6_2/1.0E+6+N6_3+cq7*f7_1*1.0593E-1-cq7*n7_2-f7_2*sq7*1.0593E-1-n7_1*sq7;

f5_1 = F5_1+cq6*f6_1-f6_2*sq6;
f5_2 = F5_2-f6_3;
f5_3 = F5_3+cq6*f6_2+f6_1*sq6;

n5_1 = F5_2*6.3883E-2-F5_3*9.431999999999999E-3+N5_1-f6_3*1.0593E-1+cq6*f6_2*1.750499999999995E-4+cq6*n6_1+f6_1*sq6*1.750499999999995E-4-n6_2*sq6;
n5_2 = F5_1*(-6.3883E-2)-F5_3/1.0E+6+N5_2-n6_3-cq6*f6_1*1.0593E-1+f6_2*sq6*1.0593E-1;
n5_3 = F5_1*9.431999999999999E-3+F5_2/1.0E+6+N5_3-cq6*f6_1*1.750499999999995E-4+cq6*n6_2+f6_2*sq6*1.750499999999995E-4+n6_1*sq6;

f4_1 = F4_1+cq5*f5_1-f5_2*sq5;
f4_2 = F4_2+f5_3;
f4_3 = F4_3-cq5*f5_2-f5_1*sq5;

n4_1 = F4_2*1.5006E-2-F4_3*7.5478E-2+N4_1+f5_3*6.375E-3+cq5*f5_2*2.0843E-1+cq5*n5_1+f5_1*sq5*2.0843E-1-n5_2*sq5;
n4_2 = F4_1*(-1.5006E-2)+F4_3*1.8E-5+N4_2+n5_3-cq5*f5_1*6.375E-3+f5_2*sq5*6.375E-3;
n4_3 = F4_1*7.5478E-2-F4_2*1.8E-5+N4_3+cq5*f5_1*2.0843E-1-cq5*n5_2-f5_2*sq5*2.0843E-1-n5_1*sq5;

f3_1 = F3_1+cq4*f4_1-f4_2*sq4;
f3_2 = F3_2-f4_3;
f3_3 = F3_3+cq4*f4_2+f4_1*sq4;

n3_1 = F3_2*1.17892E-1-F3_3*6.641E-3+N3_1-f4_3*2.1038E-1+cq4*f4_2*6.375E-3+cq4*n4_1+f4_1*sq4*6.375E-3-n4_2*sq4;
n3_2 = F3_1*(-1.17892E-1)+F3_3*4.4E-5+N3_2-n4_3-cq4*f4_1*2.1038E-1+f4_2*sq4*2.1038E-1;
n3_3 = F3_1*6.641E-3-F3_2*4.4E-5+N3_3-cq4*f4_1*6.375E-3+cq4*n4_2+f4_2*sq4*6.375E-3+n4_1*sq4;

f2_1 = F2_1+cq3*f3_1-f3_2*sq3;
f2_2 = F2_2+f3_3;
f2_3 = F2_3-cq3*f3_2-f3_1*sq3;

n2_1 = F2_2*1.3278E-2-F2_3*9.958E-2+N2_1+f3_3*6.375E-3+cq3*f3_2*2.1038E-1+cq3*n3_1+f3_1*sq3*2.1038E-1-n3_2*sq3;
n2_2 = F2_1*(-1.3278E-2)+F2_3*4.4E-5+N2_2+n3_3-cq3*f3_1*6.375E-3+f3_2*sq3*6.375E-3;
n2_3 = F2_1*9.958E-2-F2_2*4.4E-5+N2_3+cq3*f3_1*2.1038E-1-cq3*n3_2-f3_2*sq3*2.1038E-1-n3_1*sq3;

f1_1 = F1_1+cq2*f2_1-f2_2*sq2;
f1_2 = F1_2-f2_3;
f1_3 = F1_3+cq2*f2_2+f2_1*sq2;

n1_1 = F1_2*7.336E-2-F1_3*1.0364E-2+N1_1-f2_3*1.2838E-1+cq2*f2_2*5.375E-3+cq2*n2_1+f2_1*sq2*5.375E-3-n2_2*sq2;
n1_2 = F1_1*(-7.336E-2)+F1_3*2.3E-5+N1_2-n2_3-cq2*f2_1*1.2838E-1+f2_2*sq2*1.2838E-1;
n1_3 = F1_1*1.0364E-2-F1_2*2.3E-5+N1_3-cq2*f2_1*5.375E-3+cq2*n2_2+f2_2*sq2*5.375E-3+n2_1*sq2;

u = zeros(7,1);
u(1) = n1_3;
u(2) = n2_3;
u(3) = n3_3;
u(4) = n4_3;
u(5) = n5_3;
u(6) = n6_3;
u(7) = n7_3;
