#ifndef PZ_SPARSE_H
#define PZ_SPARSE_H

#include "Parameters.h"

// A specialized implementation of PZsparse
// so that more variables can be tracked
// qddae
// qdae
// qde
// sinqe
// cosqe
// k
// 6 * 7 = 42 variables in total
// qddae (maximum degree 1 bit (1))
// qdae  (maximum degree 1 bit (1))
// qde   (maximum degree 1 bit (1))
// sinqe (maximum degree 2 bits (3))
// cosqe (maximum degree 2 bits (3))
// k     (maximum degree 2 bits (3))
// 9 * 7 = 63 digits stored in one uint64_t

const uint64_t MOVE_BIT_INC[NUM_FACTORS * 6] = {2, 2, 2, 2, 2, 2, 2,
                                                2, 2, 2, 2, 2, 2, 2,
                                                2, 2, 2, 2, 2, 2, 2,
                                                1, 1, 1, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1, 1, 1};

const uint64_t DEGREE_MASK[NUM_FACTORS * 6] = {3, 3, 3, 3, 3, 3, 3,
                                               3, 3, 3, 3, 3, 3, 3,
                                               3, 3, 3, 3, 3, 3, 3,
                                               1, 1, 1, 1, 1, 1, 1,
                                               1, 1, 1, 1, 1, 1, 1,
                                               1, 1, 1, 1, 1, 1, 1};

TYPE getCenter(Interval a);

TYPE getRadius(Interval a);

typedef struct Monomial_ {
    TYPE coeff = 0;
    uint64_t degree = 0; // one hash unsigned integer for degrees for all factors, _K4_ _K3_ _K2_ _K1_

    Monomial_(TYPE coeff_inp, uint64_t degree_inp) : coeff(coeff_inp), degree(degree_inp) {};
} Monomial;

class PZsparse {
public:
    TYPE center = 0;
    vector<Monomial> polynomial;
    Interval independent = Interval(0.0, 0.0);

    uint64_t degreeArray[NUM_FACTORS * 6] = {0};

    PZsparse() {};

    PZsparse(const PZsparse& pz_inp);

    PZsparse(TYPE center_inp);

    PZsparse(TYPE center_inp, TYPE uncertainty_percent);

    PZsparse(Interval interval_inp);

    PZsparse(TYPE center_inp, Interval independent_inp);

    PZsparse(TYPE center_inp, TYPE* coeff_inp, uint64_t degree_inp[][NUM_FACTORS * 6], int num_monomials);

    PZsparse(TYPE center_inp, TYPE* coeff_inp, uint64_t degree_inp[][NUM_FACTORS * 6], int num_monomials, Interval independent_inp);

    ~PZsparse() {};

    PZsparse operator=(const TYPE& a);
    
    PZsparse operator=(const Interval& a);

    PZsparse operator=(const PZsparse& a);

    PZsparse operator-();

    PZsparse operator+(const PZsparse& a);

    PZsparse operator+(const TYPE a);

    PZsparse operator+=(const PZsparse& a);

    PZsparse operator-(const PZsparse& a);

    PZsparse operator-(const TYPE a);

    PZsparse operator*(const PZsparse& a);

    PZsparse operator*(const TYPE a);

    PZsparse operator/(const TYPE a);

    void print(bool independent_only = false);

    void simplify();

    void reduce();

    Interval slice(const TYPE* factor);

    void slice(TYPE* gradient, const TYPE* factor); // 1st-order gradient of slice

    Interval toInterval();

    void convertHashToDegree(const uint64_t degree);
};

uint64_t convertDegreeToHash(const uint64_t* degreeArray);

PZsparse operator+(const TYPE a, const PZsparse& b);

PZsparse operator-(const TYPE a, const PZsparse& b);

PZsparse operator*(const TYPE a, const PZsparse& b);

// 3 x 1 PZ vector
class vecPZsparse {
public:
    std::shared_ptr<PZsparse> elt[3];

    vecPZsparse(); // default to be all zero vector

    vecPZsparse(const TYPE* centers_inp);

    vecPZsparse(const PZsparse& elt1, const PZsparse& elt2, const PZsparse& elt3);

    ~vecPZsparse() {};

    vecPZsparse operator=(const vecPZsparse& a);

    vecPZsparse operator+(const vecPZsparse& a);

    vecPZsparse operator*(const PZsparse& a);

    void print(bool centerOnly = false);
};

vecPZsparse operator*(const TYPE a, const vecPZsparse& b);

vecPZsparse operator*(const PZsparse& a, const vecPZsparse& b);

vecPZsparse cross(const vecPZsparse& v1, const vecPZsparse& v2);

vecPZsparse cross(const vecPZsparse& v1, const TYPE* v2);

vecPZsparse cross(const TYPE* v1, const vecPZsparse& v2);

// 3 x 3 PZ matrix
// {{0, 1, 2},
//  {3, 4, 5},
//  {6, 7, 8}}
class matPZsparse {
public:
    std::shared_ptr<PZsparse> elt[9]; // row major

    matPZsparse(); // default to be identity matrix

    matPZsparse(const TYPE* centers_inp);
    
    matPZsparse(const TYPE* centers_inp, const TYPE uncertainty_percent);

    matPZsparse(const TYPE roll, const TYPE pitch, const TYPE yaw);

    matPZsparse(const PZsparse& elt1, const PZsparse& elt2, const PZsparse& elt3);

    // rotation matrix constructor, negative axis not supported now
    matPZsparse(const PZsparse& cosElt, const PZsparse& sinElt, const int axis);

    ~matPZsparse() {};

    matPZsparse operator=(const matPZsparse& a);

    void print(bool centerOnly = false);
};

vecPZsparse operator*(const matPZsparse& m, const vecPZsparse& v);

matPZsparse operator*(const matPZsparse& m1, const matPZsparse& m2);

void transpose(matPZsparse& a, const matPZsparse& b);

#endif