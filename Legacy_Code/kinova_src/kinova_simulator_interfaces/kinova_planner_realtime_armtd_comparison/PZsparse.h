#ifndef PZ_SPARSE_H
#define PZ_SPARSE_H

#include "Parameters.h"

// A specialized implementation of PZsparse
// so that more variables can be tracked
// sinqe
// cosqe
// k
// 3 * 7 = 21 variables in total
// sinqe (maximum degree 1 bits (1))
// cosqe (maximum degree 1 bits (1))
// k     (maximum degree 1 bits (1))
// 3 * 7 = 21 digits stored in one unsigned int

const unsigned int MOVE_BIT_INC[NUM_FACTORS * 3] = {1, 1, 1, 1, 1, 1, 1,
                                                    1, 1, 1, 1, 1, 1, 1,
                                                    1, 1, 1, 1, 1, 1, 1};

const unsigned int DEGREE_MASK[NUM_FACTORS * 3] = {1, 1, 1, 1, 1, 1, 1,
                                                   1, 1, 1, 1, 1, 1, 1,
                                                   1, 1, 1, 1, 1, 1, 1};

TYPE getCenter(Interval a);

TYPE getRadius(Interval a);

typedef struct Monomial_ {
    TYPE coeff = 0;
    unsigned int degree = 0; // MAXIMUM_DEGREE hash unsigned integer for degrees for all factors, _K4_ _K3_ _K2_ _K1_

    Monomial_(TYPE coeff_inp, unsigned int degree_inp) : coeff(coeff_inp), degree(degree_inp) {};
} Monomial;

class PZsparse {
public:
    TYPE center = 0;
    vector<Monomial> polynomial;
    Interval independent = Interval(0.0, 0.0);

    unsigned int degreeArray[NUM_FACTORS * 3] = {0};

    PZsparse() {};

    PZsparse(const PZsparse& pz_inp);

    PZsparse(TYPE center_inp);

    PZsparse(TYPE center_inp, TYPE uncertainty_percent);

    PZsparse(Interval interval_inp);

    PZsparse(TYPE center_inp, Interval independent_inp);

    PZsparse(TYPE center_inp, TYPE* coeff_inp, unsigned int degree_inp[][NUM_FACTORS * 3], int num_monomials);

    PZsparse(TYPE center_inp, TYPE* coeff_inp, unsigned int degree_inp[][NUM_FACTORS * 3], int num_monomials, Interval independent_inp);

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

    void print();

    void simplify();

    void reduce();

    Interval slice(const TYPE* factor);

    void slice(TYPE* gradient, const TYPE* factor); // 1st-order gradient of slice

    void convertHashToDegree(const unsigned int degree);
};

unsigned int convertDegreeToHash(const unsigned int* degreeArray);

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

    void print();
};

vecPZsparse operator*(const TYPE a, const vecPZsparse& b);

vecPZsparse operator*(const PZsparse& a, const vecPZsparse& b);

vecPZsparse cross(const vecPZsparse& v1, const vecPZsparse& v2);

vecPZsparse cross(const vecPZsparse& v1, const TYPE* v2);

vecPZsparse cross(const TYPE* v1, const vecPZsparse& v2);

// 3 x 3 PZ matrix
class matPZsparse {
public:
    std::shared_ptr<PZsparse> elt[9]; // row major

    matPZsparse(); // default to be identity matrix

    matPZsparse(const TYPE* centers_inp);
    
    matPZsparse(const TYPE* centers_inp, const TYPE uncertainty_percent);

    matPZsparse(const TYPE roll, const TYPE pitch, const TYPE yaw);

    matPZsparse(const PZsparse& elt1, const PZsparse& elt2, const PZsparse& elt3);

    matPZsparse(const PZsparse& cosElt, const PZsparse& sinElt, const int axes);

    ~matPZsparse() {};

    matPZsparse operator=(const matPZsparse& a);

    void print();
};

vecPZsparse operator*(const matPZsparse m, const vecPZsparse& v);

matPZsparse operator*(const matPZsparse& m1, const matPZsparse& m2);

void transpose(matPZsparse& a, const matPZsparse& b);

#endif
