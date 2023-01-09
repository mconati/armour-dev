#ifndef PZ_SPARSE_CPP
#define PZ_SPARSE_CPP

#include "PZsparse.h"

/*
Helper functions
*/

double getCenter(Interval a) {
    return (a.lower() + a.upper()) * 0.5;
}

double getRadius(Interval a) {
    return (a.upper() - a.lower()) * 0.5;
}

bool Monomial_sorter_degree(Monomial const& lhs, Monomial const& rhs) {
    return lhs.degree < rhs.degree;
}

/*
Initialization
*/

PZsparse::PZsparse(uint NRows_inp, uint NCols_inp) {
    NRows = NRows_inp;
    NCols = NCols_inp;
    center = Eigen::MatrixXd::Zero(NRows, NCols);
    independent = Eigen::MatrixXd::Zero(NRows, NCols);
}

PZsparse::PZsparse(const PZsparse& pz_inp) {
    NRows = pz_inp.NRows;
    NCols = pz_inp.NCols;
    center = pz_inp.center;
    polynomial = pz_inp.polynomial;
    independent = pz_inp.independent;
}

// 1x1 PZ
PZsparse::PZsparse(double center_inp) {
    NRows = 1;
    NCols = 1;
    center.resize(NRows, NCols);
    center(0) = center_inp;
    independent = Eigen::MatrixXd::Zero(NRows, NCols);
}

// NxM PZ
PZsparse::PZsparse(const Eigen::MatrixXd& center_inp) {
    NRows = center_inp.rows();
    NCols = center_inp.cols();
    center = center_inp;
    independent = Eigen::MatrixXd::Zero(NRows, NCols);
}

// 1x1 PZ
PZsparse::PZsparse(double center_inp, double uncertainty_percent) {
    NRows = 1;
    NCols = 1;
    center.resize(NRows, NCols);
    center(0) = center_inp;
    independent.resize(NRows, NCols);
    independent(0) = uncertainty_percent * fabs(center_inp);
}

// // 1x1 PZ
// PZsparse::PZsparse(Interval interval_inp) {
//     NRows = 1;
//     NCols = 1;
//     center.resize(NRows, NCols);
//     center(0) = getCenter(interval_inp);
//     independent.resize(NRows, NCols);
//     independent(0) = getRadius(interval_inp);
// }

// 1x1 PZ
PZsparse::PZsparse(double center_inp, double* coeff_inp, uint64_t degree_inp[][NUM_FACTORS * 6], uint num_monomials) {
    NRows = 1;
    NCols = 1;
    
    center.resize(NRows, NCols);
    center(0) = center_inp;

    polynomial.reserve(num_monomials);

    for (uint i = 0; i < num_monomials; i++) {
        polynomial.emplace_back(coeff_inp[i], convertDegreeToHash(degree_inp[i]));
    }

    independent = Eigen::MatrixXd::Zero(NRows, NCols);

    simplify();
}

// 1x1 PZ
PZsparse::PZsparse(double center_inp, double* coeff_inp, uint64_t degree_inp[][NUM_FACTORS * 6], uint num_monomials, Interval independent_inp) {
    NRows = 1;
    NCols = 1;
    
    center.resize(NRows, NCols);
    center(0) = center_inp;

    polynomial.reserve(num_monomials);

    for (uint i = 0; i < num_monomials; i++) {
        polynomial.emplace_back(coeff_inp[i], convertDegreeToHash(degree_inp[i]));
    }

    // assume independent_inp is centered at 0
    independent.resize(NRows, NCols);
    independent(0) = getRadius(independent_inp);

    simplify();
}

// 3x3 PZ
PZsparse::PZsparse(const double roll, const double pitch, const double yaw) {
    NRows = 3;
    NCols = 3;
    center.resize(NRows, NCols);
    
    center(0,0) = cos(pitch)*cos(yaw);
    center(0,1) = -cos(pitch)*sin(yaw);
    center(0,2) = sin(pitch);
    center(1,0) = cos(roll)*sin(yaw) + cos(yaw)*sin(pitch)*sin(roll);
    center(1,1) = cos(roll)*cos(yaw) - sin(pitch)*sin(roll)*sin(yaw);
    center(1,2) = -cos(pitch)*sin(roll);
    center(2,0) = sin(roll)*sin(yaw) - cos(roll)*cos(yaw)*sin(pitch);
    center(2,1) = cos(yaw)*sin(roll) + cos(roll)*sin(pitch)*sin(yaw);
    center(2,2) = cos(pitch)*cos(roll);

    independent = Eigen::MatrixXd::Zero(NRows, NCols);
}

// 3x3 PZ
PZsparse::PZsparse(double cos_center_inp, double* cos_coeff_inp, uint64_t cos_degree_inp[][NUM_FACTORS * 6], uint cos_num_monomials,
                   double sin_center_inp, double* sin_coeff_inp, uint64_t sin_degree_inp[][NUM_FACTORS * 6], uint sin_num_monomials,
                   const uint axis) {
    NRows = 3;
    NCols = 3;
    
    makeRotationMatrix(center, cos_center_inp, sin_center_inp, axis);

    polynomial.reserve(cos_num_monomials + sin_num_monomials);

    Eigen::MatrixXd coeff_temp;
    for (uint i = 0; i < cos_num_monomials; i++) {
        makeRotationMatrix(coeff_temp, cos_coeff_inp[i], 0, axis);
        polynomial.emplace_back(coeff_temp, convertDegreeToHash(cos_degree_inp[i]));
    }

    for (uint i = 0; i < sin_num_monomials; i++) {
        makeRotationMatrix(coeff_temp, 0, sin_coeff_inp[i], axis);
        polynomial.emplace_back(coeff_temp, convertDegreeToHash(sin_degree_inp[i]));
    }

    // assume independent_inp is centered at 0
    // makeRotationMatrix(independent, getRadius(cos_independent_inp), getRadius(sin_independent_inp), axis);
    independent = Eigen::MatrixXd::Zero(3,3);

    simplify();
}

/*
Internal functions
*/

void PZsparse::makeRotationMatrix(Eigen::MatrixXd& R, const double cosElt, const double sinElt, const uint axis, bool startFromZero) {
    if (startFromZero) {
        R = Eigen::MatrixXd::Zero(3,3);
    }
    else {
        R = Eigen::MatrixXd::Identity(3,3);
    }

    const double negSinElt = -1.0 * sinElt;

    switch (axis) {
        case 0: // fixed joints
            // NOTE:
            // This is just a simplified implementation!!!
            // The rotation matrix is just identity matrix for all fixed joints of Fetch
            // don't do anything
            return;
        case 1: // rx
            R(1,1) = cosElt;
            R(1,2) = negSinElt;
            R(2,1) = sinElt;
            R(2,2) = cosElt;
            break;
        case 2: // ry
            R(0,0) = cosElt;
            R(0,2) = sinElt;
            R(2,0) = negSinElt;
            R(2,2) = cosElt;
            break;
        case 3: // rz
            R(0,0) = cosElt;
            R(0,1) = negSinElt;
            R(1,0) = sinElt;
            R(1,1) = cosElt;
            break;
        default:
            WARNING_PRINT("Undefined axis");
            throw -1;
    }
}

bool PZsparse::checkDimensions() {
    if (center.rows() != NRows) {
        cout << "PZsparse warning: center matrix number of rows not consistent!" << endl;
        return false;
    }
    if (center.cols() != NCols) {
        cout << "PZsparse warning: center matrix number of columns not consistent!" << endl;
        return false;
    }
    if (independent.rows() != NRows) {
        cout << "PZsparse warning: independent generator matrix number of rows not consistent!" << endl;
        return false;
    }
    if (independent.cols() != NCols) {
        cout << "PZsparse warning: independent generator matrix number of columns not consistent!" << endl;
        return false;
    }
    return true;
}

void PZsparse::simplify() {
    assert(checkDimensions());

    sort(polynomial.begin(), polynomial.end(), Monomial_sorter_degree);

    Eigen::MatrixXd reduce_amount(NRows, NCols); 
    reduce_amount.setZero();

    vector<Monomial> polynomial_new;
    polynomial_new.reserve(polynomial.size());

    uint i = 0;
    while (i < polynomial.size()) {
        uint j;
        uint64_t degree = polynomial[i].degree;

        for (j = i + 1; j < polynomial.size(); j++) {
            if (polynomial[j].degree != degree) {
                break;
            }

            polynomial[i].coeff += polynomial[j].coeff;
        }

        Eigen::MatrixXd temp = polynomial[i].coeff;
        if (temp.norm() <= SIMPLIFY_THRESHOLD) {
            reduce_amount += temp.cwiseAbs();
        }
        else {
            polynomial_new.emplace_back(polynomial[i]);
        }

        i = j;
    }

    polynomial = polynomial_new;

    if (reduce_amount.norm() != 0) {
        independent = independent + reduce_amount;
    }
}

void PZsparse::reduce() {
    assert(checkDimensions());

    vector<Monomial> polynomial_new;
    polynomial_new.reserve(polynomial.size());

    for (auto it : polynomial) {
        if (it.degree <= (1 << (2 * NUM_FACTORS))) { // only dependent on k
            polynomial_new.emplace_back(it.coeff, it.degree);
        }
        else {
            independent += it.coeff.cwiseAbs();
        }
    }

    polynomial = polynomial_new;
}

MatrixXInt PZsparse::slice(const double* factor) {
    assert(checkDimensions());

    MatrixXInt res(NRows, NCols);
    Eigen::MatrixXd res_center;
    Eigen::MatrixXd res_radius;

    res_center = center;
    res_radius = independent;

    for (auto it : polynomial) {
        Eigen::MatrixXd resTemp = it.coeff;

        if (it.degree <= (1 << (2 * NUM_FACTORS))) { // only dependent on k
            convertHashToDegree(it.degree);

            for (uint j = 0; j < NUM_FACTORS; j++) {
                resTemp *= pow(factor[j], degreeArray[j]);
            }

            res_center += resTemp;
        }
        else { // this line should never be triggered if you run reduce first
            res_radius += resTemp.cwiseAbs();
        }
    }

    for (uint i = 0; i < NRows; i++) {
        for (uint j = 0; j < NCols; j++) {
            res(i,j) = Interval(res_center(i,j) - res_radius(i,j),
                                res_center(i,j) + res_radius(i,j));
        }
    }

    return res;
}

void PZsparse::slice(Eigen::Array<Eigen::MatrixXd, NUM_FACTORS, 1>& gradient, const double* factor) {
    assert(checkDimensions());

    for (uint k = 0; k < NUM_FACTORS; k++) {
        gradient[k] = Eigen::MatrixXd::Zero(NRows, NCols);
    }

    Eigen::Array<Eigen::MatrixXd, NUM_FACTORS, 1> resTemp;

    for (auto it : polynomial) {
        if (it.degree <= (1 << (2 * NUM_FACTORS))) { // only dependent on k
            for (uint k = 0; k < NUM_FACTORS; k++) {
                resTemp[k] = it.coeff;
            }

            convertHashToDegree(it.degree);

            for (uint j = 0; j < NUM_FACTORS; j++) {
                for (uint k = 0; k < NUM_FACTORS; k++) {
                    if (j == k) { // differentiate this!
                        if (degreeArray[j] == 0) { // monomial unrelated to k
                            resTemp[k] = Eigen::MatrixXd::Zero(NRows, NCols);
                        }
                        else {
                            resTemp[k] *= degreeArray[j] * pow(factor[j], degreeArray[j] - 1);
                        }
                    }
                    else {
                        resTemp[k] *= pow(factor[j], degreeArray[j]);
                    }
                }
            }

            // for (uint k = 0; k < NUM_FACTORS; k++) {
            //     gradient[k] += resTemp[k];
            // }
            gradient += resTemp;
        }
    }
}

// Interval PZsparse::toInterval() {
//     Interval res = independent + center;

//     for (auto it : polynomial) {
//         double resTemp = it.coeff;
//         res += Interval(-fabs(resTemp), fabs(resTemp));
//     }

//     return res;
// }

void PZsparse::convertHashToDegree(uint64_t degree) {
    uint64_t move_bit = 0;
    
    for (uint64_t i = 0; i < NUM_FACTORS * 6; i++) {
        degreeArray[i] = degree & DEGREE_MASK[i];    
        degree >>= MOVE_BIT_INC[i];
    }

    return;
}

uint64_t convertDegreeToHash(const uint64_t* degreeArray) {
    uint64_t degree = 0;
    uint64_t move_bit = 0;

    for (uint64_t i = 0; i < NUM_FACTORS * 6; i++) {
        if (degreeArray[i] > 1) {
            WARNING_PRINT("degree can not be larger than 1!");
            throw;
        }

        degree += (degreeArray[i] << move_bit);

        move_bit += MOVE_BIT_INC[i];
    }

    return degree;
}

std::ostream& operator<<(std::ostream& os, PZsparse& a) {
    // if (independent_only) {
    //     Interval temp = center + independent;
    //     cout << "[ " << temp.lower() << ", " << temp.upper() << " ]\n\n";
    //     return;
    // }

    os << a.center << " +...\n";

    for (auto it : a.polynomial) {
        os << '(' << it.coeff << ')';
        
        a.convertHashToDegree(it.degree);
        
        os << " * k^(";
        for (uint j = 0; j < NUM_FACTORS; j++) {
            os << a.degreeArray[j];
        }
        os << ") ";

        os << " * cosqe^(";
        for (uint j = 0; j < NUM_FACTORS; j++) {
            os << a.degreeArray[j];
        }
        os << ") ";

        os << " * sinqe^(";
        for (uint j = 0; j < NUM_FACTORS; j++) {
            os << a.degreeArray[j];
        }
        os << ") ";

        os << " * qde^(";
        for (uint j = 0; j < NUM_FACTORS; j++) {
            os << a.degreeArray[j];
        }
        os << ") ";

        os << " * qdae^(";
        for (uint j = 0; j < NUM_FACTORS; j++) {
            os << a.degreeArray[j];
        }
        os << ") ";

        os << " * qddae^(";
        for (uint j = 0; j < NUM_FACTORS; j++) {
            os << a.degreeArray[j];
        }
        os << ") ";

        os << " +...\n";
    }

    os << "[ " << -a.independent << ", " << a.independent << " ]\n\n";

    return os;
}

/*
Arithmetic
*/

PZsparse PZsparse::operator=(const double a) {
    NRows = 1;
    NCols = 1;
    center.resize(NRows, NCols);
    center(0) = a;
    polynomial.clear();
    independent = Eigen::MatrixXd::Zero(NRows, NCols);
    return *this;
}

// PZsparse PZsparse::operator=(const Interval& a) {
//     center = getCenter(a);
//     polynomial.clear();
//     independent = a - center;
//     return *this;
// }

PZsparse PZsparse::operator=(const PZsparse& a) {
    NRows = a.NRows;
    NCols = a.NCols;
    center = a.center;
    polynomial = a.polynomial;
    independent = a.independent;
    return *this;
}

PZsparse PZsparse::operator-() {
    assert(checkDimensions());

    PZsparse res(NRows, NCols);
    
    res.center = -center;

    res.polynomial.reserve(polynomial.size());

    for (auto it : polynomial) {
        res.polynomial.emplace_back(-it.coeff, it.degree);
    }

    res.independent = -independent;

    return res;
}

PZsparse PZsparse::operator+(const PZsparse& a) {
    assert(checkDimensions());

    // check if they are add-able
    assert(a.NRows == NRows || a.NCols == NCols);

    PZsparse res(NRows, NCols);

    res.center = center + a.center;

    res.polynomial.reserve(polynomial.size() + a.polynomial.size());

    res.polynomial.insert(res.polynomial.end(), polynomial.begin(), polynomial.end());

    for (auto it : a.polynomial) {
        res.polynomial.push_back(it);
    }

    res.independent = independent + a.independent;

    res.simplify();
    
    return res;
}

PZsparse PZsparse::operator+(const double a) {
    assert(checkDimensions());

    PZsparse res = *this;

    for (uint i = 0; i < NRows; i++) {
        for (uint j = 0; j < NCols; j++) {
            res.center(i,j) += a;
        }
    }

    return res;
}

PZsparse operator+(const double a, const PZsparse& b) {
    // assert(checkDimensions());

    PZsparse res = b;

    for (uint i = 0; i < b.NRows; i++) {
        for (uint j = 0; j < b.NCols; j++) {
            res.center(i,j) += a;
        }
    }

    return res;
}

PZsparse PZsparse::operator+=(const PZsparse& a) {
    assert(checkDimensions());

    assert(a.NRows == NRows || a.NCols == NCols);

    center += a.center;

    polynomial.reserve(polynomial.size() + a.polynomial.size());

    for (auto it : a.polynomial) {
        polynomial.push_back(it);
    }

    independent += a.independent;

    simplify();
    
    return *this;
}

PZsparse PZsparse::operator-(const PZsparse& a) {
    assert(checkDimensions());

    // check if they are add-able
    assert(a.NRows == NRows || a.NCols == NCols);

    PZsparse res(NRows, NCols);

    res.center = center - a.center;

    res.polynomial.reserve(polynomial.size() + a.polynomial.size());

    res.polynomial.insert(res.polynomial.end(), polynomial.begin(), polynomial.end());

    for (auto it : a.polynomial) {
        res.polynomial.emplace_back(-it.coeff, it.degree);
    }

    res.independent = independent + a.independent;

    res.simplify();
    
    return res;
}

PZsparse PZsparse::operator-(const double a) {
    assert(checkDimensions());

    PZsparse res = *this;

    for (uint i = 0; i < NRows; i++) {
        for (uint j = 0; j < NCols; j++) {
            res.center(i,j) -= a;
        }
    }

    return res;
}

PZsparse operator-(const double a, const PZsparse& b) {
    // assert(checkDimensions());

    PZsparse res = b;

    for (uint i = 0; i < b.NRows; i++) {
        for (uint j = 0; j < b.NCols; j++) {
            res.center(i,j) -= a;
        }
    }

    return res;
}

PZsparse PZsparse::operator*(const PZsparse& a) {
    assert(checkDimensions());

    assert(NCols == a.NRows);

    PZsparse res(NRows, a.NCols);

    // center * center
    res.center = center * a.center;

    cout << res.center << endl;

    res.polynomial.reserve(polynomial.size() + a.polynomial.size() + polynomial.size() * a.polynomial.size());
    // a.center * polynomial
    for (auto it : polynomial) {
        res.polynomial.emplace_back(it.coeff * a.center, it.degree);
    }

    // center * a.polynomial
    for (auto it : a.polynomial) {
        res.polynomial.emplace_back(it.coeff * center, it.degree);
    }

    // polynomial * a.polynomial (degree for each factor shouldn't be larger than 1)
    Eigen::MatrixXd reduce_amount_1 = Eigen::MatrixXd::Zero(NRows, a.NCols);

    for (auto it1 : polynomial) {
        for (auto it2 : a.polynomial) {
            Eigen::MatrixXd multiply_coeff = it1.coeff * it2.coeff;

            // Do not have to check carry
            // if we already know the maximum degree in the polynomial
            res.polynomial.emplace_back(multiply_coeff, it1.degree + it2.degree);
        }
    }

    // a.independent * (center + polynomial)
    Eigen::MatrixXd reduce_amount_2 = center.cwiseAbs();

    for (auto it : polynomial) {
        reduce_amount_2 += it.coeff.cwiseAbs();
    }

    reduce_amount_2 *= a.independent;

    // independent * (a.center + a.polynomial)
    Eigen::MatrixXd reduce_amount_3 = a.center.cwiseAbs();

    for (auto it : a.polynomial) {
        reduce_amount_3 += it.coeff.cwiseAbs();
    }

    reduce_amount_3 = independent * reduce_amount_3;

    // independent * a.independent + add reduced intervals
    Eigen::MatrixXd reduce_amount = reduce_amount_1 + reduce_amount_2 + reduce_amount_3;

    res.independent = independent * a.independent + reduce_amount;

    res.simplify();
    
    return res;
}

PZsparse PZsparse::operator*(const double a) {
    assert(checkDimensions());

    PZsparse res(NRows, NCols);

    res.center = center * a;

    res.polynomial.reserve(polynomial.size());

    for (auto it : polynomial) {
        res.polynomial.emplace_back(a * it.coeff, it.degree);
    }

    res.independent = independent * fabs(a);

    return res;
}

PZsparse operator*(const double a, const PZsparse& b) {
    // assert(b.checkDimensions());

    PZsparse res(b.NRows, b.NCols);

    res.center = b.center * a;

    res.polynomial.reserve(b.polynomial.size());

    for (auto it : b.polynomial) {
        res.polynomial.emplace_back(a * it.coeff, it.degree);
    }

    res.independent = b.independent * fabs(a);

    return res;
}

PZsparse PZsparse::operator/(const double a) {
    assert(checkDimensions());

    PZsparse res(NRows, NCols);

    res.center = center / a;

    res.polynomial.reserve(polynomial.size());

    for (auto it : polynomial) {
        res.polynomial.emplace_back(it.coeff / a, it.degree);
    }

    res.independent = independent / fabs(a);

    return res;
}

PZsparse PZsparse::transpose() {
    assert(checkDimensions());

    PZsparse res(NCols, NRows);

    res.center = center.transpose();

    res.polynomial.reserve(polynomial.size());

    for (auto it : polynomial) {
        res.polynomial.emplace_back(it.coeff.transpose(), it.degree);
    }

    res.independent = independent.transpose();

    return res;
}

#endif