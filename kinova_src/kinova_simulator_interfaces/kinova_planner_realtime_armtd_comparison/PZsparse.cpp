#ifndef PZ_SPARSE_CPP
#define PZ_SPARSE_CPP

#include "PZsparse.h"

TYPE getCenter(Interval a) {
    return (a.lower() + a.upper()) * 0.5;
}

TYPE getRadius(Interval a) {
    return (a.upper() - a.lower()) * 0.5;
}

bool Monomial_sorter_degree(Monomial_ const& lhs, Monomial_ const& rhs) {
    return lhs.degree < rhs.degree;
}

bool Monomial_sorter_coeff(Monomial_ const& lhs, Monomial_ const& rhs) {
    return fabs(lhs.degree) > fabs(rhs.degree);
}

PZsparse::PZsparse(const PZsparse& pz_inp) {
    center = pz_inp.center;
    polynomial = pz_inp.polynomial;
    independent = pz_inp.independent;
}

PZsparse::PZsparse(TYPE center_inp) {
    center = center_inp;
}

PZsparse::PZsparse(TYPE center_inp, TYPE uncertainty_percent) {
    center = center_inp;
    if (center > 0) {
        independent = Interval(-uncertainty_percent * center, uncertainty_percent * center);
    }
    else {
        independent = Interval(uncertainty_percent * center, -uncertainty_percent * center);
    }
}

PZsparse::PZsparse(Interval interval_inp) {
    center = getCenter(interval_inp);
    independent = interval_inp - center;
}

PZsparse::PZsparse(TYPE center_inp, Interval independent_inp) {
    center = center_inp;
    independent = independent_inp;
}

PZsparse::PZsparse(TYPE center_inp, TYPE* coeff_inp, unsigned int degree_inp[][NUM_FACTORS * 3], int num_monomials) {
    center = center_inp;

    polynomial.reserve(num_monomials);

    for (int i = 0; i < num_monomials; i++) {
        polynomial.emplace_back(coeff_inp[i], convertDegreeToHash(degree_inp[i]));
    }

    independent = Interval(0.0, 0.0);

    simplify();
}

PZsparse::PZsparse(TYPE center_inp, TYPE* coeff_inp, unsigned int degree_inp[][NUM_FACTORS * 3], int num_monomials, Interval independent_inp) {
    center = center_inp;

    polynomial.reserve(num_monomials);

    for (int i = 0; i < num_monomials; i++) {
        polynomial.emplace_back(coeff_inp[i], convertDegreeToHash(degree_inp[i]));
    }

    independent = independent_inp;

    simplify();
}

void PZsparse::simplify() {
    sort(polynomial.begin(), polynomial.end(), Monomial_sorter_degree);

    TYPE reduce_amount = 0; 

    vector<Monomial> polynomial_new;
    polynomial_new.reserve(polynomial.size());

    int i = 0;
    while (i < polynomial.size()) {
        int j;
        unsigned int degree = polynomial[i].degree;

        for (j = i + 1; j < polynomial.size(); j++) {
            if (polynomial[j].degree != degree) {
                break;
            }

            polynomial[i].coeff += polynomial[j].coeff;
        }

        TYPE temp = fabs(polynomial[i].coeff);
        if (temp <= SIMPLIFY_THRESHOLD) {
            reduce_amount += temp;
        }
        else {
            polynomial_new.emplace_back(polynomial[i]);
        }

        i = j;
    }

    polynomial = polynomial_new;

    if (reduce_amount != 0) {
        independent = independent + Interval(-reduce_amount, reduce_amount);
    }
}

void PZsparse::print() {
    cout << center << " +...\n";

    for (auto it : polynomial) {
        cout << '(' << it.coeff << ')';
        
        convertHashToDegree(it.degree);
        
        for (int j = 0; j < NUM_FACTORS * 3; j++) {
            cout << " * k" << j << " ^ " << degreeArray[j];
        }

        cout << " +...\n";
    }

    cout << "[ " << independent.lower() << ", " << independent.upper() << " ]\n\n";
}

void PZsparse::reduce() {
    vector<Monomial> polynomial_new;
    polynomial_new.reserve(polynomial.size());

    for (auto it : polynomial) {
        if (it.degree <= (1 << (NUM_FACTORS))) { // only dependent on k
            polynomial_new.emplace_back(it.coeff, it.degree);
        }
        else {
            independent += Interval(-fabs(it.coeff), fabs(it.coeff));
        }
    }

    polynomial = polynomial_new;
}

Interval PZsparse::slice(const TYPE* factor) {
    Interval res = independent + center;

    for (auto it : polynomial) {
        TYPE resTemp = it.coeff;

        if (it.degree <= (1 << (NUM_FACTORS))) { // only dependent on k
            convertHashToDegree(it.degree);

            for (int j = 0; j < NUM_FACTORS; j++) {
                resTemp *= pow(factor[j], degreeArray[j]);
            }

            res += resTemp;
        }
        else { // this line should never be triggered if you run reduce first
            res += Interval(-fabs(resTemp), fabs(resTemp));
        }
    }

    return res;
}

void PZsparse::slice(TYPE* gradient, const TYPE* factor) {
    std::memset(gradient, 0, NUM_FACTORS * sizeof(TYPE));

    TYPE resTemp[NUM_FACTORS];

    for (auto it : polynomial) {
        if (it.degree <= (1 << (NUM_FACTORS))) { // only dependent on k
            for (int k = 0; k < NUM_FACTORS; k++) {
                resTemp[k] = it.coeff;
            }

            convertHashToDegree(it.degree);

            for (int j = 0; j < NUM_FACTORS; j++) {
                for (int k = 0; k < NUM_FACTORS; k++) {
                    if (j == k) { // differentiate this!
                        if (degreeArray[j] == 0) { // monomial unrelated to k
                            resTemp[k] = 0;
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

            for (int k = 0; k < NUM_FACTORS; k++) {
                gradient[k] += resTemp[k];
            }
        }
    }
}

PZsparse PZsparse::operator=(const TYPE& a) {
    center = a;
    polynomial.clear();
    independent = Interval(0,0);
    return *this;
}

PZsparse PZsparse::operator=(const Interval& a) {
    center = getCenter(a);
    polynomial.clear();
    independent = a - center;
    return *this;
}

PZsparse PZsparse::operator=(const PZsparse& a) {
    center = a.center;
    polynomial = a.polynomial;
    independent = a.independent;
    return *this;
}

PZsparse PZsparse::operator-() {
    PZsparse res;

    res.center = -center;

    res.polynomial.reserve(polynomial.size());

    for (auto it : polynomial) {
        res.polynomial.emplace_back(-it.coeff, it.degree);
    }

    res.independent = -independent;

    return res;
}

PZsparse PZsparse::operator+(const PZsparse& a) {
    PZsparse res;

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

PZsparse PZsparse::operator+(const TYPE a) {
    PZsparse res = *this;

    res.center += a;

    return res;
}

PZsparse operator+(const TYPE a, const PZsparse& b) {
    PZsparse res = b;

    res.center += a;

    return res;
}

PZsparse PZsparse::operator+=(const PZsparse& a) {
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
    PZsparse res;

    res.center = center - a.center;

    res.polynomial.reserve(polynomial.size() + a.polynomial.size());

    res.polynomial.insert(res.polynomial.end(), polynomial.begin(), polynomial.end());

    for (auto it : a.polynomial) {
        res.polynomial.emplace_back(-it.coeff, it.degree);
    }

    res.independent = independent - a.independent;

    res.simplify();
    
    return res;
}

PZsparse PZsparse::operator-(const TYPE a) {
    PZsparse res = *this;

    res.center -= a;

    return res;
}

PZsparse operator-(const TYPE a, const PZsparse& b){
    PZsparse res = b;

    res.center -= a;

    return res;
}

PZsparse PZsparse::operator*(const PZsparse& a) {
    PZsparse res;

    // center * center
    res.center = center * a.center;

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
    TYPE reduce_amount_1 = 0;

    for (auto it1 : polynomial) {
        for (auto it2 : a.polynomial) {
            unsigned int multiply_degree = it1.degree + it2.degree;
            TYPE multiply_coeff = it1.coeff * it2.coeff;

            // Do not have to check carry
            // if we already know the maximum degree in the polynomial
            res.polynomial.emplace_back(multiply_coeff, multiply_degree);
        }
    }

    // a.independent * (center + polynomial)
    TYPE reduce_amount_2 = fabs(center);

    for (auto it : polynomial) {
        reduce_amount_2 += fabs(it.coeff);
    }

    reduce_amount_2 *= a.independent.upper();

    // independent * (a.center + a.polynomial)
    TYPE reduce_amount_3 = fabs(a.center);

    for (auto it : a.polynomial) {
        reduce_amount_3 += fabs(it.coeff);
    }

    reduce_amount_3 *= independent.upper();

    // compute independent interval
    TYPE reduce_amount = reduce_amount_1 + reduce_amount_2 + reduce_amount_3;

    res.independent = independent * a.independent + Interval(-reduce_amount, reduce_amount);

    res.simplify();
    
    return res;
}

PZsparse PZsparse::operator*(const TYPE a) {
    PZsparse res;

    res.center = center * a;

    res.polynomial.reserve(polynomial.size());

    for (auto it : polynomial) {
        res.polynomial.emplace_back(a * it.coeff, it.degree);
    }

    res.independent = independent * a;

    return res;
}

PZsparse operator*(const TYPE a, const PZsparse& b) {
    PZsparse res;

    res.center = b.center * a;

    res.polynomial.reserve(b.polynomial.size());

    for (auto it : b.polynomial) {
        res.polynomial.emplace_back(a * it.coeff, it.degree);
    }

    res.independent = b.independent * a;

    return res;
}

PZsparse PZsparse::operator/(const TYPE a) {
    PZsparse res;

    res.center = center / a;

    res.polynomial.reserve(polynomial.size());

    for (auto it : polynomial) {
        res.polynomial.emplace_back(it.coeff / a, it.degree);
    }

    res.independent = independent / a;

    return res;
}

void PZsparse::convertHashToDegree(unsigned int degree) {
    unsigned int move_bit = 0;
    
    for (unsigned int i = 0; i < NUM_FACTORS* 3; i++) {
        degreeArray[i] = degree & DEGREE_MASK[i];    
        degree >>= MOVE_BIT_INC[i];
    }

    return;
}

unsigned int convertDegreeToHash(const unsigned int* degreeArray) {
    unsigned int degree = 0;
    unsigned int move_bit = 0;

    for (unsigned int i = 0; i < NUM_FACTORS* 3; i++) {
        if (degreeArray[i] > 1) {
            WARNING_PRINT("degree can not be larger than 1!");
            throw;
        }

        degree += (degreeArray[i] << move_bit);

        move_bit += MOVE_BIT_INC[i];
    }

    return degree;
}

vecPZsparse::vecPZsparse() {
    for (int i = 0; i < 3; i++) {
        elt[i] = std::make_shared<PZsparse>(0.0);
    }
}

vecPZsparse::vecPZsparse(const TYPE* centers_inp) {
    for (int i = 0; i < 3; i++) {
        elt[i] = std::make_shared<PZsparse>(centers_inp[i]);
    }
}

vecPZsparse::vecPZsparse(const PZsparse& elt1, const PZsparse& elt2, const PZsparse& elt3) {
    elt[0] = std::make_shared<PZsparse>(elt1);
    elt[1] = std::make_shared<PZsparse>(elt2);
    elt[2] = std::make_shared<PZsparse>(elt3);
}

void vecPZsparse::print() {
    cout << "{\n";

    for (int i = 0; i < 3; i++) {
        if (elt[i] == nullptr) {
            cout << "0\n\n";
        }
        else {
            elt[i]->print();
        }
    }

    cout << "}\n\n";
}

vecPZsparse vecPZsparse::operator=(const vecPZsparse& a) {
    for (int i = 0; i < 3; i++) {
        elt[i] = a.elt[i];
    }

    return *this;
}

vecPZsparse vecPZsparse::operator+(const vecPZsparse& a) {
    vecPZsparse res;

    res.elt[0] = std::make_shared<PZsparse>(*(elt[0]) + *(a.elt[0]));
    res.elt[1] = std::make_shared<PZsparse>(*(elt[1]) + *(a.elt[1]));
    res.elt[2] = std::make_shared<PZsparse>(*(elt[2]) + *(a.elt[2]));

    return res;
}

vecPZsparse vecPZsparse::operator*(const PZsparse& a) {
    vecPZsparse res;

    res.elt[0] = std::make_shared<PZsparse>(*(elt[0]) * a);
    res.elt[1] = std::make_shared<PZsparse>(*(elt[1]) * a);
    res.elt[2] = std::make_shared<PZsparse>(*(elt[2]) * a);

    return res;
}

vecPZsparse operator*(const TYPE a, const vecPZsparse& b) {
    vecPZsparse res;

    res.elt[0] = std::make_shared<PZsparse>(a * *(b.elt[0]));
    res.elt[1] = std::make_shared<PZsparse>(a * *(b.elt[1]));
    res.elt[2] = std::make_shared<PZsparse>(a * *(b.elt[2]));

    return res;
}

vecPZsparse operator*(const PZsparse& a, const vecPZsparse& b) {
    vecPZsparse res;

    res.elt[0] = std::make_shared<PZsparse>(*(b.elt[0]) * a);
    res.elt[1] = std::make_shared<PZsparse>(*(b.elt[1]) * a);
    res.elt[2] = std::make_shared<PZsparse>(*(b.elt[2]) * a);

    return res;
}

vecPZsparse cross(const vecPZsparse& v1, const vecPZsparse& v2) {
    vecPZsparse res;

    res.elt[0] = std::make_shared<PZsparse>( *(v1.elt[1]) * *(v2.elt[2]) - *(v1.elt[2]) * *(v2.elt[1]) );
    res.elt[1] = std::make_shared<PZsparse>( *(v1.elt[2]) * *(v2.elt[0]) - *(v1.elt[0]) * *(v2.elt[2]) );
    res.elt[2] = std::make_shared<PZsparse>( *(v1.elt[0]) * *(v2.elt[1]) - *(v1.elt[1]) * *(v2.elt[0]) );

    return res;
}

vecPZsparse cross(const vecPZsparse& v1, const TYPE* v2) {
    vecPZsparse res;

    res.elt[0] = std::make_shared<PZsparse>( *(v1.elt[1]) * v2[2] - *(v1.elt[2]) * v2[1] );
    res.elt[1] = std::make_shared<PZsparse>( *(v1.elt[2]) * v2[0] - *(v1.elt[0]) * v2[2] );
    res.elt[2] = std::make_shared<PZsparse>( *(v1.elt[0]) * v2[1] - *(v1.elt[1]) * v2[0] );

    return res;
}

vecPZsparse cross(const TYPE* v1, const vecPZsparse& v2) {
    vecPZsparse res;

    res.elt[0] = std::make_shared<PZsparse>( v1[1] * *(v2.elt[2]) - v1[2] * *(v2.elt[1]) );
    res.elt[1] = std::make_shared<PZsparse>( v1[2] * *(v2.elt[0]) - v1[0] * *(v2.elt[2]) );
    res.elt[2] = std::make_shared<PZsparse>( v1[0] * *(v2.elt[1]) - v1[1] * *(v2.elt[0]) );

    return res;
}

matPZsparse::matPZsparse() {
    elt[0] = std::make_shared<PZsparse>(1.0);
    elt[1] = std::make_shared<PZsparse>();
    elt[2] = std::make_shared<PZsparse>();
    elt[3] = std::make_shared<PZsparse>();
    elt[4] = std::make_shared<PZsparse>(1.0);
    elt[5] = std::make_shared<PZsparse>();
    elt[6] = std::make_shared<PZsparse>();
    elt[7] = std::make_shared<PZsparse>();
    elt[8] = std::make_shared<PZsparse>(1.0);
}

matPZsparse::matPZsparse(const TYPE* centers_inp) {
    for (int i = 0; i < 9; i++) {
        elt[i] = std::make_shared<PZsparse>(centers_inp[i]);
    }
}

matPZsparse::matPZsparse(const TYPE* centers_inp, const TYPE uncertainty_percent) {
    for (int i = 0; i < 9; i++) {
        TYPE temp_lb = (1 - uncertainty_percent) * centers_inp[i];
        TYPE temp_ub = (1 + uncertainty_percent) * centers_inp[i];
        if (temp_lb > temp_ub) {
            swap(temp_lb, temp_ub);
        }
        elt[i] = std::make_shared<PZsparse>(Interval(temp_lb, temp_ub));
    }
}

matPZsparse::matPZsparse(const TYPE roll, const TYPE pitch, const TYPE yaw) {
    elt[0] = std::make_shared<PZsparse>(cos(pitch)*cos(yaw));
    elt[1] = std::make_shared<PZsparse>(-cos(pitch)*sin(yaw));
    elt[2] = std::make_shared<PZsparse>(sin(pitch));
    elt[3] = std::make_shared<PZsparse>(cos(roll)*sin(yaw) + cos(yaw)*sin(pitch)*sin(roll));
    elt[4] = std::make_shared<PZsparse>(cos(roll)*cos(yaw) - sin(pitch)*sin(roll)*sin(yaw));
    elt[5] = std::make_shared<PZsparse>(-cos(pitch)*sin(roll));
    elt[6] = std::make_shared<PZsparse>(sin(roll)*sin(yaw) - cos(roll)*cos(yaw)*sin(pitch));
    elt[7] = std::make_shared<PZsparse>(cos(yaw)*sin(roll) + cos(roll)*sin(pitch)*sin(yaw));
    elt[8] = std::make_shared<PZsparse>(cos(pitch)*cos(roll));
}

matPZsparse::matPZsparse(const PZsparse& elt1, const PZsparse& elt2, const PZsparse& elt3) {
    elt[0] = std::make_shared<PZsparse>(elt1);
    elt[1] = std::make_shared<PZsparse>();
    elt[2] = std::make_shared<PZsparse>();
    elt[3] = std::make_shared<PZsparse>();
    elt[4] = std::make_shared<PZsparse>(elt2);
    elt[5] = std::make_shared<PZsparse>();
    elt[6] = std::make_shared<PZsparse>();
    elt[7] = std::make_shared<PZsparse>();
    elt[8] = std::make_shared<PZsparse>(elt3);
}

matPZsparse::matPZsparse(const PZsparse& cosElt, const PZsparse& sinElt, const int axes) {
    PZsparse negSinElt = -1.0 * sinElt;

    switch (axes) {
        case 0: // fixed joints
            // NOTE:
            // This is just a simplified implementation!!!
            // The rotation matrix is just identity matrix for all fixed joints of Fetch
            elt[0] = std::make_shared<PZsparse>(1.0);
            elt[1] = std::make_shared<PZsparse>();
            elt[2] = std::make_shared<PZsparse>();
            elt[3] = std::make_shared<PZsparse>();
            elt[4] = std::make_shared<PZsparse>(1.0);
            elt[5] = std::make_shared<PZsparse>();
            elt[6] = std::make_shared<PZsparse>();
            elt[7] = std::make_shared<PZsparse>();
            elt[8] = std::make_shared<PZsparse>(1.0);
        case 1: // rx
            elt[0] = std::make_shared<PZsparse>(1);
            elt[1] = std::make_shared<PZsparse>();
            elt[2] = std::make_shared<PZsparse>();
            elt[3] = std::make_shared<PZsparse>();
            elt[4] = std::make_shared<PZsparse>(cosElt);
            elt[5] = std::make_shared<PZsparse>(negSinElt);
            elt[6] = std::make_shared<PZsparse>();
            elt[7] = std::make_shared<PZsparse>(sinElt);
            elt[8] = std::make_shared<PZsparse>(cosElt);
            break;
        case 2: // ry
            elt[0] = std::make_shared<PZsparse>(cosElt);
            elt[1] = std::make_shared<PZsparse>();
            elt[2] = std::make_shared<PZsparse>(sinElt);
            elt[3] = std::make_shared<PZsparse>();
            elt[4] = std::make_shared<PZsparse>(1);
            elt[5] = std::make_shared<PZsparse>();
            elt[6] = std::make_shared<PZsparse>(negSinElt);
            elt[7] = std::make_shared<PZsparse>();
            elt[8] = std::make_shared<PZsparse>(cosElt);
            break;
        case 3: // rz
            elt[0] = std::make_shared<PZsparse>(cosElt);
            elt[1] = std::make_shared<PZsparse>(negSinElt);
            elt[2] = std::make_shared<PZsparse>();
            elt[3] = std::make_shared<PZsparse>(sinElt);
            elt[4] = std::make_shared<PZsparse>(cosElt);
            elt[5] = std::make_shared<PZsparse>();
            elt[6] = std::make_shared<PZsparse>();
            elt[7] = std::make_shared<PZsparse>();
            elt[8] = std::make_shared<PZsparse>(1);
            break;
        default:
            throw;
    }
}

void matPZsparse::print() {
    cout << "{\n";

    for (int i = 0; i < 9; i++) {
        if (elt[i] == nullptr) {
            cout << "0\n\n";
        }
        else {
            elt[i]->print();
        }
    }

    cout << "}\n\n";
}

matPZsparse matPZsparse::operator=(const matPZsparse& a) {
    for (int i = 0; i < 9; i++) {
        elt[i] = a.elt[i];
    }

    return *this;
}

vecPZsparse operator*(const matPZsparse m, const vecPZsparse& v) {
    vecPZsparse res;

    res.elt[0] = std::make_shared<PZsparse>( *(m.elt[0]) * *(v.elt[0]) + *(m.elt[1]) * *(v.elt[1]) + *(m.elt[2]) * *(v.elt[2]) );
    res.elt[1] = std::make_shared<PZsparse>( *(m.elt[3]) * *(v.elt[0]) + *(m.elt[4]) * *(v.elt[1]) + *(m.elt[5]) * *(v.elt[2]) );
    res.elt[2] = std::make_shared<PZsparse>( *(m.elt[6]) * *(v.elt[0]) + *(m.elt[7]) * *(v.elt[1]) + *(m.elt[8]) * *(v.elt[2]) );

    return res;
}

matPZsparse operator*(const matPZsparse& m1, const matPZsparse& m2) {
    matPZsparse res;

    res.elt[0] = std::make_shared<PZsparse>( *(m1.elt[0]) * *(m2.elt[0]) + *(m1.elt[1]) * *(m2.elt[3]) + *(m1.elt[2]) * *(m2.elt[6]) );
    res.elt[1] = std::make_shared<PZsparse>( *(m1.elt[0]) * *(m2.elt[1]) + *(m1.elt[1]) * *(m2.elt[4]) + *(m1.elt[2]) * *(m2.elt[7]) );
    res.elt[2] = std::make_shared<PZsparse>( *(m1.elt[0]) * *(m2.elt[2]) + *(m1.elt[1]) * *(m2.elt[5]) + *(m1.elt[2]) * *(m2.elt[8]) );
    res.elt[3] = std::make_shared<PZsparse>( *(m1.elt[3]) * *(m2.elt[0]) + *(m1.elt[4]) * *(m2.elt[3]) + *(m1.elt[5]) * *(m2.elt[6]) );
    res.elt[4] = std::make_shared<PZsparse>( *(m1.elt[3]) * *(m2.elt[1]) + *(m1.elt[4]) * *(m2.elt[4]) + *(m1.elt[5]) * *(m2.elt[7]) );
    res.elt[5] = std::make_shared<PZsparse>( *(m1.elt[3]) * *(m2.elt[2]) + *(m1.elt[4]) * *(m2.elt[5]) + *(m1.elt[5]) * *(m2.elt[8]) );
    res.elt[6] = std::make_shared<PZsparse>( *(m1.elt[6]) * *(m2.elt[0]) + *(m1.elt[7]) * *(m2.elt[3]) + *(m1.elt[8]) * *(m2.elt[6]) );
    res.elt[7] = std::make_shared<PZsparse>( *(m1.elt[6]) * *(m2.elt[1]) + *(m1.elt[7]) * *(m2.elt[4]) + *(m1.elt[8]) * *(m2.elt[7]) );
    res.elt[8] = std::make_shared<PZsparse>( *(m1.elt[6]) * *(m2.elt[2]) + *(m1.elt[7]) * *(m2.elt[5]) + *(m1.elt[8]) * *(m2.elt[8]) );

    return res;
}

void transpose(matPZsparse& a, const matPZsparse& b) {
    a.elt[0] = b.elt[0];
    a.elt[1] = b.elt[3];
    a.elt[2] = b.elt[6];
    a.elt[3] = b.elt[1];
    a.elt[4] = b.elt[4];
    a.elt[5] = b.elt[7];
    a.elt[6] = b.elt[2];
    a.elt[7] = b.elt[5];
    a.elt[8] = b.elt[8];
}

#endif