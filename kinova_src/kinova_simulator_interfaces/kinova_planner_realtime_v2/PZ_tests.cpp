#include "Trajectory.h"

int main() {
    Eigen::Matrix<double,1,3> center1;
    center1 << 1,2,3;
    PZsparse a(center1);

    Eigen::Matrix<double,3,1> center2;
    center2 << 4,5,6;
    PZsparse b(center2);

    cout << a << endl;

    cout << b << endl;

    PZsparse c = a * b;

    // double factor[NUM_FACTORS] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7};
    // MatrixXInt res = a.slice(factor);

    // for (int i = 0; i < res.rows(); i++) {
    //     for (int j = 0; j < res.cols(); j++) {
    //         cout << "[ " << res(i,j).lower() << ", " << res(i,j).upper() << " ], ";
    //     }
    //     cout << endl;
    // }

    cout << c << endl;
    
    return 0;
}