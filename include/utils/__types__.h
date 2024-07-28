#ifndef __TYPES__H
#define __TYPES__H

#include <Eigen/Dense>
#include <complex>

// uint
typedef unsigned int uint;

// linalg vector operations
typedef Eigen::VectorXd vector;           // 1D dynamic vector type
typedef Eigen::VectorXcd complex_vector;  //     ~(complex)~

// linear algebra matrix operations
typedef Eigen::MatrixXd matrix;           // 2D dynamic matrix type
typedef Eigen::MatrixXcd complex_matrix;  //     ~(complex)~

// component wise operations
typedef Eigen::ArrayXXd array;           // 2D dynamic array type
typedef Eigen::ArrayXXcd complex_array;  //     ~(complex)~


// fixed array more performant when sizes are known at compile time
template <typename type, int n, int m>
using ndarray = Eigen::Array<type, n, m>;

#endif  // __TYPES__H