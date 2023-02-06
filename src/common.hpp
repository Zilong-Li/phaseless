#ifndef COMMON_H_
#define COMMON_H_

#include <Eigen/Dense>
#include <random>
#include <vector>

using MyMat2D = Eigen::MatrixXf; // use MatrixXf if no accuracy drop
using MyArr2D = Eigen::ArrayXXf; // use ArrayXf if no accuracy drop
using MatDouble2D = Eigen::MatrixXd; // use matrix for linear algebra operations
using ArrDouble2D = Eigen::ArrayXXd; // use array for element-wise operations
using ArrDouble1D = Eigen::ArrayXd;
using ArrFloat2D = Eigen::ArrayXXf;
using ArrFloat1D = Eigen::ArrayXf;
using DoubleVec1D = std::vector<double>;
using FloatVec1D = std::vector<float>;

template<typename MatrixType, typename RandomEngineType>
inline MatrixType RandomUniform(const Eigen::Index numRows,
                                const Eigen::Index numCols,
                                RandomEngineType & engine,
                                double a,
                                double b)
{
    std::uniform_real_distribution<typename MatrixType::Scalar> uniform_real_distribution{a, b}; // or using 0.05, 0.95
    const auto uniform{[&](typename MatrixType::Scalar) { return uniform_real_distribution(engine); }};
    return MatrixType::NullaryExpr(numRows, numCols, uniform);
};

#endif // COMMON_H_
