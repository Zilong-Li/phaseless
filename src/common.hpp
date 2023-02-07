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

/*
** @param gli  genotype likelihoods of current individual i,3 x nsnps
*/
inline ArrDouble2D emissionCurIterInd(const ArrDouble2D & gli, const ArrDouble2D & F, bool use_log)
{
    int k1, k2, g1, g2;
    int M = F.rows();
    int C = F.cols();
    ArrDouble2D emitDip(M, C * C); // emission probabilies, nsnps x (C x C)
    for(k1 = 0; k1 < C; k1++)
    {
        for(k2 = 0; k2 < C; k2++)
        {
            emitDip.col(k1 * C + k2).setZero();
            for(g1 = 0; g1 <= 1; g1++)
            {
                for(g2 = 0; g2 <= 1; g2++)
                {
                    emitDip.col(k1 * C + k2) += gli.row(g1 + g2).transpose()
                                                * (g1 * F.col(k1) + (1 - g1) * (1 - F.col(k1)))
                                                * (g2 * F.col(k2) + (1 - g2) * (1 - F.col(k2)));
                }
            }
        }
    }
    if(use_log)
        return emitDip.log();
    else
    {
        // be careful with underflow
        const double maxEmissionMatrixDifference = 1e-10;
        auto x = emitDip.rowwise().maxCoeff();
        emitDip = emitDip.colwise() / x;
        emitDip = (emitDip < maxEmissionMatrixDifference).select(maxEmissionMatrixDifference, emitDip);

        return emitDip;
    }
}

#endif // COMMON_H_
