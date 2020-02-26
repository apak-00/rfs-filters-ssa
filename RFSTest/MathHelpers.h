#pragma once
#include <Eigen/Core>

using namespace Eigen;

namespace MathHelpers {

	double mahalanobis(VectorXd _v1, VectorXd _v2, MatrixXd _S);
	double laplace(const double& _x, const double& _mu, const double& _b);
	double gaussianLikelihood(const VectorXd _z, const VectorXd _zPred, const MatrixXd _S);

	template<typename T>
	double hellinger(const T& _c1, const T & _c2);

	VectorXd randVec(const VectorXd & _lowerBound, const VectorXd & upperBound);

	VectorXd esf(const VectorXd& _z);
}

/**
* <summary> Hellinger distance. </summary>
* <param name = "_v1"> First mean. </param>
* <param name = "_v1"> First mean's covariance. </param>
* <param name = "_v1"> Second mean. </param>
* <param name = "_v1"> Second mean's covariance. </param>
* <returns> Hellinger distance. </returns>
*/
template<typename T>
inline double MathHelpers::hellinger(const T& _c1, const T & _c2)
{
	VectorXd vDiff = _c1.m - _c2.m;
	MatrixXd pSum = _c1.P + _c2.P;
	double epsilon = (-0.25 * vDiff.transpose() * pSum.inverse() * vDiff)(0, 0);
	return 1 - sqrt(sqrt((_c1.P * _c2.P).determinant()) / (0.5 * pSum).determinant()) * exp(epsilon);
}