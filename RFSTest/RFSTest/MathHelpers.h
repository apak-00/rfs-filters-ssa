#pragma once
#include <Eigen/Core>

using namespace Eigen;

namespace MathHelpers {

	double mahalanobis(VectorXd _v1, VectorXd _v2, MatrixXd _S);
	double laplace(const double& _x, const double& _mu, const double& _b);
	double gaussianLikelihood(const VectorXd _z, const VectorXd _zPred, const MatrixXd _S);

	VectorXd esf(const VectorXd& _z);
	
}