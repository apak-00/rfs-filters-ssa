#pragma once
#include <Eigen/Core>

using namespace Eigen;

namespace MathHelpers {

	double mahalanobis(VectorXd _v1, VectorXd _v2, MatrixXd _S);

	VectorXd esf(const VectorXd& _z);
	
}