#include "MathHelpers.h"
#include <Eigen/Dense>


double MathHelpers::mahalanobis(VectorXd _v1, VectorXd _v2, MatrixXd _S)
{
	VectorXd d = _v1 - _v2;
	return d.transpose() * _S.llt().solve(d);
}