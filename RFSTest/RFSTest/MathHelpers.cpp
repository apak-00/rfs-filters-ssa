#include "MathHelpers.h"
#include <Eigen/Dense>


double MathHelpers::mahalanobis(VectorXd _v1, VectorXd _v2, MatrixXd _S)
{
	VectorXd d = _v1 - _v2;
	return d.transpose() * _S.llt().solve(d);
}

VectorXd MathHelpers::esf(const VectorXd & _Z)
{
	if (!_Z.size())
		return VectorXd::Ones(1);

	size_t nz = _Z.size();
	MatrixXd F = MatrixXd::Zero(2, nz);
	int in = 0, inm = 1, temp;

	for (size_t n = 0; n < nz; n++) 
	{
		F(in, 0) = F(inm, 0) + _Z(n);
		for (size_t k = 1; k < n; k++)
		{
			if (k == n)
				F(in, k) = _Z(n) * F(inm, k - 1);
			else
				F(in, k) = F(inm, k) + _Z(n) * F(inm, k - 1);
		}
		temp = in;
		in = inm;
		inm = temp;
	}
	
	// TODO: Check // Works (?)
	VectorXd result(nz + 1);
	result << 1, F.block(0, 1, nz, 1);
	
	return result;
}
