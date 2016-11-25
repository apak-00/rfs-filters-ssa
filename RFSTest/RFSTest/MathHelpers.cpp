#include "MathHelpers.h"
#include <Eigen/Dense>

/*
 * <summary> Computer Mahalanobis distance for a given pair of vectors and a common covariance. </summary>
 */
double MathHelpers::mahalanobis(VectorXd _v1, VectorXd _v2, MatrixXd _S)
{
	VectorXd d = _v1 - _v2;
	return d.transpose() * _S.llt().solve(d);
}

/*
 * <summary> Return the value of the univariate Laplace distribution for a give set of values. </summary>
 * <param name = "_x"> Vaiable. </param>
 * <param name = "_mu"> Location parameter (e.g. mean). </param>
 * <param name = "_b"> Scale parameter. </param>
 */
double MathHelpers::laplace(const double & _x, const double & _mu, const double & _b)
{
	return 0.5 / _b * exp(-abs((_x - _mu)) / _b);
}

/*
 * <summary> Compute Gaussian likelihood for a given pair of vectors and a common covariance. </summary>
 */
double MathHelpers::gaussianLikelihood(const VectorXd _z, const VectorXd _zPred, const MatrixXd _S)
{
	return (1.0 / sqrt(pow(2.0 * M_PI, _z.size()) * _S.determinant())) * exp(-0.5 * mahalanobis(_z, _zPred, _S));
}

 /**
 * <summary> [Old, Temporary] Returns a VectorXd containing random values within the specified range. </summary>
 * <param name = "_lowerBound"> </param>
 * <param name = "_lowerBound"> </param>
 * <returns> A VectorXd with random values. </returns>
 */
 VectorXd MathHelpers::randVec(const VectorXd & _lowerBound, const VectorXd & upperBound)
 {
	 size_t l = _lowerBound.size();
	 VectorXd result = VectorXd::Zero(l);

	 for (size_t i = 0; i < l; i++)
		 result(i) = (double)(upperBound(i) - _lowerBound(i)) * (double)rand() / (double)RAND_MAX - (upperBound(i) - _lowerBound(i)) / 2;

	 return result;
 }

 /*
 * <summary> Computes elementary symmetric functions. </summary>
 * <par> For CPHD filter. </par>
 * TODO: Complete
 */
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
