#include "UnscentedKalmanFilter.h"
#include "unscented_sampler.hpp"
#include <Eigen/Core>

#include <iostream>

using namespace Eigen;
using namespace std;

/*
 * <summary> An empty constructor of the UnscentedKalmanFilter class. </summary>
 */
UnscentedKalmanFilter::UnscentedKalmanFilter() : KalmanFilter(), sigmaSamplingW(0.5) {}

/*
 * <summary> Default constructor of the UnscentedKalmanFilter class. </summary>
 * <param name = "_Q"> State transition noise matrix. </param>
 * <param name = "_dt"> Timestep. Defaults to zero. </param>
 */
UnscentedKalmanFilter::UnscentedKalmanFilter(const decltype(Q) _Q, const decltype(dt) _dt) : 
	KalmanFilter(MatrixXd::Identity(_Q.rows(), _Q.cols()), _Q, _dt), sigmaSamplingW(0.5) {}

UnscentedKalmanFilter::UnscentedKalmanFilter(const decltype(Q) _Q, const decltype(sigmaSamplingW) _sigmaSamplingW, const decltype(dt) _dt) :
	KalmanFilter(MatrixXd::Identity(_Q.rows(), _Q.cols()), _Q, _dt), sigmaSamplingW(_sigmaSamplingW) {}

/*
 * <summary> Prediction step of the UnscentedKalmanFilter. </summary>
 * <param name = "_gc"> A Gaussian component to be predicted. </param>
 */ 
void UnscentedKalmanFilter::predict(gaussian_component & _gc)
{
	std::vector<VectorXd> sigmaPoints, sigmaPointsPredicted;
	std::vector<double> sigmaWeights;
	icl::standard_unscented_sampler<6, double> sampler;
	int md = _gc.m.size(), zd = 3;

	size_t L = _gc.m.size() * 2;
	VectorXd recM = VectorXd::Zero(_gc.m.size()), d(_gc.m.size());
	MatrixXd recP = MatrixXd::Zero(_gc.P.rows(), _gc.P.cols());

	double kappa = 0.0, alpha = 1e-3, beta = 2;
	double lambda = alpha * alpha * (L + kappa) - L;

	// Noise
	size_t nSize = 3;
	MatrixXd noiseP = MatrixXd::Identity(nSize, nSize);
	noiseP *= 0.001;

	// Augmentation
	VectorXd augM = VectorXd::Zero(md + nSize);
	MatrixXd augP = MatrixXd::Zero(augM.size(), augM.size());
	augM.head(md) = _gc.m;
	augP.block(0, 0, md, md) = _gc.P;
	augP.block(md, md, nSize, nSize) = noiseP;

	// Derive the sigma points
	getSigmaPoints(augM, augP, 0.5, sigmaPoints, sigmaWeights);

	sigmaPointsPredicted.resize(sigmaPoints.size());

	// Propagate the points
	for (size_t i = 0; i < sigmaPoints.size(); i++)
	{
		auto mean = sigmaPoints[i].head(md);
		mean.tail(nSize) += sigmaPoints[i].tail(nSize);
		sigmaPointsPredicted[i] = Astro::integrationPrediction(mean, dt);		// First 6 values of the sigma point (mean)
	}

	// Reconstruct the mean and covariance
	// Mean
	for (size_t i = 0; i < sigmaPoints.size(); i++)
		recM += sigmaPointsPredicted[i] * sigmaWeights[i];

	// Covariance
	// First term for covariance
	//sigmaWeights[0] += 1 - alpha * alpha + beta;

	for (size_t i = 0; i < sigmaPoints.size(); i++)
	{
		d = sigmaPointsPredicted[i] - recM;
		recP += d * d.transpose() * sigmaWeights[i];
	}

	//std::cout << "P1: " << std::endl << _gc.P << std::endl;
	//std::cout << "P2: " << std::endl << recP << std::endl;

	// Reassign the values
	_gc.m = recM;
	_gc.P = recP;
}

/*
 * <summary> Update step of the UnscentedKalmanFilter. </summary>
 * <param name = "_gc"> A Gaussian component to be updated. </param>
 * <param name = "_sensor"> Reference to the sensor class containing the measurements. </param>
 * <param name = "_zNum"> Mesruement number (for multiple-measurement case). </param>
 */
void UnscentedKalmanFilter::update(gaussian_component & _gc, Sensor& _sensor, const size_t& _zNum)
{
	debug = false;
	vector<VectorXd> sigmaPoints, sigmaPointsProjected;
	vector<double> sigmaWeights;
	icl::standard_unscented_sampler<6, double> sampler;
	int md = _gc.m.size(), zd = _sensor.getZDim();

	size_t L = _gc.m.size() + 1;
	VectorXd recM = VectorXd::Zero(_gc.m.size()), dm(_gc.m.size()),			// Reconstructed mean and temporary difference vector
		recZ = VectorXd::Zero(_sensor.getZDim()), dz(_sensor.getZDim());	// Reconstructed measuremnt and temporary difference vector
	MatrixXd recPZZ = MatrixXd::Zero(_sensor.getZDim(), _sensor.getZDim()),	// S
		recPXZ = MatrixXd::Zero(_gc.m.size(), _sensor.getZDim());

	double kappa = 0.0, alpha = 1e-3, beta = 2;

	MatrixXd R = _sensor.getR(), K;

	VectorXd n(1);
	n << 0.08;
	MatrixXd noiseP = MatrixXd::Identity(n.size(), n.size());
	noiseP.diagonal() = n;

	// Augmentation
	VectorXd augM = VectorXd::Zero(md + n.size());
	MatrixXd augP = MatrixXd::Zero(augM.size(), augM.size());

	augM.head(md) = _gc.m;
	augP.block(0, 0, md, md) = _gc.P;
	augP.block(md, md, n.size(), n.size()) = noiseP;

	// Get sigma points
	getSigmaPoints(augM, augP, sigmaSamplingW, sigmaPoints, sigmaWeights);

	sigmaPointsProjected.resize(sigmaPoints.size());

	// Project sigma points onto the measurement state
	for (size_t i = 0; i < sigmaPoints.size(); i++)
	{
		auto mean = sigmaPoints[i].head(md);
		sigmaPointsProjected[i] = Astro::temeToRAZEL(mean, _sensor.getPosition(), _sensor.getDateJD(), 
			_sensor.getLOD(), _sensor.getXp(), _sensor.getYp());
		sigmaPointsProjected[i].head(n.size()) += sigmaPoints[i].tail(n.size());
		recZ += sigmaPointsProjected[i].head(3) * sigmaWeights[i];
	}

	if (debug)
	{
		cout << _gc.P << endl;
		cout << "Sigma Points: " << endl;
		for (size_t i = 0; i < sigmaPoints.size(); i++)
			cout << sigmaPoints[i].transpose() << endl;
	
		cout << "Sigma Points Projected: " << endl;
		for (size_t i = 0; i < sigmaPointsProjected.size(); i++)
			cout << sigmaPointsProjected[i].transpose() << endl;
		
	}

	//sigmaWeights[0] += 1 - alpha * alpha + beta;

	for (size_t i = 0; i < sigmaPoints.size(); i++) 
	{
		dz = sigmaPointsProjected[i].head(zd) - recZ;		// Projected - reconstructed measurement
		dm = sigmaPoints[i].head(md) - _gc.m;				// Sigma points - mean

		recPZZ += sigmaWeights[i] * dz * dz.transpose();
		recPXZ += sigmaWeights[i] * dm * dz.transpose();
	}

	// The S matrix
	_sensor.setS(recPZZ);

	// Kalman Update
	// TODO: Use solve
	MatrixXd KInv, KSol, tempKDiff;
	KInv = recPXZ * recPZZ.inverse();
	KSol = recPZZ.transpose().llt().solve(recPXZ.transpose()).transpose();
	tempKDiff = (KInv - KSol).cwiseAbs();
	
	if (debug)
	{
		cout << "K: (inverse): " << endl << KInv << endl;
		cout << "K: (solve): " << endl << KSol << endl;
		cout << "Diff: " << endl << tempKDiff << endl;
	}
	
	K = KSol;

	if (debug)
	{
		cout << "recZ: " << recZ.transpose() << endl;
		cout << "recPZZ: " << endl << recPZZ << endl;
		cout << "recPXZ: " << endl << recPXZ << endl;
		cout << "K: " << endl << K << endl;
		cout << "recPZZ (inverse): " << endl << recPZZ.inverse() << endl;
	}

	VectorXd oldM = _gc.m;
	MatrixXd oldP = _gc.P;

	_gc.m += K * (_sensor.z[_zNum] - recZ);
	_gc.P -= K * recPZZ * K.transpose();

	//std::cout << "P_U_1: " << std::endl << oldP << endl;
	//std::cout << "P_U_2: " << std::endl << _gc.P << endl;
}

/* 
* <summary> Get the sigma points and the corresponding weights for the specified state vector and its covariance matrix. </summary>
* <param name = "_mean"> State vector. </mean>
* <param name = "_cov"> State vector covariance </mean>
* <param name = "_w0"> Zero weight. </mean>
* <param name = "_points"> Output sigma points. </param>
* <param name = "_weights"> Output sigma points' weights. </param>
*/
void UnscentedKalmanFilter::getSigmaPoints(const VectorXd &_mean, const MatrixXd &_cov, const double &_w0, 
	 vector<Eigen::VectorXd>& _points, vector<double>& _weights)
{
	size_t dim_point = _mean.rows();
	size_t num_sigma = 2 * dim_point + 1;
	_points.resize(num_sigma);
	_weights.resize(num_sigma);

	// Fill the sigma weights
	double w1 = (1.0 - _w0) / (2.0 * (double)dim_point);
	_weights[0] = _w0;
	_points[0] = _mean;
	fill(_weights.begin() + 1, _weights.end(), w1);
	MatrixXd sqS = (dim_point / (1.0 - _w0) * _cov).llt().matrixL();
	for (size_t i = 0; i < dim_point; i++)
	{
		_points[1 + i] = _mean + sqS.col(i);
	}
	for (size_t i = 0; i < dim_point; i++)
	{
		_points[1 + i + dim_point] = _mean - sqS.col(i);
	}
}

