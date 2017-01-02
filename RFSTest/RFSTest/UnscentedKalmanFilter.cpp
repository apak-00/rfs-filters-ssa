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
	size_t stateSize = _gc.m.size(), noiseSize = 3;

	VectorXd recM = VectorXd::Zero(stateSize), d(stateSize);
	MatrixXd recP = MatrixXd::Zero(stateSize, stateSize);

	// Noise (acceleration)
	MatrixXd noiseP = MatrixXd::Identity(noiseSize, noiseSize);
	if (!noiseSize)
		noiseP *= 1e-2;

	// Augmentation
	VectorXd augM = VectorXd::Zero(stateSize + noiseSize);
	MatrixXd augP = MatrixXd::Zero(augM.size(), augM.size());
	augM.head(stateSize) = _gc.m;
	augP.block(0, 0, stateSize, stateSize) = _gc.P;

	if (!noiseSize)
		augP.block(stateSize, stateSize, noiseSize, noiseSize) = noiseP;

	// Derive the sigma points
	getSigmaPoints(augM, augP, 0.5, sigmaPoints, sigmaWeights);
	sigmaPointsPredicted.resize(sigmaPoints.size());
	auto mean = sigmaPoints[0].head(stateSize);

	// Propagate the points
	for (size_t i = 0; i < sigmaPoints.size(); i++)
	{
		// First 6 values of the sigma point (mean)
		mean = sigmaPoints[i].head(stateSize);
		if (!noiseSize)
			sigmaPointsPredicted[i] = Astro::integrationPrediction(mean, dt);
		else
			// Last 3 values of the sigma point -> acceleration noise
			sigmaPointsPredicted[i] = Astro::integrationPrediction(mean, dt, sigmaPoints[i].tail(3));
	}

	//sigmaWeights[0] += 1 - alpha * alpha + beta;

	// Reconstruct the mean and covariance
	// Mean
	for (size_t i = 0; i < sigmaPoints.size(); i++)
		recM += sigmaPointsPredicted[i] * sigmaWeights[i];

	// Covariance
	for (size_t i = 0; i < sigmaPoints.size(); i++)
	{
		d = sigmaPointsPredicted[i] - recM;
		recP += d * d.transpose() * sigmaWeights[i];
	}

	// Reassign the values
	_gc.m = recM;
	_gc.P = recP;

	MatrixXd tempVarianceMatrix = MatrixXd::Identity(stateSize, stateSize) * 0.1;
	//tempVarianceMatrix.diagonal() = tempVarianceVector;

	_gc.P += tempVarianceMatrix;
}

/*
 * <summary> Update step of the UnscentedKalmanFilter. </summary>
 * <param name = "_gc"> A Gaussian component to be updated. </param>
 * <param name = "_sensor"> Reference to the sensor class containing the measurements. </param>
 * <param name = "_zNum"> Mesruement number (for multiple-measurement case). </param>
 */
void UnscentedKalmanFilter::update(gaussian_component & _gc, Sensor& _sensor, const size_t& _zNum)
{
	vector<VectorXd> sigmaPoints, sigmaPointsProjected;
	vector<double> sigmaWeights;
	icl::standard_unscented_sampler<6, double> sampler;
	size_t stateSize = _gc.m.size(), noiseSize = 3, zSize = _sensor.getZ(0).size();

	VectorXd recM = VectorXd::Zero(_gc.m.size()), dm(_gc.m.size()),			// Reconstructed mean and temporary difference vector
		recZ = VectorXd::Zero(_sensor.getZDim()), dz(_sensor.getZDim());	// Reconstructed measuremnt and temporary difference vector
	MatrixXd recPZZ = MatrixXd::Zero(_sensor.getZDim(), _sensor.getZDim()),	// S
		recPXZ = MatrixXd::Zero(_gc.m.size(), _sensor.getZDim());

	MatrixXd R = _sensor.getR(), K;
	VectorXd noise(3);
	noise << 0.05, 1e-4, 1e-4;
	MatrixXd noiseP = MatrixXd::Identity(noiseSize, noiseSize);
	noiseP.diagonal() = noise;

	// Augmentation
	VectorXd augM = VectorXd::Zero(stateSize + noiseSize);
	MatrixXd augP = MatrixXd::Zero(augM.size(), augM.size());

	augM.head(stateSize) = _gc.m;
	augP.block(0, 0, stateSize, stateSize) = _gc.P;
	augP.block(stateSize, stateSize, noiseSize, noiseSize) = noiseP;

	// Get sigma points
	getSigmaPoints(augM, augP, sigmaSamplingW, sigmaPoints, sigmaWeights);

	sigmaPointsProjected.resize(sigmaPoints.size());

	// Project sigma points onto the measurement state
	auto mean = sigmaPoints[0].head(stateSize);
	for (size_t i = 0; i < sigmaPoints.size(); i++)
	{
		mean = sigmaPoints[i].head(stateSize);
		sigmaPointsProjected[i] = Astro::temeToRAZEL(mean, _sensor.getPosition(), _sensor.getDateJD(), 
			_sensor.getLOD(), _sensor.getXp(), _sensor.getYp());
		sigmaPointsProjected[i].head(noiseSize) += sigmaPoints[i].tail(noiseSize);
		recZ += sigmaPointsProjected[i].head(3) * sigmaWeights[i];
	}

	//sigmaWeights[0] += 1 - alpha * alpha + beta;

	for (size_t i = 0; i < sigmaPoints.size(); i++) 
	{
		dz = sigmaPointsProjected[i].head(zSize) - recZ;		// Projected - reconstructed measurement
		dm = sigmaPoints[i].head(stateSize) - _gc.m;				// Sigma points - mean

		recPZZ += sigmaWeights[i] * dz * dz.transpose();
		recPXZ += sigmaWeights[i] * dm * dz.transpose();
	}

	// The S matrix
	_sensor.setS(recPZZ);
	// Temporary fix for the update measurement
	// TODO:Change in the future
	_sensor.setPredictedZ(recZ);

	// Kalman Update
	K = recPZZ.transpose().llt().solve(recPXZ.transpose()).transpose();

	VectorXd oldM = _gc.m;
	MatrixXd oldP = _gc.P;

	_gc.m += K * (_sensor.z[_zNum] - recZ);
	_gc.P -= K * recPZZ * K.transpose();
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

