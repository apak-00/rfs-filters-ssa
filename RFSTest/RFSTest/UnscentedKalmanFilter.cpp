#include "UnscentedKalmanFilter.h"
#include "unscented_sampler.hpp"
#include <Eigen/Core>

#include <iostream>

using namespace Eigen;
using namespace std;

/*
 * <summary> An empty constructor of the UnscentedKalmanFilter class. </summary>
 */
UnscentedKalmanFilter::UnscentedKalmanFilter() : KalmanFilter() {}

/*
 * <summary> Default constructor of the UnscentedKalmanFilter class. </summary>
 * <param name = "_Q"> State transition noise matrix. </param>
 * <param name = "_dt"> Timestep. Defaults to zero. </param>
 */
UnscentedKalmanFilter::UnscentedKalmanFilter(const decltype(Q) _Q, const decltype(dt) _dt) : 
	KalmanFilter(MatrixXd::Identity(_Q.rows(), _Q.cols()), _Q, _dt) {}

/*
 * <summary> Prediction step of the UnscentedKalmanFilter. </summary>
 * <param name = "_gc"> A Gaussian component to be predicted. </param>
 */ 
void UnscentedKalmanFilter::predict(gaussian_component & _gc)
{
	vector<VectorXd> sigmaPoints, sigmaPointsPredicted;
	vector<double> sigmaWeights;
	icl::standard_unscented_sampler<6, double> sampler;

	size_t L = _gc.m.size() * 2;
	VectorXd recM = VectorXd::Zero(_gc.m.size()), d(_gc.m.size());
	MatrixXd recP = MatrixXd::Zero(_gc.P.rows(), _gc.P.cols());

	double kappa = 0.0, alpha = 1e-3, beta = 0;
	double lambda = alpha * alpha * (L + kappa) - L;

	// Derive the sigma points
	getSigmaPoints(_gc.m, _gc.P, 0.5, sigmaPoints, sigmaWeights);

	sigmaPointsPredicted.resize(sigmaPoints.size());

	// Propagate the points
	for (size_t i = 0; i < sigmaPoints.size(); i++)
	{
		//sigmaPointsPredicted[i] = F * sigmaPoints[i].head(_gc.m.size());
		sigmaPointsPredicted[i] = Astro::integrationPrediction(sigmaPoints[i], dt);
	}
		
	// Reconstruct the mean and covariance
	// Mean
	for (size_t i = 0; i < sigmaPoints.size(); i++)
		recM += sigmaPointsPredicted[i] * sigmaWeights[i];

	// Covariance
	// First term for covariance

	sigmaWeights[0] += 1 - alpha * alpha + beta;
	
	for (size_t i = 0; i < sigmaPoints.size(); i++)
	{
		d = sigmaPointsPredicted[i] - recM;
		recP += d * d.transpose() * sigmaWeights[i];
	}

	//std::cout << "P1: " << _gc.m.transpose() << std::endl;
	//std::cout << "P2: " << recM.transpose() << std::endl;

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
	vector<VectorXd> sigmaPoints, sigmaPointsProjected;
	vector<double> sigmaWeights;
	icl::standard_unscented_sampler<6, double> sampler;

	size_t L = _gc.m.size() + _sensor.getZDim();
	VectorXd recM = VectorXd::Zero(_gc.m.size()), dm(_gc.m.size()),				// Reconstructed mean and temporary difference vector
		recZ = VectorXd::Zero(_sensor.getZDim()), dz(_sensor.getZDim());	// Reconstructed measuremnt and temporary difference vector
		 
	MatrixXd recPZZ = MatrixXd::Zero(_sensor.getZDim(), _sensor.getZDim()),
		recPXZ = MatrixXd::Zero(_gc.m.size(), _sensor.getZDim());

	double kappa = 0.0, alpha = 1e-3, beta = 0;

	MatrixXd R = _sensor.getR();

	// Get sigma points
	getSigmaPoints(_gc.m, _gc.P, 0.5, sigmaPoints, sigmaWeights);

	sigmaPointsProjected.resize(sigmaPoints.size());

	// Project sigma points onto the measurement state
	for (size_t i = 0; i < sigmaPoints.size(); i++)
	{
		sigmaPointsProjected[i] = Astro::temeToRAZEL(sigmaPoints[i], _sensor.getPosition(), _sensor.getDateJD(), 
			_sensor.getLOD(), _sensor.getXp(), _sensor.getYp());
		recZ += sigmaPointsProjected[i].head(3) * sigmaWeights[i];
		
		//cout << sigmaPoints[i].transpose() << endl;
		//cout << sigmaPointsProjected[i].transpose() << endl << endl;
	}

	//sigmaWeights[0] += 1 - alpha * alpha + beta;

	for (size_t i = 0; i < sigmaPoints.size(); i++) 
	{
		dz = sigmaPointsProjected[i].head(3) - recZ;		// Projected - reconstructed measurement
		dm = sigmaPoints[i].head(6) - _gc.m;				// Sigma points - mean

		recPZZ += sigmaWeights[i] * dz * dz.transpose();
		recPXZ += sigmaWeights[i] * dm * dz.transpose();
	}

	// The S matrix
	_sensor.setS(recPZZ);

	// Kalman Update
	// TODO: Use solve
	MatrixXd K = recPXZ * recPZZ.inverse();

	VectorXd updM = _gc.m + K * (_sensor.z[_zNum] - recZ);
	
	//std::cout << "Z : " << recZ.transpose() << std::endl;
	//std::cout << "U1: " << _gc.m.transpose() << std::endl;
	//std::cout << "U2: " << updM.transpose() << std::endl;

	//std::cout << "recPXZ : " << std::endl << recPXZ << endl;
	//std::cout << "recPZZ : " << std::endl << recPZZ << endl;
	//std::cout << "recPZZi : " << std::endl << recPZZ.inverse() << endl;
	//std::cout << "mult : " << std::endl << recPZZ * recPZZ.inverse() << endl;
	
	//cout << "Before: "  << endl <<_gc.m.transpose() << endl;
	//cout << _gc.P << endl;

	_gc.m += K * (_sensor.z[_zNum] - recZ);
	_gc.P -= K * recPZZ * K.transpose();

	//cout << "After: " << endl << _gc.m.transpose() << endl;
	//cout << _gc.P << endl;
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

