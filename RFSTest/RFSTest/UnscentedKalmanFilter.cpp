#include "UnscentedKalmanFilter.h"
#include "unscented_sampler.hpp"

using namespace Eigen;
using namespace std;

UnscentedKalmanFilter::UnscentedKalmanFilter() : KalmanFilter()
{
}

UnscentedKalmanFilter::UnscentedKalmanFilter(const decltype(Q) _Q, const decltype(dt) _dt) : KalmanFilter(MatrixXd::Identity(_Q.rows(), _Q.cols()), _Q, _dt) 
{
	
}

void UnscentedKalmanFilter::predict(gaussian_component & _gc)
{
	vector<VectorXd> sigmaPoints, sigmaPointsPredicted;
	vector<double> sigmaWeights;
	icl::standard_unscented_sampler<6, double> sampler;

	size_t L = _gc.m.size() * 2;
	VectorXd augM(L), recM = VectorXd::Zero(_gc.m.size()), d(_gc.m.size());
	MatrixXd augP(L, L), recP = VectorXd::Zero(_gc.P.rows(), _gc.P.cols());

	double kappa = 0.0, alpha = 1e-3, beta = 0;
	double lambda = alpha * alpha * (L + kappa) - L;

	// Augment the state
	augM << _gc.m, getCVq(_gc.m.size(), dt);
	augP << _gc.P, MatrixXd::Zero(_gc.P.rows(), _gc.P.cols()), MatrixXd::Zero(_gc.P.rows(), _gc.P.cols()), Q;

	// Derive the sigma points
	sampler.get_points(augM, augP, 0.5, sigmaPoints, sigmaWeights);

	sigmaPointsPredicted.resize(sigmaPoints.size());

	// Propagate the points
	for (size_t i = 0; i < sigmaPoints.size(); i++)
		sigmaPointsPredicted[i] = Astro::integrationPrediction(sigmaPoints[i].head(_gc.m.size()), dt);
	
	// Reconstruct the mean and covariance
	// Mean
	for (size_t i = 0; i < 2 * L; i++)
		recM += sigmaPointsPredicted[i] * sigmaWeights[i];

	// Covariance
	// First term for covariance

	sigmaWeights[0] += 1 - alpha * alpha + beta;
	
	for (size_t i = 0; sigmaPoints.size(); i++)
	{
		d = sigmaPointsPredicted[i] - recM;
		recP += d * d.transpose() * sigmaWeights[i];
	}

	// Reassign the values
	_gc.m = recM;
	_gc.P = recP;
}

void UnscentedKalmanFilter::update(gaussian_component & _gc, Sensor& _sensor, const size_t& _zNum)
{
	vector<VectorXd> sigmaPoints, sigmaPointsProjected;
	vector<double> sigmaWeights;
	icl::standard_unscented_sampler<6, double> sampler;

	size_t L = _gc.m.size() * + _sensor.getZDim();
	VectorXd augM(L), r(3),
		recM = VectorXd::Zero(_gc.m.size()), dm(_gc.m.size()),				// Reconstructed mean and temporary difference vector
		recZ = VectorXd::Zero(_sensor.getZDim()), dz(_sensor.getZDim());	// Reconstructed measuremnt and temporary difference vector
		 
	MatrixXd augP(L, L), 
		recPZZ = VectorXd::Zero(_sensor.getZDim(), _sensor.getZDim()),
		recPXZ = VectorXd::Zero(_gc.m.size(), _sensor.getZDim());

	double kappa = 0.0, alpha = 1e-3, beta = 0;

	r << 0.75, 0.05, 0.05; r *= 5;		// Measurement noise vector

	// Augment the state
	augM << _gc.m, r;
	augP << _gc.P, 
		MatrixXd::Zero(_gc.P.rows(), _sensor.getR().cols()), 
		MatrixXd::Zero(_sensor.getR().rows(), _gc.P.cols()), 
		_sensor.getR();

	// Get sigma points
	sampler.get_points(_gc.m, _sensor.getR(), 0.5, sigmaPoints, sigmaWeights);

	sigmaPointsProjected.resize(sigmaPoints.size());

	// Project sigma points onto the measurement state
	for (size_t i = 0; i < sigmaPoints.size(); i++)
	{
		sigmaPointsProjected[i] = Astro::temeToRAZEL(sigmaPoints[i].head(6), _sensor.getPosition(), _sensor.getDateJD(), _sensor.getLOD(), _sensor.getXp(), _sensor.getYp());
		recZ += sigmaPointsProjected[i].head(3) * sigmaWeights[i];
	}

	sigmaWeights[0] += 1 - alpha * alpha + beta;

	for (size_t i = 0; i < sigmaPoints.size(); i++) 
	{
		dz = sigmaPointsProjected[i].head(3) - recZ;		// Projected - reconstructed measurement
		dm = sigmaPoints[i].head(6) - _gc.m;				// Sigma points - mean

		recPZZ += sigmaWeights[i] * dz * dz.transpose();
		recPXZ += sigmaWeights[i] * dm * dz.transpose();
	}

	// Kalman Update
	MatrixXd K = recPXZ * recPZZ.inverse();

	_gc.m += K * (_sensor.z[_zNum] - recZ);
	_gc.P -= K * recPZZ * K.transpose();
}

