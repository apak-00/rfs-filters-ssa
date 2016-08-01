#include "ExtendedKalmanFilter.h"
#include "unscented_sampler.hpp"
#include <iostream>

/**
 * <summary> Empty constructor of the ExtendedKalmanFilter class. </summary>
 */
ExtendedKalmanFilter::ExtendedKalmanFilter() : KalmanFilter() {}

/**
 * <summary> Standard constructor of the ExtendedKalmanFilter class. </summary>
 * <par> Temporary: F is present. </par>
 * <param name = "_F"> Transition matrix. </param>
 * <param name = "_Q"> Process noise matrix. </param>
 * <param name = "_dt"> Timestep. </param>
 */
ExtendedKalmanFilter::ExtendedKalmanFilter(const decltype(F)& _F, const decltype(Q) _Q, const decltype(dt) _dt) : KalmanFilter(_F, _Q, _dt) {}

/**
 * <summary> Copy constructor of the ExtendedKalmanFilter class. </summary>
 * <param name = "_ExtendedKalmanFilter"> An instance of the ExtendedKalmanFilter to copy from. </param> 
 */
ExtendedKalmanFilter::ExtendedKalmanFilter(const ExtendedKalmanFilter & _ExtendedKalmanFilter) : KalmanFilter(_ExtendedKalmanFilter) {}

/**
 * <summary> ExtendedKalmanFilter prediciton. </summary>
 * <param name = "_gc"> A gaussian_component to be predicted. </param>
 */
void ExtendedKalmanFilter::predict(gaussian_component & _gc)
{
	if (false)
	{
		if (_gc.P.block<3, 3>(0, 0).determinant() < 56)
			_gc.kindaConverged = true;

		if (false && _gc.kindaConverged)
		{
			VectorXd temp(_gc.m.size());
			MatrixXd shepperd = Astro::getShepperdMatrix(_gc.m, dt, temp, Astro::MU_E);
			_gc.m = temp;
			//_gc.P = shepperd * _gc.P * shepperd.transpose() + Q;
		}
		else
		{
			//_gc.m = F * _gc.m;		// Constant Velocity Mean Update
			_gc.m = Astro::integrationPrediction(_gc.m, dt);
		}

		_gc.P = F * _gc.P * F.transpose() + Q;
	}
	else
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
			noiseP *= 1e-6;

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
	}
}

/**
* <summary> ExtendedKalmanFilter update. </summary>
* <param name = "_gc"> A gaussian_component to be updated. </param>
* <param name = "_sensor"> A reference to the sensor containing the measurements. </param>
* <param name = "_zNum"> Measurement number. </param>
*/
void ExtendedKalmanFilter::update(gaussian_component & _gc, Sensor & _sensor, const size_t & _zNum)
{
	MatrixXd tf = Astro::getSEZToTEMECovTfMat(_sensor.getPosition(), _sensor.getDateJD(), _sensor.getXp(), _sensor.getYp(), _gc.m.size());
	MatrixXd pSEZ = tf.transpose() * _gc.P * tf;
	VectorXd mSEZ = Astro::temeToSEZ(_gc.m, _sensor.getPosition(), _sensor.getDateJD(), _sensor.getLOD(), _sensor.getXp(), _sensor.getYp());

	MatrixXd H = Astro::getSEZToRAZELJacobian(mSEZ, _sensor.getZDim());			// Hardcoded Jacobian
	//MatrixXd H = Astro::getSEZToRAZELJacobianFADBAD(mSEZ, _sensor.getZDim());
	
	MatrixXd S = H * pSEZ * H.transpose() + _sensor.R;
	MatrixXd K = pSEZ * H.transpose() * S.inverse();

	VectorXd predZ = Astro::sezToRAZEL(mSEZ).head(3);

	mSEZ += K * (_sensor.z[_zNum] - predZ);
	pSEZ = (MatrixXd::Identity(_sensor.sDim, _sensor.sDim) - K * H) * pSEZ;

	VectorXd oldM = _gc.m;
	MatrixXd oldP = _gc.P;

	_gc.m = Astro::sezToTEME(mSEZ, _sensor.getPosition(), _sensor.getDateJD(), _sensor.getLOD(), _sensor.getXp(), _sensor.getYp());
	_gc.P = tf * pSEZ * tf.transpose();

	_sensor.setH(H);
	_sensor.setS(S);
	
	// Temporay fix for the predicted measurement
	// TODO: Change in the future
	_sensor.setPredictedZ(predZ);
}

/*
* <summary> Get the sigma points and the corresponding weights for the specified state vector and its covariance matrix. </summary>
* <param name = "_mean"> State vector. </mean>
* <param name = "_cov"> State vector covariance </mean>
* <param name = "_w0"> Zero weight. </mean>
* <param name = "_points"> Output sigma points. </param>
* <param name = "_weights"> Output sigma points' weights. </param>
*/
void ExtendedKalmanFilter::getSigmaPoints(const VectorXd &_mean, const MatrixXd &_cov, const double &_w0,
	std::vector<Eigen::VectorXd>& _points, std::vector<double>& _weights)
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