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

	mSEZ += K * (_sensor.z[_zNum] - Astro::sezToRAZEL(mSEZ).head(3));
	pSEZ = (MatrixXd::Identity(_sensor.sDim, _sensor.sDim) - K * H) * pSEZ;

	VectorXd oldM = _gc.m;
	MatrixXd oldP = _gc.P;

	_gc.m = Astro::sezToTEME(mSEZ, _sensor.getPosition(), _sensor.getDateJD(), _sensor.getLOD(), _sensor.getXp(), _sensor.getYp());
	_gc.P = tf * pSEZ * tf.transpose();

	//std::cout << "P_U_1: " << std::endl << oldP << std::endl;
	//std::cout << "P_U_2: " << std::endl << _gc.P << std::endl;

	_sensor.setH(H);
	_sensor.setS(S);
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