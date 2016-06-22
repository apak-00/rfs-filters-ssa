#include "KalmanFilter.h"

// Temp
#include <iostream>

KalmanFilter::KalmanFilter(): SingleTargetFilter() {}

/**
 * <summary> Basic constructor of the KalmanFilter class. </summary>
 * <param name = "_F"> Transition matrix. </param>
 * <param name = "_Q"> Transition noise matrix. </param>
 * <param name = "_t"> Timestep </param>
 */
KalmanFilter::KalmanFilter(const decltype(F)& _F, const decltype(Q) _Q, const decltype(dt) _dt) 
	: SingleTargetFilter(_dt) ,F(_F), Q(_Q)
{
	assert(F.rows() == F.cols() && "F Matrix is not square.");
	assert(Q.rows() == Q.cols() && "Q Matrix is not square.");
}

/**
 * <summary> Copy constructor of the KalmanFilter class. </summary>
 * <param name = "_kf"> An instance of Kalman Filter to copy from. </param>
 */
KalmanFilter::KalmanFilter(const KalmanFilter & _kf)
	: F(_kf.F), Q(_kf.Q), SingleTargetFilter(_kf.dt) {}

/**
 * <summary> Assignment operator overloading. </summary>
 * <param name = "_kf"> An instance of Kalman Filter to copy from. </param>
 */
KalmanFilter& KalmanFilter::operator=(const KalmanFilter & _kf)
{
	F = _kf.F;
	Q = _kf.Q;
	dt = _kf.dt;

	return *this;
}

/**
 * <summary> Kalman Filter prediciton. </summary>
 * <param name = "_gc"> A Gaussian component to predict. </param>
 */
void KalmanFilter::predict(gaussian_component & _gc)
{
	_gc.m = F * _gc.m;
	_gc.P = F * _gc.P * F.transpose() + Q; 
}

/**
 * <summary> Kalman Filter update. </summary>
 * <param name = "_gc"> A Gaussian component to predict. </param>
 * <param name = "_sensor"> A sensor to take measurement from. </param>
 * <param name = "_zNum"> Measurement number. </param>
 */
void KalmanFilter::update(gaussian_component& _gc, Sensor& _sensor, const size_t& _zNum)
{
	MatrixXd S = MatrixXd::Zero(_sensor.zDim, _sensor.zDim);
	MatrixXd K = MatrixXd::Zero(_sensor.sDim, _sensor.zDim);

	S = _sensor.H * _gc.P * _sensor.H.transpose() + _sensor.R;     // Innovation (residual) covariance
	_sensor.setS(S);

	K = _gc.P * _sensor.H.transpose() * S.inverse();             // Kalman gain

	_gc.m = _gc.m + K * (_sensor.z.at(_zNum) - _sensor.H * _gc.m);						// Updated mean
	_gc.P = (MatrixXd::Identity(_sensor.sDim, _sensor.sDim) - K * _sensor.H) * _gc.P;	// Updated covariance
}

/**
 * <summary> Timestep modifier. </summary>
 * <param name = "_t"> New timestep. </t>
 */
void KalmanFilter::setT(const decltype(dt)& _dt)
{
	dt = _dt;
	F = getCVF((size_t)F.rows(), _dt);
}

/**
 * <summary> Get constant velocity state transition matrix of the required dimension. </summary>
 * <param name = "_dim"> The dimensionality of the state vector. </param>
 * <param name = "_dt"> The timestep. </param>
 * <returns> State transition MatrixXd of the specified dimension. </returns>
 */
MatrixXd KalmanFilter::getCVF(const size_t & _dim, const double& _dt)
{
	assert(_dim % 2 == 0 && "Works only for even-sized matrices");
	MatrixXd F = MatrixXd::Identity(_dim, _dim);

	for (size_t i = 0, j = _dim / 2; j < _dim; i++, j++)
		F(i, j) = _dt;

	return F;
}

/**
* <summary> Get state transition covariance matrix for constant velocity model. </summary>
* <param name = "_dim"> The dimensionality of the state vector. </param>
* <param name = "_dt"> The timestep. </param>
* <returns> Process noise MatrixXd of the specified dimension. </returns>
*/
MatrixXd KalmanFilter::getCVQ(const size_t& _dim, const double& _dt)
{
	assert(_dim % 2 == 0 && "Works only for even-sized matrices");
	VectorXd q(_dim);
	double dt2 = _dt * _dt / 2;
	size_t _dim2 = _dim / 2;

	for (size_t i = 0; i < _dim / 2; i++) {
		q(i) = dt2;
		q(i + _dim2) = _dt * 0.01;
	}

	return q * q.transpose();
}

/**
 * <summary> Get auxialry vector for state transition covariance matrix calculation. </summary>
 * <param name = "_dim"> The dimensionality of the state vector. </param>
 * <param name = "_dt"> The timestep. </param>
 * <returns> Auxilary vector of the specified dimension. </returns>
 */
VectorXd KalmanFilter::getCVq(const size_t & _dim, const double & _dt)
{
	assert(_dim % 2 == 0 && "Works only for even-sized matrices");
	VectorXd q(_dim);
	double dt2 = _dt * _dt / 2;
	size_t _dim2 = _dim / 2;

	for (size_t i = 0; i < _dim / 2; i++) {
		q(i) = dt2;
		q(i + _dim2) = _dt * 0.01;
	}

	return q;
}
