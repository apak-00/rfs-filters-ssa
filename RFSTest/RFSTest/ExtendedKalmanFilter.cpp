#include "ExtendedKalmanFilter.h"
#include <iostream>

/**
 * <summary> Empty constructor of th ExtendedKalmanFilter class. </summary>
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
	
	if (_gc.kindaConverged)
	{
		VectorXd temp(_gc.m.size());
		MatrixXd shepperd = Astro::getShepperdMatrix(_gc.m, dt, temp, Astro::MU_E);

		_gc.m = temp;
		//_gc.P = shepperd * _gc.P * shepperd.transpose() + Q;
		_gc.P = F * _gc.P * F.transpose() + Q;
	}
	else 
	{
		//_gc.m = F * _gc.m;		// Constant Velocity Mean Update
		_gc.m = Astro::integrationPrediction(_gc.m, dt);
		_gc.P = F * _gc.P * F.transpose() + Q;
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

	//MatrixXd H = Astro::getSEZToRAZELJacobian(mSEZ, _sensor.getZDim());			// Hardcoded Jacobian
	MatrixXd H = Astro::getSEZToRAZELJacobianFADBAD(mSEZ, _sensor.getZDim());
	
	MatrixXd S = H * pSEZ * H.transpose() + _sensor.R;
	MatrixXd K = pSEZ * H.transpose() * S.inverse();

	mSEZ += K * (_sensor.z[_zNum] - Astro::sezToRAZEL(mSEZ).head(3));
	pSEZ = (MatrixXd::Identity(_sensor.sDim, _sensor.sDim) - K * H) * pSEZ;

	_gc.m = Astro::sezToTEME(mSEZ, _sensor.getPosition(), _sensor.getDateJD(), _sensor.getLOD(), _sensor.getXp(), _sensor.getYp());
	_gc.P = tf * pSEZ * tf.transpose();

	_sensor.setH(H);
	_sensor.setS(S);
}