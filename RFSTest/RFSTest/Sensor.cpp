#include <math.h>
#include "Sensor.h"

#include <iostream>

Sensor::Sensor() : zDim(0), sDim(0), pD(0), lambda(0), V(0), kappa(0), t(0),
	signalStrength (0), signalStrengthThreshold (0), zDate(), xp(0), yp(0), lod(0) {}

/**
* <summary> Constructor of the Sensor class (I). </summary>
* <para> Takes a number of basic parameters and initializes everything else to zero. </para>
* <para> Measurement covariance (R), observation matrx (H), position are to be set afterwards. </para>
* <para> kappa is initialized as lambda * V / sqrt(pow(2 * M_PI, _zDim)).
*
* <param name = "_zDim"> Dimensionality of the observation. </param>
* <param name = "_sDim"> Dimensionality of the state. </param>
* <param name = "_pD"> Probability of detection. </param>
* <param name = "_lambda"> Expected numer of false alarms. </param>
* <param name = "_V"> Clutter distribution. </param>
*/
Sensor::Sensor(const size_t& _zDim, const size_t& _sDim, const double& _pD, const double& _lambda, const double& _V) :
	zDim(_zDim), sDim(_sDim), pD(_pD), lambda(_lambda), V(_V), kappa(lambda * V / sqrt(pow(2 * M_PI, _zDim))),
	z(), R(MatrixXd::Zero(_zDim, _zDim)), H(MatrixXd::Zero(_zDim, _sDim)), S(MatrixXd::Zero(_zDim, _zDim)),
	position(VectorXd::Zero(3)), bearing(VectorXd::Zero(4)), t(0), signalStrength(0), signalStrengthThreshold(0), 
	zDate() ,xp(0), yp(0), lod(0) {}

/**
* <summary> Constructor of the Sensor class (I). </summary>
* <para> Takes a number of parameters and initializes all variable (e.g. measurement) to zero. </para>
* <para> kappa is initialized as lambda * V / sqrt(pow(2 * M_PI, _zDim)).
*
* <param name = "_zDim"> Dimensionality of the observation. </param>
* <param name = "_sDim"> Dimensionality of the state. </param>
* <param name = "_pD"> Probability of detection. </param>
* <param name = "_lambda"> Expected numer of false alarms. </param>
* <param name = "_V"> Clutter distribution. </param>
* <param name = "_R"> Observation noise matrix. </param>
* <param name = "_H"> Observation matrix </param>
* <param name = "_position"> Position of the senor in WGS-84 coordinates
* (latitude (deg), longitude (deg), altitude (m)) </param>
* <param name = "_sst"> Minimal signal strength threshold value. </param>
*/
Sensor::Sensor(const size_t& _zDim, const size_t& _sDim, const double& _pD, const double& _lambda, const double& _V,
	const MatrixXd& _R, const MatrixXd& _H, const VectorXd& _position, const double& _sst) :
	zDim(_zDim), sDim(_sDim), pD(_pD), lambda(_lambda), V(_V), kappa(lambda * V / sqrt(pow(2 * M_PI, _zDim))),
	z(), R(_R), H(_H), S(MatrixXd::Zero(_zDim, _zDim)), position(_position),
	bearing(VectorXd::Zero(4)), t(0), signalStrength(0), signalStrengthThreshold(_sst), 
	zDate(), xp(0), yp(0), lod(0) {}

/**
* <summary> Copy constructor of the Sensor class. </summary>
* <param name = "_s"> A constant reference to the instance of the Sensor class to copy from. </param>
*/
Sensor::Sensor(const Sensor& _s) : zDim(_s.zDim), sDim(_s.sDim), pD(_s.pD), lambda(_s.lambda), V(_s.V),
z(_s.z), R(_s.R), H(_s.H), S(_s.S), position(_s.position), bearing(_s.bearing), t(_s.t),
signalStrength(_s.signalStrength), signalStrengthThreshold(_s.signalStrengthThreshold), zDate(_s.zDate), xp(_s.xp), yp(_s.yp), lod(_s.lod) {}

/**
* <summary> Assignment operator overload. </summary>
* <param name = "_s"> A constant referece to the assigned Sensor instance. </param>
*/
Sensor& Sensor::operator = (const Sensor& _s)
{
	pD = _s.pD;
	lambda = _s.lambda;
	kappa = _s.kappa;
	zDim = _s.zDim;
	sDim = _s.sDim;
	z = _s.z;
	R = _s.R;
	H = _s.H;
	S = _s.S;
	position = _s.position;
	bearing = _s.bearing;
	t = _s.t;
	signalStrength = _s.signalStrength;
	signalStrengthThreshold = _s.signalStrengthThreshold;
	zDate = _s.zDate;
	xp = _s.xp;
	yp = _s.yp;
	lod = _s.lod;

	return *this;
}

/**
 * <summary> Returns the Julian Day number for the current sensor date. </summary>
 * <returns> The Julian Day number. </returns>
 */
double Sensor::getDateJD() const
{
	return Astro::getJulianDay(zDate);
}
