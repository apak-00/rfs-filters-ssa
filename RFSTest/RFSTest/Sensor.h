#pragma once
#ifndef SENSOR_H
#define SENSOR_H

#include <vector>
#include <Eigen/Eigen>
#include "Astro.h"

using namespace Eigen;

/**
* <summary> Sensor class emulates a sensor profile used to obtain measurements of the orbiting objects. </summary>
* <para> Initially based on D. Clark's code SensorProfile class. </para>
* <date> 2016-01-12 </date>
* <author> Andrey Pak </author>
*/
class Sensor
{
	friend class KalmanFilter;  // For ease of access during update step
	friend class ExtendedKalmanFilter;
	friend class UnscentedKalmanFilter;
	friend class GMPHDFilter;
	friend class GMCPHDFilter;

protected:
	size_t zDim;  // Observation dimension
	size_t sDim;  // State dimension
	double pD;          // Probability of detection
	double lambda;      // Expected number of false alarms, Poisson
	double V;           // ???
	double kappa;       // Lambda * clutter distribution, Poisson

	std::vector<VectorXd> z;		// Measurements' vector
	std::vector<VectorXd> zPrev;	// Previous measurements
	MatrixXd R;         // Observation noise
	MatrixXd H;         // Observation matrix (state space -> observation space)
	MatrixXd S;

	// TODO: Remove in the future
	VectorXd predictedZ;

	VectorXd position;  // Sensor positionin WGS-84
	VectorXd bearing;   // Sensor bearing and bearing rates

	Astro::date zDate;  // Measurement Date

	double t;           // Timestep
	double signalStrength;              // Radar-specific: returned signal strength
	double signalStrengthThreshold;     // Radar-specific: the lowest signal strength
										// when the proper measurement can be obtained

	double xp;			// Polar motion coefficient x
	double yp;			// Polar motion coefficient y
	double lod;			// Length-of-day (see Vallado)

public:
	Sensor();
	Sensor(const size_t& _zDim, const size_t& _sDim, const double& _pD, const double& _lambda, const double& _V);
	Sensor(const size_t& _zDim, const size_t& _sDim, const double& _pD, const double& _lambda, const double& _V,
		const MatrixXd& _R, const MatrixXd& _H, const VectorXd& _position = VectorXd::Zero(3),
		const double& signalStrengthThreshold = 0);
	Sensor(const Sensor& _sensor);
	Sensor& operator=(const Sensor& _sensor);

	double getDateJD() const;

	/* Accessors and mutators */
	const auto getPD() const { return pD; }
	const auto getLambda() const { return lambda; }
	const auto getKappa() const { return kappa; }
	const auto getZDim() const { return zDim; }
	const auto getSDim() const { return sDim; }
	const auto getZ() const { return z; }
	const auto getZPrev() const { return zPrev; }
	const auto getZ(const size_t& _id) { return z[_id]; }
	const auto getR() const { return R; }
	const auto getH() const { return H; }
	const auto getS() const { return S; }
	const auto getPosition() const { return position; }
	const auto getBearing() const { return bearing; }
	const auto getT() const { return t; }
	const auto getSignalStrength() const { return signalStrength; }
	const auto getSignalStrengthThreshold() const { return signalStrengthThreshold; }
	const auto getDate() const { return zDate; }
	const auto getXp() const { return xp; }
	const auto getYp() const { return yp; }
	const auto getLOD() const { return lod; }

	// TODO: Remove in the future
	const auto getPredictedZ() const { return predictedZ; }
 
	void setPD(const decltype(pD)& _pD) { pD = _pD; }
	void setLambda(const decltype(lambda) & _lambda) { lambda = _lambda; }
	void setKappa(const decltype(kappa) & _kappa) { kappa = _kappa; }
	void setZDim(const decltype(zDim) & _zDim) { zDim = _zDim; }
	void setSDim(const decltype(sDim) & _sDim) { sDim = _sDim; }
	void setZ(const decltype(z) & _z) { z = _z; }
	void setZPrev(const decltype(zPrev) _zPrev) { zPrev = _zPrev; };
	void setR(const decltype(R) & _R) { R = _R; }
	void setH(const decltype(H) & _H) { H = _H; }
	void setS(const decltype(S) & _S) { S = _S; }
	void setPosition(const decltype(position) & _position) { position = _position; }
	void setBearing(const decltype(bearing) & _bearing) { bearing = _bearing; }
	void setT(const decltype(t) & _t) { t = _t; }
	void setSignalStrength(const decltype (signalStrength) & _ss) { signalStrength = _ss; }
	void setSignalStrengthThreshold(const decltype (signalStrengthThreshold) & _sst) { signalStrengthThreshold = _sst; }
	void setDate(const decltype (zDate) & _zDate) { zDate = _zDate; }
	void setXp(const decltype (xp)& _xp) { xp = _xp; }
	void setYp(const decltype (yp)& _yp) { yp = _yp; }
	void setLOD(const decltype(lod)& _lod) { lod = _lod; }

	// TODO: Remove in the future
	void setPredictedZ(const decltype(predictedZ)& _pz) { predictedZ = _pz; }
};

#endif // SENSOR_H
