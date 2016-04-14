#pragma once
#include <Eigen/Dense>
#include "Sensor.h"
#include "gmm.h"

using namespace Eigen;

/**
 * <summary> Kalman Filter class. </summary>
 */
class KalmanFilter
{

protected:
	MatrixXd F;		// State transition matrix
	MatrixXd Q;		// Process noise
	double dt;		// Timestep

public:
	KalmanFilter();
	KalmanFilter(const decltype(F)& _F, const decltype(Q) _Q, const decltype(dt) _dt);
	KalmanFilter(const KalmanFilter& _kf);

	KalmanFilter& operator= (const KalmanFilter& _kf);

	virtual void predict(gaussian_component& _gc);
	virtual void update(gaussian_component& _gc, Sensor& _sensor, const size_t& _zNum);

	void setT(const decltype(dt)& _dt);

	auto getTimestep() { return dt; };

	static MatrixXd getCVF(const size_t& _dim, const double& _dt);
	static MatrixXd getCVQ(const size_t& _dim, const double& _dt);

};

