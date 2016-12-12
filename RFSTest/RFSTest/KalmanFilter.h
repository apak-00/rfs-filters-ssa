#pragma once
#include <Eigen/Dense>
#include "Sensor.h"
#include "MixtureModels.h"
#include "SingleTargetFilter.h"

using namespace Eigen;

/**
 * <summary> Kalman Filter class. </summary>
 */
class KalmanFilter : public SingleTargetFilter<gaussian_component>
{

protected:
	MatrixXd F;		// State transition matrix
	MatrixXd Q;		// Process noise

public:
	KalmanFilter();
	KalmanFilter(const decltype(F)& _F, const decltype(Q) _Q, const decltype(dt) _dt = 0);
	KalmanFilter(const KalmanFilter& _kf);

	KalmanFilter& operator= (const KalmanFilter& _kf);

	virtual void predict(gaussian_component& _gc);
	virtual void update(gaussian_component& _gc, Sensor& _sensor, const size_t& _zNum = 0);

	virtual void setT(const decltype(dt)& _dt);

	static MatrixXd getCVF(const size_t& _dim, const double& _dt);
	static MatrixXd getCVQ(const size_t& _dim, const double& _dt);
	static VectorXd getCVq(const size_t& _dim, const double& _dt);

	bool debug = false;
};

