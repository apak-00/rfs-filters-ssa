#pragma once
#include <Eigen/Dense>
#include "Sensor.h"
#include "gmm.h"

/**
 * <summary> Template class for single-target filter. </summary>
 */
template <typename C>
class SingleTargetFilter
{
protected:
	double dt;

public:
	SingleTargetFilter() { dt = 0; }
	SingleTargetFilter(const double& _dt) : dt(_dt) {}
	~SingleTargetFilter() {}

	virtual void predict(C& _c) = 0;
	virtual void update(C& _c, Sensor& _sensor, const size_t& _zNum = 0) = 0;

	// Timestep
	virtual void setT(const decltype(dt)& _dt) { dt = _dt; }
	virtual double getT() { return dt; };
};

