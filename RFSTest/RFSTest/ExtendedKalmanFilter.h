#pragma once
#include "KalmanFilter.h"

class ExtendedKalmanFilter:
	public KalmanFilter
{
public:
	ExtendedKalmanFilter();
	ExtendedKalmanFilter(const decltype(F)& _F, const decltype(Q) _Q, const decltype(dt) _dt);
	ExtendedKalmanFilter(const ExtendedKalmanFilter& _ekf);

	void predict(gaussian_component& _gc);
    void update(gaussian_component& _gc, Sensor& _sensor, const size_t& _zNum);
};

