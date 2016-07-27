#pragma once
#include "UnscentedKalmanFilter.h"
#include "ExtendedKalmanFilter.h"

class TestFilter :
	public UnscentedKalmanFilter, public ExtendedKalmanFilter, public KalmanFilter
{
public:
	TestFilter(); 
	TestFilter(const decltype(ExtendedKalmanFilter::Q) _Q, const decltype(sigmaSamplingW) _w, const decltype(UnscentedKalmanFilter::dt) _dt = 0);

	void predict(gaussian_component& _gc);
	void update(gaussian_component& _gc, Sensor& _sensor, const size_t& _zNum);

	void setT(const decltype(KalmanFilter::dt)& _dt);
};

