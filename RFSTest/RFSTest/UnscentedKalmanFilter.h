#pragma once
#include "KalmanFilter.h"
#include <vector>

class UnscentedKalmanFilter :
	public KalmanFilter
{
private:

public:

	UnscentedKalmanFilter();
	UnscentedKalmanFilter(const decltype(Q) _Q, const decltype(dt) _dt = 0);

	void predict(gaussian_component& _gc);
	void update(gaussian_component& _gc, Sensor& _sensor, const size_t& _zNum);

};

