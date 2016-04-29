#pragma once
#include "KalmanFilter.h"

class UnscentedKalmanFilter :
	public KalmanFilter
{
public:

	UnscentedKalmanFilter();

	void predict();
	void update();

};

