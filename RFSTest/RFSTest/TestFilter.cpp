#include "TestFilter.h"
#include "unscented_sampler.hpp"

TestFilter::TestFilter() : KalmanFilter () {}

TestFilter::TestFilter(const decltype(ExtendedKalmanFilter::Q) _Q, const decltype(sigmaSamplingW) _w, 
	const decltype(UnscentedKalmanFilter::dt) _dt) : UnscentedKalmanFilter (_Q, _w, _dt) {}

void TestFilter::predict(gaussian_component & _gc)
{
	UnscentedKalmanFilter::predict(_gc);
}

void TestFilter::update(gaussian_component & _gc, Sensor & _sensor, const size_t & _zNum)
{
	ExtendedKalmanFilter::update(_gc, _sensor, _zNum);
}

void TestFilter::setT(const decltype(KalmanFilter::dt)& _dt)
{
	KalmanFilter::setT(_dt);
}

