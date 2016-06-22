#pragma once
#include <memory>
#include "KalmanFilter.h"
 
/**
 * <summary> Gaussian Mixture Random Finite Set Filter class. </summary>
 */
template <typename Mixture>
class GMRFSFilter
{
protected:
	std::shared_ptr<KalmanFilter> filter;

public:
	GMRFSFilter() {};
	~GMRFSFilter() {};

	virtual void predict(Mixture& _m, Sensor& _sensor) = 0;
	virtual void update(Mixture& _m, Sensor& _sensor) = 0;

	auto getT() { return filter->getT(); }
	virtual void setT(const double & _t) { filter->setT(_t); }
};
