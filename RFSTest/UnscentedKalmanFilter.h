#pragma once
#include "KalmanFilter.h"
#include <vector>

class UnscentedKalmanFilter :
	public KalmanFilter
{

protected:
	double sigmaSamplingW;

public:
	UnscentedKalmanFilter();
	UnscentedKalmanFilter(const decltype(Q) _Q, const decltype(dt) _dt = 0);
	UnscentedKalmanFilter(const decltype(Q) _Q, const decltype(sigmaSamplingW) _w, const decltype(dt) _dt = 0);

	virtual void predict(gaussian_component& _gc);
	virtual void update(gaussian_component& _gc, Sensor& _sensor, const size_t& _zNum);

	static void getSigmaPoints(const VectorXd &_mean, const MatrixXd &_cov, const double &_w0,
		std::vector<Eigen::VectorXd>& _points, std::vector<double>& _weights);
};

