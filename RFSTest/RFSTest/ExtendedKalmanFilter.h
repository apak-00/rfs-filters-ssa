#pragma once
#include "KalmanFilter.h"

/*
 * <summary> Extended Kalman Filter class. </summary>
 */
class ExtendedKalmanFilter:
	public KalmanFilter
{
public:
	ExtendedKalmanFilter();
	ExtendedKalmanFilter(const decltype(F)& _F, const decltype(Q) _Q, const decltype(dt) _dt);
	ExtendedKalmanFilter(const ExtendedKalmanFilter& _ekf);

	void predict(gaussian_component& _gc);
    void update(gaussian_component& _gc, Sensor& _sensor, const size_t& _zNum);
private:
	void getSigmaPoints(const VectorXd &_mean, const MatrixXd &_cov, const double &_w0,
		std::vector<Eigen::VectorXd>& _points, std::vector<double>& _weights);
};

