#pragma once
#include <Eigen/Dense>
#include "KalmanFilter.h"
#include "Sensor.h"
#include "gmm.h"

/*
* <summary> Gaussian Mixture Probability Hypothesis Density Filter class. </summary>
*/
class GMPHDFilter
{
protected:
	KalmanFilter kf;
	unsigned int nBirthComponents;	// Number of birth components
	double birthIntensity;			// Birth intensity
	double pS;						// Probability of survival

	MatrixXd initialCovariance;		// Initial component covariance

	VectorXd lowerBound;			// Lower birth bound
	VectorXd upperBound;			// Upper birth bound

public:
	GMPHDFilter(const KalmanFilter& _kf, const unsigned int& _nBirthComponents, const double& _birthIntensity,
		const double& _pS, const MatrixXd& _iCov, const VectorXd& _lBound, const VectorXd& _uBound);

	virtual void predict(gaussian_mixture& _gmm);
	virtual void update(gaussian_mixture& _gmm, Sensor& sensor);

	void updateKFTimestep(const double& _t);
};

