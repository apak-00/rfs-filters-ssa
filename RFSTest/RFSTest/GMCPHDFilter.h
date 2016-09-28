#pragma once
#include <Eigen/Dense>
#include "GMRFSFilter.h"
#include "Sensor.h"
#include "gmm.h"

class GMCPHDFilter :
	public GMRFSFilter <gaussian_mixture>
{
private:
	unsigned int nBirthComponents;	// Number of birth components
	double birthIntensity;			// Birth intensity
	double pS;						// Probability of survival

	MatrixXd initialCovariance;		// Initial component covariance

	VectorXd lowerBound;			// Lower birth bound
	VectorXd upperBound;			// Upper birth bound

	// CPHD
	size_t NMax;					// Maximum expected number of targets	

	VectorXd cardinality;

public:
	GMCPHDFilter(std::shared_ptr<KalmanFilter> _kf, const unsigned int& _nBirthComponents, const double& _birthIntensity,
		const double& _pS, const MatrixXd& _iCov, const VectorXd& _lBound, const VectorXd& _uBound, const size_t& _NMax);

	virtual void predict(gaussian_mixture& _gmm, Sensor& _sensor);
	virtual void update(gaussian_mixture& _gmm, Sensor& _sensor);

	auto getQ() { return -1; }

};

