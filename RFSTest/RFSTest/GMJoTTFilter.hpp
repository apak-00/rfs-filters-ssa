#pragma once
#include <Eigen/Dense>
#include "GMRFSFilter.h"
#include "Sensor.h"
#include "gmm.h"
#include "MathHelpers.h"

/*
* <summary> Gaussian Mixture Joint Target Detection and Tracking Filter class. </summary>
*/
class GMJoTTFilter : public GMRFSFilter<gaussian_mixture>
{
protected:
	double q;						// Probability of target existence

	size_t nBirthComponents;	    // Number of birth components
	double birthIntensity;			// Birth intensity
	double pS;						// Probability of target survival
	double pB;						// Probabilirt of target birth

	MatrixXd initialCovariance;		// Initial component covariance

	VectorXd lowerBound;			// Lower birth bound
	VectorXd upperBound;			// Upper birth bound

public:
	auto getQ() { return q; }
	auto getT() { return filter->getT(); }
	
	/**
	* <summary> Main constructor of the JoTT filter class. </summary>
	*
	* <param name = "_kf"> An instance of the Kalman Filter for single target state propagation. </param>
	* <param name = "_nBirthComponents"> Number of birth components for the Gaussian Mixture during the prediction step. </param>
	* <param name = "_birthIntensity"> Intensity of the birth components. </param>
	* <param name = "_pS"> Probability of target survival. </param>
	* <param name = "_iCov"> Initial target state uncertainty (covariance). </param>
	* <param name = "_lBound"> Lower state bound for random state generation. </param>
	* <param name = "_uBound"> Upper state bound for random state generation. </param>
	* <param name = "_q"> Initial probability of target existence. </param>
	* <param name = "_pB"> Probability of target birth. </param>
	*/
	GMJoTTFilter(std::shared_ptr<KalmanFilter> _kf, const size_t& _nBirthComponents, const double & _birthIntensity,
		const double & _pS, const MatrixXd & _iCov, const VectorXd & _lBound, const VectorXd & _uBound, const double & _q, const double& _pB) :
		nBirthComponents(_nBirthComponents), birthIntensity(_birthIntensity), pS(_pS), initialCovariance(_iCov),
		lowerBound(_lBound), upperBound(_uBound), q(_q), pB(_pB) 
	{
		filter = _kf;
	}

	/**
	* <summary> Prediction step of the GM JoTT filter. </summary>
	* <par> Updates the probability of target existence q, state of the Gaussian Mixture,
	* stores previous GM components' weights for the update step. </par>
	*
	* <param name = "_gmm"> A reference to the Gaussian Mixture to be precited. </param>
	*/
	void predict(gaussian_mixture & _gmm, Sensor& _sensor)
	{
		// Probability of target existence
		double qPred = pB * (1 - q) + pS * q, range;
		double initialWeight = (birthIntensity / nBirthComponents) * pB * (1 - q) / qPred;
		std::vector<double> birthRanges;
		VectorXd birth(_gmm.dim());

		// Predict existing components
		for (auto &gc : _gmm.components) {
			filter->predict(gc);
			gc.w *= pS * q / qPred;
		}

		q = qPred;

		if (_sensor.getZ().size() != 0)
		{
			for (size_t i = 0; i < nBirthComponents; i++)
				birthRanges.push_back((double)(i + 2) * 200);

			if (nBirthComponents == 1)
			{
				birthRanges[0] = 1000 + rand() % 500 - 250;
			}

			// Birth
			for (size_t i = 0; (i < nBirthComponents) && (_gmm.size() < _gmm.nMax); i++)
			{
				// Uniform birth test
				range = birthRanges[i];

				if (_gmm.dim() == 2)
					birth << range, 0;
				else if (_gmm.dim() == 6)
				{
					VectorXd m(_gmm.dim());
					m << range, _sensor.getBearing(), 0, 0, 0;
					birth = Astro::razelToTEME(m, _sensor.getPosition(), _sensor.getDateJD(), _sensor.getLOD(), _sensor.getXp(), _sensor.getYp());
				}
				_gmm.addComponent(gaussian_component(birth, initialCovariance, initialWeight, _gmm.idCounter++));
			}
		}

	}

	/**
	* <summary> Update step of the GM JoTT filter. </summary>
	*
	* <param name = "_gmm"> A reference to the Gaussian Mixture to be updated. </param>
	* <param name = "_sensor"> A reference to the sensor to read the measurements from. </sensor>
	*/
	void update(gaussian_mixture & _gmm, Sensor & _sensor)
	{
		double cz = 1.0 / 57903 * 200; // 200 UKF 42 EKF

		auto pD = _sensor.getPD();
		size_t n0 = _gmm.size();
		double delta_k = 0; 

		// [2] Second term first, avoiding temporary variables
		for (size_t i = 0; i < _sensor.getZ().size(); i++) {

			size_t n = _gmm.size();

			for (size_t j = 0; j < n0; j++) {

				gaussian_component gct(_gmm[j]);
				filter->update(gct, _sensor, i);
				
				auto qk = (1.0 / sqrt(pow(2.0 * M_PI, _sensor.getZDim()) * _sensor.getS().determinant())) 
					* exp(-0.5 * MathHelpers::mahalanobis(_sensor.getZ(i), _sensor.getPredictedZ(), _sensor.getS()));
				gct.w *= qk / (_sensor.getLambda() * cz);

				delta_k += gct.w;

				//if (gct.w > 1e-7)vt
				_gmm.addComponent(gct);
			}
		}

		delta_k = pD * (1 - delta_k);

		// [1] First term weight update (0 -> n0)
		// (Predicted mixture x (1 - pD) / (1 - delta_k))
		for (size_t j = 0; j < n0; j++)
			_gmm[j].w *= (1 - pD) / (1 - delta_k);

		// [2] Second term weight update (n0 -> end)
		for (size_t j = n0; j < _gmm.size(); j++)
			_gmm[j].w *= pD / (1 - delta_k);

		q = (1 - delta_k) / (1 - q * delta_k) * q;

		// If the track got split, the component with the highest weight keeps the track
		for (size_t i = 0; i < _gmm.components.size(); i++)
			for (size_t j = 0; j < _gmm.components.size(); j++)
				if (_gmm.components[i].tag[0] == _gmm.components[j].tag[0] && i != j)
					if (_gmm.components[i].w > _gmm.components[j].w)
						_gmm.components[j].initTag(_gmm.idCounter++);
					else
						_gmm.components[i].initTag(_gmm.idCounter++);
	}
};

