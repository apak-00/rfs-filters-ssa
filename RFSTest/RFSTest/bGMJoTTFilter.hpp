#pragma once
#include <Eigen/Dense>
#include "GMRFSFilter.h"
#include "Sensor.h"
#include "MixtureModels.h"

/*
* <summary> Beta-Gaussian Mixture Joint Target Detection and Tracking class. </summary>
*/
class BGMJoTTFilter : public GMRFSFilter<beta_gaussian_mixture>
{
protected:
	double q;						// Probability of target existence

	size_t nBirthComponents;	// Number of birth components
	double birthIntensity;			// Birth intensity
	double pS;						// Probability of target survival
	double pB;						// Probabilirt of target birth

	MatrixXd initialCovariance;		// Initial component covariance

	VectorXd lowerBound;			// Lower birth bound
	VectorXd upperBound;			// Upper birth bound

	double epsilon;					// Epsilon for beta component update

public:
	auto getQ() { return q; }

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
	BGMJoTTFilter(std::shared_ptr<KalmanFilter> _kf, const size_t & _nBirthComponents, const double & _birthIntensity,
		const double & _pS, const MatrixXd & _iCov, const VectorXd & _lBound, const VectorXd & _uBound, const double & _q, const double& _pB, 
		const double& _epsilon) : nBirthComponents(_nBirthComponents), birthIntensity(_birthIntensity), pS(_pS), initialCovariance(_iCov),
		lowerBound(_lBound), upperBound(_uBound), q(_q), pB(_pB), epsilon(_epsilon) 
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
	void predict(beta_gaussian_mixture & _bgmm, Sensor& _sensor)
	{
		// Probability of target existence
		double qPred = pB * (1 - q) + pS * q, range;
		double initialWeight = (birthIntensity / nBirthComponents) * pB * (1 - q) / qPred;
		std::vector<double> birthRanges;
		VectorXd birth(_bgmm.dim());

		// Predict existing components
		for (auto &bgc : _bgmm.components) 
		{
			filter->predict(bgc);
			bgc.w *= pS * q / qPred;
		}

		for (size_t i = 0; i < nBirthComponents; i++)
			birthRanges.push_back((double)(i + 2) * 200);

		if (nBirthComponents == 1)
			birthRanges[0] = 1000 + rand() % 500 - 250;

		// Birth
		for (size_t i = 0; (i < nBirthComponents); i++)
		{
			// Uniform birth test
			range = birthRanges[i];
			std::cout << " B: " << range << " '";

			if (_bgmm.dim() == 2)
				birth << range, 0;
			else if (_bgmm.dim() == 6)
			{
				VectorXd m(_bgmm.dim());
				m << range, _sensor.getBearing(), 0, 0, 0;
				birth = Astro::razelToTEME(m, _sensor.getPosition(), _sensor.getDateJD(), _sensor.getLOD(), _sensor.getXp(), _sensor.getYp());
			}

			_bgmm.addComponent(beta_gaussian_component(birth, initialCovariance, initialWeight, _bgmm.idCounter++, 10, 90));
		}

		q = qPred;

		// Update the beta components
		for (auto &bgc : _bgmm.components)
			updateBetaComponent(bgc, filter->getT() * 0.01);
	}

	/**
	 * <summary> Updates beta distribution with the specified epsilon for the variance increase. </summary>
	 * <param name = "_bgc"> A beta-Gaussian component containing the distribution to update. </param>
	 * <param name = "_epsilon"> Amount of variance increase during the update (between 0 and 1). </param>
	 */
	void updateBetaComponent(beta_gaussian_component& _bgc, const double& _epsilon) 
	{
		double uvrec = 1.0 / (_bgc.u + _bgc.v);
		double variance = (uvrec + _epsilon) * _bgc.u * _bgc.v * uvrec / (_bgc.u + _bgc.v + 1);
		double theta = uvrec * (uvrec * uvrec * _bgc.u * _bgc.v / variance - 1);

		_bgc.u *= theta;
		_bgc.v *= theta;

		if (_bgc.u > 10e10) _bgc.u = 10e4;
		if (_bgc.u < 10e-10) _bgc.u = 10e-4;
		if (_bgc.v > 10e10) _bgc.v = 10e4;
		if (_bgc.v < 10e-10) _bgc.u = 10e-4;
	}

	/**
	* <summary> Update step of the GM JoTT filter. </summary>
	*
	* <param name = "_gmm"> A reference to the Gaussian Mixture to be updated. </param>
	* <param name = "_sensor"> A reference to the sensor to read the measurements from. </sensor>
	*/
	void update(beta_gaussian_mixture & _bgmm, Sensor & _sensor)
	{
		double cz = 1.0 / 57903 * 400;

		size_t n0 = _bgmm.size();
		double delta_k = 0, pDWeightedSum = 0;
		std::vector<double> pD;		// Vector instead of single value for GMJoTT
		pD.resize(n0);

		for (size_t i = 0; i < n0; i++) 
		{
			pD[i] = _bgmm[i].u / (_bgmm[i].u + _bgmm[i].v);
			// Sum for the expected probability detection
			pDWeightedSum += pD[i] * _bgmm[i].w;
		}

		// [2] Second term first, avoiding temporary variables
		for (size_t i = 0; i < _sensor.getZ().size(); i++) 
		{
			size_t n = _bgmm.size();

			for (size_t j = 0; j < n0; j++) {

				beta_gaussian_component bgct(_bgmm[j]);
				filter->update(bgct, _sensor, i);

				
				auto qk = (1.0 / sqrt(pow(2.0 * M_PI, _sensor.getZDim()) * _sensor.getS().determinant()))
					* exp(-0.5 * MathHelpers::mahalanobis(_sensor.getZ(i), _sensor.getPredictedZ(), _sensor.getS()));

				// Update weight
				bgct.w *= pD[j] * qk / (_sensor.getLambda() * cz);

				// New variable pD
				delta_k +=  bgct.w;

				//if (bgct.w > 1e-7)vt
				_bgmm.addComponent(bgct);
			}
		}
		
		// New delta_k for the beta-Gaussian Mixture
		delta_k = pDWeightedSum - delta_k;
		
		// [1] First term weight update (0 -> n0)
		// (Predicted mixture x (1 - pD) / (1 - delta_k))
		for (size_t j = 0; j < n0; j++)
		{
			_bgmm[j].w *= (1 - pD[j]) / (1 - delta_k);
			_bgmm[j].v++;
		}
			
		// [2] Second term weight update (n0 -> end)
		for (size_t j = n0; j < _bgmm.size(); j++)
		{
			_bgmm[j].w *= 1 / (1 - delta_k);
			_bgmm[j].u++;
		}
			
		q = (1 - delta_k) / (1 - q * delta_k) * q;

		// If the track got split, the component with the highest weight keeps the track
		for (size_t i = 0; i < _bgmm.components.size(); i++)
			for (size_t j = 0; j < _bgmm.components.size(); j++)
				if (_bgmm.components[i].tag[0] == _bgmm.components[j].tag[0] && i != j)
					if (_bgmm.components[i].w > _bgmm.components[j].w)
						_bgmm.components[j].initTag(_bgmm.idCounter++);
					else
						_bgmm.components[i].initTag(_bgmm.idCounter++);
	}

};
