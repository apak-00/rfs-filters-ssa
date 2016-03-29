#include "GMJoTTFilter.h"

// Temporary
#include <iostream>
using namespace std;

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
template<typename T>
GMJoTTFilter<T>::GMJoTTFilter(const T& _kf, const unsigned int & _nBirthComponents, const double & _birthIntensity,
	const double & _pS, const MatrixXd & _iCov, const VectorXd & _lBound, const VectorXd & _uBound, const double & _q, const double& _pB) : 
	kf(_kf), nBirthComponents(_nBirthComponents), birthIntensity(_birthIntensity), pS(_pS), initialCovariance(_iCov),
	lowerBound(_lBound), upperBound(_uBound), q(_q), pB(_pB) {}

/**
 * <summary> Prediction step of the GM JoTT filter. </summary>
 * <par> Updates the probability of target existence q, state of the Gaussian Mixture, 
 * stores previous GM components' weights for the update step. </par>
 * 
 * <param name = "_gmm"> A reference to the Gaussian Mixture to be precited. </param>
 */


template<typename T>
void GMJoTTFilter<T>::predict(gaussian_mixture & _gmm)
{
	// Probability of target existence
	double qPred = pB * (1 - q) + pS * q;

	// Previous weights
	previousWeights = _gmm.getWeights();

	// Predict existing components
	for (auto &gc : _gmm.components) {
		kf.predict(gc);
		gc.w *= pS * q / qPred;
	}

	double initialWeight = (birthIntensity / nBirthComponents) * pB * (1 - q) / qPred;
	for (size_t i = 0; (i < nBirthComponents) && (_gmm.size() < _gmm.nMax); i++) {

		//_gmm.addComponent(gaussian_mixture::randVec(lowerBound, upperBound), initialCovariance, initialWeight);

		// Uniform birth test
		VectorXd nbMean(2);
		nbMean << 1000, 0;
		MatrixXd nbCov(2, 2);
		nbCov << 500, 0, 0, 10;

		_gmm.addComponent(gaussian_component(nbMean, nbCov, initialWeight, _gmm.idCounter++));
		previousWeights.push_back(initialWeight);
	}

	q = qPred;
}

/**
 * <summary> Update step of the GM JoTT filter. </summary>
 *
 * <param name = "_gmm"> A reference to the Gaussian Mixture to be updated. </param>
 * <param name = "_sensor"> A reference to the sensor to read the measurements from. </sensor>
 */
template <typename T>
void GMJoTTFilter<T>::update(gaussian_mixture & _gmm, Sensor & _sensor)
{
	double cz = 0.1;// 1.0 / 231609.0;			// Temporary fix for cz

	// Multiply everything by (1 - pD) (First term of the JoTT update)
	auto pD = _sensor.getPD();
	size_t n0 = _gmm.size();
	double delta_k = 0;

	// [2] Second term first, avoiding temporary variables
	for (size_t i = 0; i < _sensor.getZ().size(); i++) {
	
		size_t n = _gmm.size();

		for (size_t j = 0; j < n0; j++) {

			gaussian_component gct(_gmm[j]);
			kf.update(gct, _sensor, i);

			//if (gct.tag[0] == )

			auto qk = (1 / sqrt(pow(2 * M_PI, _sensor.getZDim()) * _sensor.getS().determinant())) * exp(-0.5 * _sensor.zMahalanobis(gct.m, i));
			gct.w *= qk / (_sensor.getLambda() * cz);		

			//delta_k += previousWeights[j] * qk / (_sensor.getLambda() * cz);
			delta_k += gct.w;

			//if (gct.w > 1e-7)
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

/**
 * <summary> Updates the timestep of the Kalman Filter. </summary>
 */
template <typename T>
void GMJoTTFilter<T>::updateKFTimestep(const double & _t)
{
	kf.setT(_t);
}
