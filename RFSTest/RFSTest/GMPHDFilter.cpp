// TODO: Complete

#include "GMPHDFilter.h"
#include "MathHelpers.h"

/**
 * <summary> Standard constructor of the GMPHDFilter. </summary>
 *
 * <param name = "_kf"> An instance of the Kalman Filter for single target state propagation. </param>
 * <param name = "_nBirthComponents"> Number of birth components for the Gaussian Mixture during the prediction step. </param>
 * <param name = "_birthIntensity"> Intensity of the birth components. </param>
 * <param name = "_pS"> Probability of target survival. </param>
 * <param name = "_iCov"> Initial target state uncertainty (covariance). </param>
 * <param name = "_lBound"> Lower state bound for random state generation. </param>
 * <param name = "_uBound"> Upper state bound for random state generation. </param>
 */
GMPHDFilter::GMPHDFilter(std::shared_ptr<KalmanFilter> _kf, const unsigned int & _nBirthComponents, const double & _birthIntensity,
	const double & _pS, const MatrixXd& _iCov, const VectorXd & _lBound, const VectorXd & _uBound) : nBirthComponents(_nBirthComponents),
	birthIntensity(_birthIntensity), pS(_pS), initialCovariance(_iCov), lowerBound(_lBound), upperBound(_uBound) 
	{
		filter = _kf;
	}

/**
 * <summary> Prediction step of the GM PHD filter. </summary>
 * <param name = "_gmm"> A reference to the Gaussian Mixture to be precited. </param>
 */
void GMPHDFilter::predict(gaussian_mixture & _gmm, Sensor& _sensor)
{
	double range;
	std::vector<double> birthRanges;
	VectorXd birth(_gmm.dim());

	// Prediction for all of the components
	for (auto &gc : _gmm.components) {
		filter->predict(gc);
		gc.w *= pS;
	}

	// New target birth
	double initialWeight = birthIntensity / nBirthComponents;

	if (_sensor.getZ().size() != 0)
	{
		for (size_t i = 0; i < nBirthComponents; i++)
			birthRanges.push_back((double)(i + 2) * 200);

		if (nBirthComponents == 1)
		{
			birthRanges[0] = 1000; // + rand() % 500 - 250;
		}

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
	
		/*
	for (size_t i = 0; (i < nBirthComponents) && (_gmm.size() < _gmm.nMax); i++)
		// Changed 22/4/2016
		_gmm.addComponent(gaussian_component(gaussian_mixture::randVec(lowerBound, upperBound), initialCovariance, initialWeight, _gmm.idCounter++));
		*/
}

	/**
 * <summary> Update step of the GM PHD Filter. </summary>
 * <param name = "_sensor"> A reference to the sensor to read the measurements from. </sensor>
 */
void GMPHDFilter::update(gaussian_mixture& _gmm, Sensor & _sensor)
{
	auto pD = _sensor.getPD();
	size_t n0 = _gmm.size();

	// [2] Compute second term first, avoiding temporary variables
	// Per measurement
	for (size_t i = 0; i < _sensor.z.size(); i++) {

		double weightSum = 0.0;
		size_t n = _gmm.size();
		
		// Perform updates to the predicted set
		for (size_t j = 0; j < n0; j++) {

			gaussian_component gct(_gmm[j]);
			filter->update(gct, _sensor, i);

			// Gaussian likelihood
			auto q = (1.0 / sqrt(pow(2.0 * M_PI, _sensor.getZDim()) * _sensor.getS().determinant()))
				* exp(-0.5 * MathHelpers::mahalanobis(_sensor.getZ(i), _sensor.getPredictedZ(), _sensor.getS()));
			gct.w *= pD * q;
			weightSum += gct.w;

			//if (gct.w > 1e-7) 
				_gmm.addComponent(gct);
		}

		// Weight normalization
		for (size_t j = n; j < _gmm.size(); j++) 
			_gmm[j].w /= weightSum + _sensor.kappa;
	}

	// [1] Compute first term (Predicted mixture x (1 - pD))
	for (size_t i = 0; i < n0; i++)
		_gmm[i].w *= (1 - pD);

	// If the track got split, the component with the highest weight keeps the track
	for (size_t i = 0; i < _gmm.components.size(); i++)
		for (size_t j = 0; j < _gmm.components.size(); j++)
			if (_gmm.components[i].tag[0] == _gmm.components[j].tag[0] && i != j)
				if (_gmm.components[i].w > _gmm.components[j].w)
					_gmm.components[j].initTag(_gmm.idCounter++);
				else
					_gmm.components[i].initTag(_gmm.idCounter++);

	/*
	for (auto &gci : _gmm.components)
		for (auto &gcj : _gmm.components)
			if (gci.tag[0] == gcj.tag[0] && gci.w != gcj.w) {
				//(gci.w > gcj.w) ? (gcj.initTag(_gmm.componentCounter++)) : (gci.initTag(_gmm.componentCounter++));
				if (gci.w > gcj.w)
					gcj.initTag(_gmm.componentCounter++);
				else
					gci.initTag(_gmm.componentCounter++);
			}
	*/
}
